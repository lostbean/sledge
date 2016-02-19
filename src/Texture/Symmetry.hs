{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE TypeFamilies               #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

-- |
-- Module      : Texture.Symmetry
-- Copyright   : (c) 2013 Edgar Gomes
-- License     : Privete-style (see LICENSE)
--
-- Maintainer  : Edgar Gomes <talktoedgar@gmail.com>
-- Stability   : experimental
-- Portability : tested on GHC only
--
-- Module to calculate symmetry in crystal orientations.
--
module Texture.Symmetry
       ( Symm (..)
       , toFZ
       , toFZGeneric
       , toFZDirect
       , getInFZ
       , isInFZ
       , getMinDistFZPlanes
          -- * Vector operations
       , getAllSymmVec
       , getUniqueSymmVecs
         -- * Misorientation
       , getMisoAngle
         -- * Fundamental zone in Frank-Rodrigues
       , isInRodriFZ
        -- * Custom symmetries
       , getSymmOps
       , SymmOp
       , symmOp
       , mkSymmOps
       , SymmAxis
       , mkSymmAxis
       ) where

import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Generic.Base    as G
import qualified Data.Vector.Generic.Mutable as M

import           Control.Monad               (liftM)
import           Data.Vector                 (Vector)
import           Data.Function               (on)

import           Control.DeepSeq

import           Hammer.Math.Algebra
import           Texture.Orientation

--import           Debug.Trace
--dbg s x = trace (s ++ show x) x

-- =======================================================================================

-- | Defines the an axis of symmetry given an axis 'Vec3' and the number of folds.
newtype SymmAxis
  = SymmAxis (Vec3, Int)
  deriving (Eq, Show)

-- | Defines a fundamental zone plane in Frank-Rodrigues space
newtype FZPlane
  = FZPlane  (Normal3, Double)
  deriving (Eq, Show)

-- | Defines the symmetric group of a certain material.
data Symm = Cubic
          | Hexagonal
          | Custom String [ SymmAxis ]
          deriving (Show, Eq)

-- | Symmetric operators that transforms one orientation to its equivalent.
newtype SymmOp =
  SymmOp
  { symmOp :: Quaternion
  } deriving (Eq, Show, NFData)

-- | Creates a symmetric axis.
mkSymmAxis :: Vec3 -> Int -> SymmAxis
mkSymmAxis v i = SymmAxis (v, i)

-- | Generates a vector of symmetric operators for a given symmetry group.
getSymmOps :: Symm -> U.Vector SymmOp
getSymmOps = U.concatMap mkSymmOps . getSymmAxes

getSymmAxes :: Symm -> U.Vector SymmAxis
getSymmAxes sym = case sym of
  (Custom _ xs) -> U.fromList xs

  Hexagonal -> U.fromList
    [ mkSymmAxis (Vec3   1    0    0 ) 1  -- id rotation

    , mkSymmAxis (Vec3   0    0    1 ) 6

    , mkSymmAxis (Vec3   1    1    1 ) 3
    , mkSymmAxis (Vec3   1  (-1)   1 ) 3
    , mkSymmAxis (Vec3   1    1  (-1)) 3
    , mkSymmAxis (Vec3 (-1)   1    1 ) 3 ]

  Cubic -> U.fromList
    [ mkSymmAxis (Vec3   1    0    0 ) 1  -- id rotation

    , mkSymmAxis (Vec3   0    0    1 ) 4
    , mkSymmAxis (Vec3   0    1    0 ) 4
    , mkSymmAxis (Vec3   1    0    0 ) 4

    , mkSymmAxis (Vec3   1    1    0 ) 2
    , mkSymmAxis (Vec3   1    0    1 ) 2
    , mkSymmAxis (Vec3   0    1    1 ) 2
    , mkSymmAxis (Vec3 (-1)   1    0 ) 2
    , mkSymmAxis (Vec3 (-1)   0    1 ) 2
    , mkSymmAxis (Vec3   0  (-1)   1 ) 2

    , mkSymmAxis (Vec3   1    1    1 ) 3
    , mkSymmAxis (Vec3   1    1  (-1)) 3
    , mkSymmAxis (Vec3   1  (-1)   1 ) 3
    , mkSymmAxis (Vec3 (-1)   1    1 ) 3
    ]

-- | Calculates the symmetric operators for given symmetric axis and the number of
-- symmetric equivalents (number of folds). Number of folds lower than two results
-- in identity operator.
mkSymmOps :: SymmAxis -> U.Vector SymmOp
mkSymmOps (SymmAxis (v, nfold))
  | nfold < 2 = U.singleton (SymmOp zerorot)
  | otherwise = let
    omegaStep = 2 * pi / (fromIntegral nfold)
    mkOne i = let
      omega = ((fromIntegral $ i + 1) * omegaStep)
      in SymmOp $ toQuaternion $ mkRodrigues v (Rad omega)
    in U.generate (nfold-1) mkOne

-- | Given a set of symmetric operators and one orientation find its
-- equivalent orientation that has the minimum rotation angle. This function
-- uses 'composeOmega' that is a specialized and optimized function for calculate
-- rotation angle after rotation composition.
getInFZ :: U.Vector SymmOp -> Quaternion -> Quaternion
getInFZ vs q
  | U.null vs = q
  | otherwise = let
    q0s  = U.map (abs . composeQ0 q . symmOp) vs
    -- Max. |q0| means min. rotation angle.
    minO = U.maxIndex q0s
    in q #<= (symmOp (vs U.! minO))

-- | Finds the symmetric equivalent orientation that belongs to the Fundamental
-- Zone given a symmetric group and one orientation. The orientation in the
-- fundamental zone is the one that has the minimum rotation angle. Note
-- that some rotation representation has antipodal equivalency ('AxisPair'
-- and 'Quaternion'), which means that a rotation of @270 deg@ clockwise
-- is the same as a rotation of @90 deg@ anti-clockwise on the opposite
-- direction. It is faster than 'toFZDirect'.
toFZ :: Symm -> Quaternion -> Quaternion
toFZ symm = getInFZ (getSymmOps symm)

-- | The same as in 'toFZ' but for/to any orientation type. Internally, the orientations
-- are converted to/from 'Quaternion's before the fundamental zone calculation by 'toFZ'.
-- The operations in 'Quaternion's are optimized to reduce the number of calculation
-- operations.
toFZGeneric :: (Rot a, Rot b) => Symm -> a -> b
toFZGeneric symm = fromQuaternion . toFZ symm . toQuaternion

-- | Finds the symmetric equivalent orientation that belongs to the Fundamental
-- Zone given a symmetric group and one orientation.
--
-- Internally the rotation composition ('+\@>') is performed directly for its
-- rotation type without any conversion to 'Quaternion'. Note this operation isn't
-- optimized since it calculates the whole rotation composition before calculating
-- the rotation angle. In most of the cases the function 'toFZ' should be faster.
toFZDirect :: (Rot a)=> V.Vector a -> a -> a
toFZDirect symmOps x = let
  minW  = V.minIndex rs
  rs    = V.map (getAbsShortOmega . (x #<=)) symmOps
  in x #<= (symmOps V.! minW)

-- | Calculates all the symmetric equivalents of a given vector. The calculation is done
-- by passive rotations (changes of base)
getAllSymmVec :: U.Vector SymmOp -> Vec3 -> U.Vector Vec3
getAllSymmVec symOps = \vec -> U.map (passiveVecRotation vec . symmOp) symOps

-- | Calculates all the /unique/ (non-repeated) symmetric equivalents of a given vector.
-- The calculation is done by filtering the output of 'getAllSymmVec'.
getUniqueSymmVecs :: U.Vector SymmOp -> Vec3 -> U.Vector Vec3
getUniqueSymmVecs symOps = nubVec . getAllSymmVec symOps

-- | Get a list of /unique/ vectors by filtering the repeated ones.
nubVec :: (U.Unbox a, DotProd a)=> U.Vector a -> U.Vector a
nubVec l = let
  func (acc, v)
    | U.null v  = acc
    | isInLoop  = func (h `U.cons` acc, U.tail v)
    | otherwise = func (h `U.cons` acc, rest)
    where
      h    = U.head v
      lh   = normsqr h
      rest = U.filter notEqual v
      isInLoop   = U.length v == U.length rest
      notEqual x = abs (h &. x - lh) > 1e-10
  in func (U.empty, l)

-- | Determine whether or not a given quaternion is in the fundamental zone.
isInFZ :: Symm -> Quaternion -> Bool
isInFZ symm = \q -> let
  w = abs (fst $ splitQuaternion q)
  -- test if at least one of the symmetric equivalents has
  -- rotation angle smaller than the input orientation. Therefore the
  -- number of calculations varies between 1 and the number of symmetric operators
  in maybe True (const False) (U.find ((> w) . abs . composeQ0 q . symmOp) os)
  where os = getSymmOps symm

-- =======================================================================================

-- | Get the minimum angular distance touching the closest symmetric plane in the
-- Rodrigues-Frank space.
getMinDistFZPlanes :: Symm -> Double
getMinDistFZPlanes symm
  | U.null as = pi
  | n <= 0    = pi
  | otherwise = pi / (fromIntegral n)
  where
    as = getSymmAxes symm
    SymmAxis (_, n) = U.maximumBy (compare `on` (\(SymmAxis (_, x)) -> x)) as

-- | Get the set of Frank-Rodrigues symmetric planes that form the fundamental zone.
getFZPlanes :: Symm -> U.Vector FZPlane
getFZPlanes = U.map mkFZPlane . getSymmAxes

-- | Creates a Frank-Rodigues symmetric plane.
mkFZPlane :: SymmAxis -> FZPlane
mkFZPlane (SymmAxis (dir, n)) = FZPlane (normal, dist)
  where
    normal = mkNormal dir
    dist   = abs $ tan (pi / (2 * (fromIntegral n)))

-- | Determine whether or not a given quaternion is in the fundamental zone by using
-- Rodrigues-Frank space and its symmetry planes. It should be faster than using all
-- the symmetric operators i.e. 'isInFZ'
isInRodriFZ :: Symm -> Quaternion -> Bool
isInRodriFZ symm = \q -> let
  rq = rodriVec $ fromQuaternion q
  -- test until find the orientation is not inside
  in maybe True (const False) (U.find (not . isInsideRFplane rq) planes)
  where
    -- memorizing planes (eta expression)
    planes = getFZPlanes symm
    -- Test if a given orientation is inside a fundamental zone's plane.
    isInsideRFplane rq (FZPlane (n, dist))
      | proj < -dist = False
      | proj >  dist = False
      | otherwise    = True
      where proj = (fromNormal n) &. rq

-- =======================================================================================

-- | Calculates the minimum misorientation angle between two 'Quaternion's taking into
-- account the crystal symmetry.
--
-- > Gb == Gab * Ga  => Gab == Gb * inv(Ga) == Gb -#- Ga
--
-- By symmetry:
--
-- >  Gb -#- Ga == inv (Ga -#- Gb)
--
-- In order to reduce calculation time, this function uses 24 symmetric operators instead
-- of the theoretical 1152 (24*24*2). At first, it calculates the misorientation and then
-- finds the lowest rotation angle (fundamental zone) among the 24 symmetrical equivalents
-- of the misorientation. Although the mathematical prove is unknown, it works quite well
-- in practice, at least for misorientation /angles/.
--
-- This function should be used for /any/ orientation, even
-- between those in the fundamental zone, due the discontinuity in the fundamental zone.
-- For example, @Euler 44 0 0@ and @Euler 46 0 0@ are @2 degrees@ apart from each other
-- but in the fundamental zone they are @88 degrees@ apart because @Euler 46 0 0@ goes to
-- @Euler (-44) 0 0@ in the fundamental zone.
getMisoAngle :: Symm -> Quaternion -> Quaternion -> Double
getMisoAngle symm = let
  foo = getAbsShortOmega . toFZ symm
  -- avoiding eta expansion of q1 and q2 to memorize
  in \q1 q2 -> foo (q2 -#- q1)

-- -------------------------------------------- Unbox SymmOp ----------------------------------------------------

newtype instance U.MVector s SymmOp = MV_SymmOp (U.MVector s Quaternion)
newtype instance U.Vector    SymmOp = V_SymmOp  (U.Vector    Quaternion)

instance U.Unbox SymmOp

instance M.MVector U.MVector SymmOp where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_SymmOp v)                     = M.basicLength v
  basicUnsafeSlice i n (MV_SymmOp v)            = MV_SymmOp $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_SymmOp v1) (MV_SymmOp v2)   = M.basicOverlaps v1 v2
  basicInitialize (MV_SymmOp v)                 = M.basicInitialize v
  basicUnsafeNew n                              = MV_SymmOp `liftM` M.basicUnsafeNew n
  basicUnsafeReplicate n (SymmOp x)             = MV_SymmOp `liftM` M.basicUnsafeReplicate n x
  basicUnsafeRead (MV_SymmOp v) i               = M.basicUnsafeRead v i >>= (return . SymmOp)
  basicUnsafeWrite (MV_SymmOp v) i (SymmOp x)   = M.basicUnsafeWrite v i x
  basicClear (MV_SymmOp v)                      = M.basicClear v
  basicSet (MV_SymmOp v) (SymmOp x)             = M.basicSet v x
  basicUnsafeCopy (MV_SymmOp v1) (MV_SymmOp v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeGrow (MV_SymmOp v) n               = MV_SymmOp `liftM` M.basicUnsafeGrow v n

instance G.Vector U.Vector SymmOp where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_SymmOp v)             = V_SymmOp `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_SymmOp v)                = MV_SymmOp `liftM` G.basicUnsafeThaw v
  basicLength (V_SymmOp v)                    = G.basicLength v
  basicUnsafeSlice i n (V_SymmOp v)           = V_SymmOp $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_SymmOp v) i            = G.basicUnsafeIndexM v i >>= (return . SymmOp)
  basicUnsafeCopy (MV_SymmOp mv) (V_SymmOp v) = G.basicUnsafeCopy mv v
  elemseq _ (SymmOp x) t                      = G.elemseq (undefined :: Vector a) x t

-- -------------------------------------------- Unbox SymmAxis ----------------------------------------------------

newtype instance U.MVector s SymmAxis = MV_SymmAxis (U.MVector s (Vec3, Int))
newtype instance U.Vector    SymmAxis = V_SymmAxis  (U.Vector    (Vec3, Int))

instance U.Unbox SymmAxis

instance M.MVector U.MVector SymmAxis where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_SymmAxis v)                         = M.basicLength v
  basicUnsafeSlice i n (MV_SymmAxis v)                = MV_SymmAxis $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_SymmAxis v1) (MV_SymmAxis v2)     = M.basicOverlaps v1 v2
  basicInitialize (MV_SymmAxis v)                     = M.basicInitialize v
  basicUnsafeNew n                                    = MV_SymmAxis `liftM` M.basicUnsafeNew n
  basicUnsafeReplicate n (SymmAxis (x,y))             = MV_SymmAxis `liftM` M.basicUnsafeReplicate n (x, y)
  basicUnsafeRead (MV_SymmAxis v) i                   = M.basicUnsafeRead v i >>= (\(x, y) -> return $ SymmAxis (x,y))
  basicUnsafeWrite (MV_SymmAxis v) i (SymmAxis (x,y)) = M.basicUnsafeWrite v i (x, y)
  basicClear (MV_SymmAxis v)                          = M.basicClear v
  basicSet (MV_SymmAxis v) (SymmAxis (x,y))           = M.basicSet v (x, y)
  basicUnsafeCopy (MV_SymmAxis v1) (MV_SymmAxis v2)   = M.basicUnsafeCopy v1 v2
  basicUnsafeGrow (MV_SymmAxis v) n                   = MV_SymmAxis `liftM` M.basicUnsafeGrow v n

instance G.Vector U.Vector SymmAxis where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_SymmAxis v)               = V_SymmAxis `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_SymmAxis v)                  = MV_SymmAxis `liftM` G.basicUnsafeThaw v
  basicLength (V_SymmAxis v)                      = G.basicLength v
  basicUnsafeSlice i n (V_SymmAxis v)             = V_SymmAxis $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_SymmAxis v) i              = G.basicUnsafeIndexM v i >>= (\(x, y) -> return $ SymmAxis (x,y))
  basicUnsafeCopy (MV_SymmAxis mv) (V_SymmAxis v) = G.basicUnsafeCopy mv v
  elemseq _ (SymmAxis (x,y)) t                    = G.elemseq (undefined :: Vector a) x $
                                                    G.elemseq (undefined :: Vector a) y t

-- -------------------------------------------- Unbox FZPlane ----------------------------------------------------

newtype instance U.MVector s FZPlane = MV_FZPlane (U.MVector s (Normal3, Double))
newtype instance U.Vector    FZPlane = V_FZPlane  (U.Vector    (Normal3, Double))

instance U.Unbox FZPlane

instance M.MVector U.MVector FZPlane where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_FZPlane v)                        = M.basicLength v
  basicUnsafeSlice i n (MV_FZPlane v)               = MV_FZPlane $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_FZPlane v1) (MV_FZPlane v2)     = M.basicOverlaps v1 v2
  basicInitialize (MV_FZPlane v)                    = M.basicInitialize v
  basicUnsafeNew n                                  = MV_FZPlane `liftM` M.basicUnsafeNew n
  basicUnsafeReplicate n (FZPlane (x,y))            = MV_FZPlane `liftM` M.basicUnsafeReplicate n (x, y)
  basicUnsafeRead (MV_FZPlane v) i                  = M.basicUnsafeRead v i >>= (\(x, y) -> return $ FZPlane (x,y))
  basicUnsafeWrite (MV_FZPlane v) i (FZPlane (x,y)) = M.basicUnsafeWrite v i (x, y)
  basicClear (MV_FZPlane v)                         = M.basicClear v
  basicSet (MV_FZPlane v) (FZPlane (x,y))           = M.basicSet v (x, y)
  basicUnsafeCopy (MV_FZPlane v1) (MV_FZPlane v2)   = M.basicUnsafeCopy v1 v2
  basicUnsafeGrow (MV_FZPlane v) n                  = MV_FZPlane `liftM` M.basicUnsafeGrow v n

instance G.Vector U.Vector FZPlane where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_FZPlane v)              = V_FZPlane `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_FZPlane v)                 = MV_FZPlane `liftM` G.basicUnsafeThaw v
  basicLength (V_FZPlane v)                     = G.basicLength v
  basicUnsafeSlice i n (V_FZPlane v)            = V_FZPlane $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_FZPlane v) i             = G.basicUnsafeIndexM v i >>= (\(x, y) -> return $ FZPlane (x,y))
  basicUnsafeCopy (MV_FZPlane mv) (V_FZPlane v) = G.basicUnsafeCopy mv v
  elemseq _ (FZPlane (x,y)) t                   = G.elemseq (undefined :: Vector a) x $
                                                  G.elemseq (undefined :: Vector a) y t
