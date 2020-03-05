{-# LANGUAGE
    NamedFieldPuns
  , RecordWildCards
  , GeneralizedNewtypeDeriving
  , TypeFamilies
  , MultiParamTypeClasses
  , TemplateHaskell
  #-}
-- |
-- Module      : Texture.Symmetry
-- Copyright   : (c) 2013 Edgar Gomes
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
  , getMisoAngleFaster 
    -- * Fundamental zone in Frank-Rodrigues
  , isInRodriFZ
   -- * Custom symmetries
  , getSymmOps
  , SymmOp
  , symmOp
  , mkSymmOps
  , SymmAxis
  , mkSymmAxis
    -- * Averaging
  , averageQuaternionWithSymm
  ) where

import Control.DeepSeq
import Control.Monad
import Data.Function (on)
import Data.Maybe (isNothing)
import Data.Vector.Unboxed.Deriving
import Linear.Vect
import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import Texture.Orientation

-- =======================================================================================

-- | Defines the an axis of symmetry given an axis 'Vec3' and the number of folds.
newtype SymmAxis
  = SymmAxis (Vec3D, Int)
  deriving (Eq, Show)

-- | Defines a fundamental zone plane in Frank-Rodrigues space
newtype FZPlane
  = FZPlane  (Normal3D, Double)
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
mkSymmAxis :: Vec3D -> Int -> SymmAxis
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
  | nfold < 2 = U.singleton (SymmOp mempty)
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
getAllSymmVec :: U.Vector SymmOp -> Vec3D -> U.Vector Vec3D
getAllSymmVec symOps = \vec -> U.map (passiveVecRotation vec . symmOp) symOps

-- | Calculates all the /unique/ (non-repeated) symmetric equivalents of a given vector.
-- The calculation is done by filtering the output of 'getAllSymmVec'.
getUniqueSymmVecs :: U.Vector SymmOp -> Vec3D -> U.Vector Vec3D
getUniqueSymmVecs symOps = nubVec . getAllSymmVec symOps

-- TODO: Is there a better complexity for this size?
-- | Get a list of /unique/ vectors by filtering the repeated ones.
nubVec :: (Fractional a, Ord a, U.Unbox (v a), DotProd a v)=> U.Vector (v a) -> U.Vector (v a)
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
  in isNothing $ U.find ((> w) . abs . composeQ0 q . symmOp) os
  where os = getSymmOps symm

-- =======================================================================================

-- | Get the minimum angular distance touching the closest symmetric plane in the
-- Rodrigues-Frank space.
getMinDistFZPlanes :: Symm -> Double
getMinDistFZPlanes symm
  | U.null as = pi
  | n <= 0    = pi
  | otherwise = pi / fromIntegral n
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
    dist   = abs . tan $ pi / (2 * fromIntegral n)

-- | Determine whether or not a given quaternion is in the fundamental zone by using
-- Rodrigues-Frank space and its symmetry planes. It should be faster than using all
-- the symmetric operators i.e. 'isInFZ'
isInRodriFZ :: Symm -> Quaternion -> Bool
isInRodriFZ symm = \q -> let
  rq = rodriVec $ fromQuaternion q
  -- test until find the orientation is not inside
  in isNothing $ U.find (not . isInsideRFplane rq) planes
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
  ops = getSymmOps symm
  -- avoiding eta expansion of q1 and q2 to memorize
  in \q1 q2 -> getMisoAngleFaster ops q1 q2

getMisoAngleFaster :: U.Vector SymmOp -> Quaternion -> Quaternion -> Double
getMisoAngleFaster ops q1 q2 = getAbsShortOmega $ getInFZ ops (q2 -#- q1)

-- ================================== Averaging quaternions regarding symmetry ===========================================

averageQuaternionWithSymm :: Foldable t => Symm -> t Quaternion -> Quaternion
averageQuaternionWithSymm symm vq
  | null vq   = mempty
  | otherwise = averageWeightedQuaternion . bringTogether $ avgSites
  where
    os       = getSymmOps symm
    avgSites = foldl addOnBucket [] vq
    omega    = cos $ 0.5 * fromAngle (Deg 5)

    addOnBucket :: [(Quaternion, Vec4D)] -> Quaternion -> [(Quaternion, Vec4D)]
    addOnBucket buckets q = go [] buckets
      where
        go acc (x : xs) = maybe (go (x : acc) xs) ((acc ++) .  (: xs)) (simpleAvg x q)
        go acc _        = (q, quaterVec q) : acc

    simpleAvg :: (Quaternion, Vec4D) -> Quaternion -> Maybe (Quaternion, Vec4D)
    simpleAvg (ref, acc) q = let
      q0 = composeQ0 (invert ref) q
      -- Doesn't need to check which representation (q or antipodal q) is the closest to
      -- the accumulator (acc &. q >= 0) because Quaternion enforces the use of half the
      -- space (only q0 > 0).
      in guard (q0 > omega) >> pure (ref, quaterVec q &+ acc)

    bringTogether xs = case L.sortBy (flip compare `on` (vlen . snd)) xs of
      ((_, q) : qs) -> (vlen q, mkQuaternion q) : map (findClosest (mkQuaternion q) . snd) qs
      _             -> []

    findClosest :: Quaternion -> Vec4D -> (Double, Quaternion)
    findClosest ref v = let
      vl  = vlen v
      q   = mkQuaternion v
      qs  = U.map ((q #<=) . symmOp) os
      ms  = U.map (composeQ0 (invert ref)) qs
      i   = U.maxIndex ms
      vi  = qs U.! i
      in (vl, vi)

-- ============================================ Unbox SymmOp ====================================================

derivingUnbox "SymmOp"
    [t| SymmOp -> Quaternion |]
    [| \(SymmOp q) -> q |]
    [| \q -> SymmOp q |]

derivingUnbox "SymmAxis"
    [t| SymmAxis -> (Vec3D, Int) |]
    [| \(SymmAxis a) -> a |]
    [| \a -> SymmAxis a |]

derivingUnbox "FZPlane"
    [t| FZPlane -> (Normal3D, Double) |]
    [| \(FZPlane a) -> a |]
    [| \a -> FZPlane a |]
