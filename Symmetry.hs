{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

-- |
-- Module      : Hammer.Texture.Symmetry
-- Copyright   : (c) 2013 Edgar Gomes
-- License     : Privete-style (see LICENSE)
--
-- Maintainer  : Edgar Gomes <talktoedgar@gmail.com>
-- Stability   : experimental
-- Portability : tested on GHC only
--
-- Module to calculate symmetry in crystal orientations.
--
module Hammer.Texture.Symmetry
       ( Symm (..)
       , toFZ
       , toFZGeneric
       , toFZDirect
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

import qualified Data.Vector as V

import           Data.Vector    (Vector)

import           Hammer.Math.Algebra
import           Hammer.Texture.Orientation

--import           Debug.Trace
--dbg s x = trace (s ++ show x) x

-- =======================================================================================

-- | Defines the an axis of symmetry given an axis 'Vec3' and the number of folds.
newtype SymmAxis = SymmAxis (Vec3, Int)       deriving (Show, Eq)

-- | Defines a fundamental zone plane in Frank-Rodrigues space
newtype FZPlane  = FZPlane  (Normal3, Double) deriving (Show, Eq)

-- | Defines the symmetric group of a certain material.
data Symm = Cubic
          | Hexagonal
          | Custom String [ SymmAxis ]
          deriving (Show, Eq)

-- | Symmetric operators that transforms one orientation to its equivalent.
newtype SymmOp =
  SymmOp
  { symmOp :: Quaternion
  } deriving (Show, Eq)

-- | Creates a symmetric axis.
mkSymmAxis :: Vec3 -> Int -> SymmAxis
mkSymmAxis v i = SymmAxis (v, i)

-- | Generates a vector of symmetric operators for a given symmetry group.
getSymmOps :: Symm -> Vector SymmOp
getSymmOps = V.concatMap mkSymmOps . getSymmAxes

getSymmAxes :: Symm -> Vector SymmAxis
getSymmAxes sym = case sym of
  (Custom _ xs) -> V.fromList xs

  Hexagonal -> V.fromList
    [ mkSymmAxis (Vec3   1    0    0 ) 1  -- id rotation

    , mkSymmAxis (Vec3   0    0    1 ) 6

    , mkSymmAxis (Vec3   1    1    1 ) 3
    , mkSymmAxis (Vec3   1  (-1)   1 ) 3
    , mkSymmAxis (Vec3   1    1  (-1)) 3
    , mkSymmAxis (Vec3 (-1)   1    1 ) 3 ]

  Cubic -> V.fromList
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
mkSymmOps :: SymmAxis -> Vector SymmOp
mkSymmOps (SymmAxis (v, nfold))
  | nfold < 2 = V.singleton (SymmOp zerorot)
  | otherwise = let
    omegaStep = 2 * pi / (fromIntegral nfold)
    mkOne i = let
      omega = ((fromIntegral $ i + 1) * omegaStep)
      in SymmOp $ toQuaternion $ mkRodrigues v (Rad omega)
    in V.generate (nfold-1) mkOne

-- | Given a set of symmetric operators and one orientation find its
-- equivalent orientation that has the minimum rotation angle. This function
-- uses 'composeOmega' that is a specialized and optimized function for calculate
-- rotation angle after rotation composition.
getInFZ :: Vector SymmOp -> Quaternion -> Quaternion
getInFZ vs q
  | V.null vs = q
  | otherwise = let
    minO = V.maxIndex q0s
    comp = (q #<=) . symmOp
    q0s  = V.map (abs . fst . splitQuaternion . comp) vs
    in comp (vs V.! minO)

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
toFZDirect :: (Rot a)=> Vector a -> a -> a
toFZDirect symmOps x = let
  minW  = V.minIndex rs
  rs    = V.map (getOmegaRange . getOmega . (x #<=)) symmOps
  in x #<= (symmOps V.! minW)

-- | Finds the minimum absolute rotation angle between both antipodal equivalent
-- rotations. Note that some rotation representation has antipodal equivalency
-- ('AxisPair' and 'Quaternion'), which means that a rotation of @271 deg@
-- clockwise is the same as a rotation of @90 deg@ anti-clockwise on the opposite
-- direction.
getOmegaRange :: Double -> Double
getOmegaRange o
  | o > pi    = 2 * pi - o
  | otherwise = o

-- | Calculates all the symmetric equivalents of a given vector. The calculation is done
-- by passive rotations (changes of base)
getAllSymmVec :: Symm -> Vec3 -> Vector Vec3
getAllSymmVec symm vec = let
  syops = getSymmOps symm
  in V.map (passiveVecRotation vec . symmOp) syops

-- | Calculates all the /unique/ (non-repeated) symmetric equivalents of a given vector.
-- The calculation is done by filtering the output of 'getAllSymmVec'.
getUniqueSymmVecs :: Symm -> Vec3 -> Vector Vec3
getUniqueSymmVecs symm v = nubVec $ getAllSymmVec symm v

-- | Get a list of /unique/ vectors by filtering the repeated ones.
nubVec :: (DotProd a)=> Vector a -> Vector a
nubVec l = let
  func (acc, v)
    | V.null v  = acc
    | isInLoop  = func (h `V.cons` acc, V.tail v)
    | otherwise = func (h `V.cons` acc, rest)
    where
      h    = V.head v
      lh   = normsqr h
      rest = V.filter notEqual v
      isInLoop   = V.length v == V.length rest
      notEqual x = abs (h &. x - lh) > 1e-10
  in func (V.empty, l)

-- =======================================================================================

-- | Get the set of Frank-Rodrigues symmetric planes that form the fundamental zone.
getFZPlanes :: Symm -> Vector FZPlane
getFZPlanes = V.map mkFZPlane . getSymmAxes

-- | Creates a Frank-Rodigues symmetric plane.
mkFZPlane :: SymmAxis -> FZPlane
mkFZPlane (SymmAxis (dir, n)) = FZPlane (normal, dist)
  where
    normal = mkNormal dir
    dist   = abs $ tan (pi / (2 * (fromIntegral n)))

-- | Test if a given orientation is inside the fundamental zone.
isInRodriFZ :: Symm -> Quaternion -> Bool
isInRodriFZ symm q = V.and . V.map (isInRodriFZPlane q) $ getFZPlanes symm

-- | Test if a given orientation is between a fundamental zone symmetric plane.
isInRodriFZPlane :: Quaternion -> FZPlane -> Bool
isInRodriFZPlane q (FZPlane (n, dist))
  | proj < -dist = False
  | proj >  dist = False
  | otherwise    = True
  where
    rq = rodriVec $ fromQuaternion q
    proj = (fromNormal n) &. rq

-- =======================================================================================

-- | Calculates the minimum misorientation angle between two 'Quaternion's taking into
-- account the crystal symmetry.
--
-- > Gb == Gab * Ga  => Gab == Gb * inv(Ga) == Gb -#- Ga
-- By symmetry:
--
-- >  Gb -#- Ga == inv (Ga -#- Gb)
--
-- In order to reduce calculation time, this function uses 24 symmetric operators instead of
-- the theoretical 1152 (24*24*2). At first, it calculates the misorientation and then finds
-- the lowest rotation angle (fundamental zone) among the 24 symmetrical equivalents of the
-- misorientation. Although the mathematical prove is unknown, it works quite well in
-- practice, at least for misorientation /angles/.
--
-- This function should be used for /any/ orientation, even
-- between those in the fundamental zone, due the discontinuity in the fundamental zone.
-- For example, @Euler 44 0 0@ and @Euler 46 0 0@ are @2 degrees@ apart in the reality but
-- the fundamental zone they @88 degrees@ apart because @Euler 46 0 0@ goes to
-- @Euler (-44) 0 0@ in the fundamental zone.
getMisoAngle :: Symm -> Quaternion -> Quaternion -> Double
getMisoAngle symm q1 q2 = getOmegaRange . getOmega $ toFZ symm (q2 -#- q1)

