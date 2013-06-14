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
       , toFZQuaternion
       , toFZCustom
       , mkSymmOps
       , getSymmOps
       , SymmOp (symmOp)
       , getAllSymmVec
       , getUniqueSymmVecs
       ) where

import qualified Data.Vector as V

import           Data.Vector    (Vector)

import           Data.List
import           Hammer.Math.Algebra
import           Hammer.Texture.Orientation

import           Debug.Trace
dbg s x = trace (s ++ show x) x

-- | Defines the symmetric group of a certain material.
data Symm = Cubic
          | Hexagonal
          | Custom String [ SymmOp ]
          deriving (Show, Eq)

-- | Symmetric operators that transforms one orientation to its equivalent.
newtype SymmOp =
  SymmOp
  { symmOp :: Quaternion
  } deriving (Show, Eq)

-- | Generates a vector of symmetric operators for a given symmetry group.
getSymmOps :: Symm -> Vector SymmOp
getSymmOps sym = case sym of
  (Custom _ xs) -> V.fromList xs

  Hexagonal -> V.concat
    [ mkSymmOps (Vec3   1    0    0 ) 1  -- id rotation

    , mkSymmOps (Vec3   0    0    1 ) 6

    , mkSymmOps (Vec3   1    1    1 ) 3
    , mkSymmOps (Vec3   1  (-1)   1 ) 3
    , mkSymmOps (Vec3   1    1  (-1)) 3
    , mkSymmOps (Vec3 (-1)   1    1 ) 3 ]

  Cubic -> V.concat
    [ mkSymmOps (Vec3   1    0    0 ) 1  -- id rotation

    , mkSymmOps (Vec3   0    0    1 ) 4
    , mkSymmOps (Vec3   0    1    0 ) 4
    , mkSymmOps (Vec3   1    0    0 ) 4

    , mkSymmOps (Vec3   1    1    0 ) 2
    , mkSymmOps (Vec3   1    0    1 ) 2
    , mkSymmOps (Vec3   0    1    1 ) 2
    , mkSymmOps (Vec3 (-1)   1    0 ) 2
    , mkSymmOps (Vec3 (-1)   0    1 ) 2
    , mkSymmOps (Vec3   0  (-1)   1 ) 2

    , mkSymmOps (Vec3   1    1    1 ) 3
    , mkSymmOps (Vec3   1    1  (-1)) 3
    , mkSymmOps (Vec3   1  (-1)   1 ) 3
    , mkSymmOps (Vec3 (-1)   1    1 ) 3
    ]

-- | Calculates the symmetric operators for given symmetric axis and the number of
-- symmetric equivalents (number of folds). Number of folds lower than two results
-- in identity operator.
mkSymmOps :: Vec3 -> Int -> Vector SymmOp
mkSymmOps v nfold
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
    symRot = (q +@>) . symmOp
    omegas = V.map (composeOmega q . symmOp) vs
    minOme = V.minIndex omegas
    in symRot (vs V.! minOme)

-- | Fundamental zone calculation for 'Quaternion'. Preferable use 'toFZ' since its
-- works for any 'Rot'ation.
toFZQuaternion :: Symm -> Quaternion -> Quaternion
toFZQuaternion symm = getInFZ (getSymmOps symm)

-- | Finds the symmetric equivalent orientation that belongs to the Fundamental
-- Zone given a symmetric group and one orientation. The orientation in the
-- fundamental zone is the one that has the minimum rotation angle. Note
-- that some rotation representation has antipodal equivalency ('AxisPair'
-- and 'Quaternion'), which means that a rotation of @271 deg@ clockwise
-- is the same as a rotation of @90 deg@ anti-clockwise on the opposite
-- direction.
--
-- Internally any orientation is converted to 'Quaternion' for the fundamental
-- zone calculation. Operations in 'Quaternion' is optimized to reduce
-- calculation operations.
toFZ :: (Rot a, Rot b) => Symm -> a -> b
toFZ symm = fromQuaternion . toFZQuaternion symm . toQuaternion


-- | Finds the symmetric equivalent orientation that belongs to the Fundamental
-- Zone given a symmetric group and one orientation.
--
-- Internally the rotation composition ('+\@>') is performed directly for its
-- rotation type without any conversion to 'Quaternion'. Note this operation isn't
-- optimized since it calculates the whole rotation composition before calculating
-- the rotation angle. In most of the cases the function 'toFZ' should be faster.
toFZCustom :: (Rot a)=> Vector a -> a -> a
toFZCustom symmOps x = let
  minW  = V.minIndex rs
  rs    = V.map (getOmegaMinPodal . (x +@>)) symmOps
  in x +@> (symmOps V.! minW)

-- | Finds the minimum absolute rotation angle between both antipodal equivalent
-- rotations. Note that some rotation representation has antipodal equivalency
-- ('AxisPair' and 'Quaternion'), which means that a rotation of @271 deg@
-- clockwise is the same as a rotation of @90 deg@ anti-clockwise on the opposite
-- direction.
getOmegaMinPodal :: (Rot a)=> a -> Double
getOmegaMinPodal x
  | w > pi    = abs (w - 2*pi)
  | otherwise = w
  where w = getOmega x

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