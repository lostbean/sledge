{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

-- |
-- Module      : Hammer.Texture.IPF
-- Copyright   : (c) 2013 Edgar Gomes
-- License     : Privete-style (see LICENSE)
--
-- Maintainer  : Edgar Gomes <talktoedgar@gmail.com>
-- Stability   : experimental
-- Portability : tested on GHC only
--
-- This module defines IPF (Inverse Pole Figure) representation of orientations.
--
module Hammer.Texture.IPF
       ( IPFSchema
       , getIPFSchema
       , invPole
       , ipfColor
       , getIPFColor
       , getCubicIPFColor
       ) where

import qualified Data.Vector as V

import           Data.Vector (Vector)

import           Hammer.Math.Algebra
import           Hammer.Texture.Orientation
import           Hammer.Texture.Symmetry

import           Debug.Trace
dbg s x = trace (s ++ show x) x

-- =======================================================================================

-- | Defines the IPF color schema. Each color (Red, Blue and Green) is associated with a
-- list family of crystalographic directions.
data IPFSchema =
  IPFSchema
  { redBase   :: Vector Normal3
  , blueBase  :: Vector Normal3
  , greenBase :: Vector Normal3
  } deriving (Show)

-- | Get the 'IPFSchema' for a given group symmetry. It should be evaluated only once
-- since the schema is constant for a certain symmetry group.
getIPFSchema :: Symm -> IPFSchema
getIPFSchema symm = case symm of
  Cubic -> IPFSchema
           (cubicFamilyDir (Vec3 0 0 1))
           (cubicFamilyDir (Vec3 0 1 1))
           (cubicFamilyDir (Vec3 1 1 1))
  x     -> error $ "[IPF] IPF schema not defined for " ++ show x
  where
    cubicFamilyDir = V.map mkNormal . getUniqueSymmVecs Cubic

-- | Calculate the inverse poles of given orientation. The inverse pole is the
-- crystalographic direction parallel to one of the external frame directions.
invPole :: RefFrame -> Quaternion -> Normal3
invPole refdir = mkNormal . case refdir of
  RD -> _1 . m
  TD -> _2 . m
  ND -> _3 . m
  where
    m = inverse . rotMat . fromQuaternion

-- | Calculates the IPF color of a direction in a certain symmetry group. The color is
-- represented in a RBG color schema where which one of the color components is associated
-- with one direction family. Then the minimum angle is calculated between each direction
-- family and the given direction. The color associated with lowest angle value (let call
-- it @alpha1@) is set to @1@ and the other two colors are given by ratio between its
-- angle and @alpha1@ (e.g. @red = alpha2/alpha1@).
ipfColor :: Symm -> Normal3 -> (Double, Double, Double)
ipfColor symm n
  | angleR < angleB && angleR < angleG = dbg2 (1, angleR / angleB, angleR / angleG)
  | angleB < angleR && angleB < angleG = dbg2 (angleB / angleR, 1, angleB / angleG)
  | otherwise                          = dbg2 (angleG / angleR, angleG / angleB, 1)
  where
    dbg2 = traceShow ( toDeg angleR,toDeg angleB,toDeg angleG )
    IPFSchema{..} =  getIPFSchema symm
    angleR = findMinAngle n redBase
    angleB = findMinAngle n blueBase
    angleG = findMinAngle n greenBase

-- | Calculates the IPF color for one of reference frame directions of a given quaternion
-- and symmetry group.
getIPFColor :: Symm -> RefFrame -> Quaternion -> (Double, Double, Double)
getIPFColor symm ref = ipfColor symm .  invPole ref

-- | Calculates the IPF color for one of reference frame directions of a given quaternion
-- in the /Cubic/ group.
getCubicIPFColor :: RefFrame -> Quaternion -> (Double, Double, Double)
getCubicIPFColor ref = ipfColor Cubic .  invPole ref . toFZ Cubic

-- | Finds the minimum angle in radians between one vector and a list of vectors.
findMinAngle :: Normal3 -> Vector Normal3 -> Double
findMinAngle t = V.minimum . V.map (angle' t)
