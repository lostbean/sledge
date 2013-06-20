{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}

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
-- Calculation of IPF color
-- The IPF color is represented in the RGB color schema where which one of the color
-- components, named Red, Green and Blue is associated with one direction family
-- (e.g. Red ~> [0 0 1]) and those direction must be linear independent. The three
-- directions form a set of three planes where each plane contains two of the three
-- directions. For example, if @vR@ ,@vG@ and @vB@ are the directions associated with the
-- RGB color then the @planeRG@ is defined by its normal @planeRG = vG x vR@. The same is
-- truth for @planeGB = vB x vG@ and @planeBR = vR x vB@. The set of three planes is
-- represented by 'IPFBase'. Then a set of three angles is calculated between the planes
-- and the given pole for each symmetric equivalent of 'IPFBase' (set of planes) and the
-- set containing the minimum values of angles is selected.
--
-- The color associated with largest angle value is set to @1@. (e.g. If the angle between
-- the @planeRG@ and the pole is the largest than the /blue/ is set to @1@) and the other
-- two colors are given by ratio between its correspondent angle and the largest angle
-- (e.g. @red = alphaBR/alphaRG@ and @blue = alphaGB/alphaRG@).
--
module Hammer.Texture.IPF
       ( IPFBase
       , getIPFSchema
       , invPole
       , ipfColor
       , getIPFColor
       , getIPFColorNoFZ
       , genIPFLegend
       ) where

import qualified Data.Vector                     as V

import           System.Random

import           Hammer.Math.Algebra
import           Hammer.Texture.Orientation
import           Hammer.Texture.SphereProjection
import           Hammer.Texture.Symmetry

--import           Debug.Trace
--dbg s x = trace (s ++ show x) x

-- =======================================================================================

type RGBColor = (Double, Double, Double)

-- | Each color (Red, Green and Blue) is associated with a linear independent direction.
-- The 'IPFBase' is the set of three planes that contains two of the three reference
-- directions. Each plane is defined by its normal direction. For example, let @vR@, @vG@
-- and @vB@ be the direction vectors of red, green and blue respectively. Then one of
-- the planes, @planeRG@, is defined by @planeRG = vG x vR@.
data IPFBase =
  IPFBase
  { planeRG :: Normal3
  , planeGB :: Normal3
  , planeBR :: Normal3
  } deriving (Show)

-- | Calculates a 'IPFSchema' for a given symmetry group and three reference directions
-- associated with the colors in the following sequence: @Red -> Green -> Blue@. In order
-- to find the right fundamental zone, the set of reference directions have to be adjusted
-- by choosing right members of the direction's family (e.g. <0 1 0> from {1 1 0} family)
mkIPFSchema :: Vec3 -> Vec3 -> Vec3 -> IPFBase
mkIPFSchema vecR vecG vecB = let
  normalRG = mkNormal $ vecG &^ vecR
  normalGB = mkNormal $ vecB &^ vecG
  normalBR = mkNormal $ vecR &^ vecB
  in IPFBase normalRG normalGB normalBR

-- | Get the 'IPFSchema' for a given group symmetry. It should be evaluated only once
-- since the schema is constant for a certain symmetry group.
getIPFSchema :: Symm -> IPFBase
getIPFSchema symm = case symm of
  Cubic -> mkIPFSchema (Vec3 0 0 1) (Vec3 1 0 1) (Vec3 1 1 1)
  x     -> error $ "[IPF] IPF schema not defined for " ++ show x

-- | Calculate the inverse poles of given orientation. The inverse pole is the
-- crystalographic direction parallel to one of the external frame directions.
invPole :: RefFrame -> Quaternion -> Normal3
invPole refdir = mkNormal . case refdir of
  RD -> _1 . m
  TD -> _2 . m
  ND -> _3 . m
  where
    -- The columns of the orientation matrix are the directions parallel to RD, TD, and ND
    -- (m = [RD, TD, ND]) and the rows represent the directions of the rotated frame (
    -- crystal frame ) axis.
    m = inverse . rotMat . fromQuaternion

colorAdjust :: RGBColor -> RGBColor
colorAdjust (r,g,b) = (adj r, adj g, adj b)
  where adj x = 0.15 + 0.85 * x

-- | Calculates the IPF color of a pole direction in a certain symmetry group.
ipfColor :: (Double, Double, Double) -> RGBColor
ipfColor (angleB, angleR, angleG)
  | angleR > angleG && angleR > angleB = (1, angleG / angleR, angleB / angleR)
  | angleG > angleR && angleG > angleB = (angleR / angleG, 1, angleB / angleG)
  | otherwise                          = (angleR / angleB, angleG / angleB, 1)

-- | Calculates the IPF color in one reference frame directions for a given orientation
--  represented in quaternion and its symmetry group.
getIPFColor :: Symm -> RefFrame -> Quaternion -> (Normal3, RGBColor)
getIPFColor symm ref q = let
  n         = invPole ref q
  (fzn, as) = findMinAngleBase symm n
  in (fzn, colorAdjust $ ipfColor as)

-- | Calculates the IPF color in one reference frame directions for a given orientation
--  represented in quaternion and its symmetry group.
getIPFColorNoFZ :: Symm -> RefFrame -> Quaternion -> (Normal3, RGBColor)
getIPFColorNoFZ symm ref q = let
  n       = invPole ref q
  (_, as) = findMinAngleBase symm n
  in (n, colorAdjust $ ipfColor as)

-- | Finds the set of angles between a given direction and the planes of an 'IPFBase'. The
-- returned angles between the direction @v@ and the planes are in radians and in the
-- following order: @(v <-> plane RG, v <-> plane GB, v <-> plane BR)@
findAnglesBase :: Normal3 -> IPFBase -> (Double, Double, Double)
findAnglesBase v IPFBase{..} = let
  func plane x
    | alpha > 0.5 * pi = alpha - 0.5 * pi
    | otherwise        = 0.5 * pi - alpha
    where alpha = angle' plane x
  in (func planeRG v, func planeGB v, func planeBR v)

-- | Finds the minimum set of angles between a given direction and all symmetric planes of
-- 'IPFBase'. The returned angles between the direction @v@ and the planes are in radians
-- and in the following order: @(v <-> plane RG, v <-> plane GB, v <-> plane BR)@
findMinAngleBase :: Symm -> Normal3 -> (Normal3, (Double, Double, Double))
findMinAngleBase symm x = let
  base   = getIPFSchema symm
  symmX  = V.map mkNormal $ getAllSymmVec symm $ fromNormal x
  ibase  = V.minIndex $ V.map (angDist . (flip findAnglesBase) base) symmX
  finalX = symmX V.! ibase
  alphas = findAnglesBase finalX base
  angDist (a1, a2, a3) = a1*a1 + a2*a2 + a3*a3
  in (finalX , alphas)

genIPFLegend :: FilePath -> Int -> IO ()
genIPFLegend name n = do
  qs <- V.replicateM n randomIO
  let
    -- (ns, cs) = V.unzip $ V.map (getIPFColorNoFZ Cubic ND) qs
    (ns, cs) = V.unzip $ V.map (getIPFColor Cubic ND) qs
    sps      = getBothProj $ V.map lambertSO3Proj ns
    out      = V.zipWith func sps cs
    func (Vec2 x y) (r, g, b) = show x ++ " " ++ show y ++ " " ++
                                show r ++ " " ++ show g ++ " " ++ show b
  writeFile name (unlines $ V.toList out)
