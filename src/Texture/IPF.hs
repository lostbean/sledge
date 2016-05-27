{-# LANGUAGE
    NamedFieldPuns
  , RecordWildCards
  #-}
-- |
-- Module      : Texture.IPF
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
module Texture.IPF
 ( invPole
 , getIPFColor
 , genIPFLegend
 , getRGBColor
 , RGBColor       (..)
 , RGBDoubleColor (..)
 ) where

import Data.Word           (Word8)
import Data.Vector.Unboxed (Vector)
import Linear.Vect
import System.Random
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector         as V

import Texture.Orientation
import Texture.SphereProjection
import Texture.Symmetry

-- =======================================================================================

newtype RGBDoubleColor = RGBDoubleColor (Double, Double, Double)

newtype RGBColor  = RGBColor  (Word8, Word8, Word8)

newtype AngularDist = AngularDist (Double, Double, Double)

-- | Each color (Red, Green and Blue) is associated with a linear independent direction.
-- The 'IPFBase' is the set of three planes that contains two of the three reference
-- directions. Each plane is defined by its normal direction. For example, let @vR@, @vG@
-- and @vB@ be the direction vectors of red, green and blue respectively. Then one of
-- the planes, @planeRG@, is defined by @planeRG = vG x vR@.
data IPFBase =
  IPFBase
  { planeRG :: Normal3D
  , planeGB :: Normal3D
  , planeBR :: Normal3D
  } deriving (Show)

-- | Calculates a 'IPFSchema' for a given symmetry group and three reference directions
-- associated with the colors in the following sequence: @Red -> Green -> Blue@. In order
-- to find the right fundamental zone, the set of reference directions have to be adjusted
-- by choosing right members of the direction's family (e.g. <0 1 0> from {1 1 0} family)
mkIPFSchema :: Vec3D -> Vec3D -> Vec3D -> IPFBase
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
invPole :: RefFrame -> Quaternion -> Normal3D
invPole refdir = mkNormal . case refdir of
  RD -> _R1 . m
  TD -> _R2 . m
  ND -> _R3 . m
  where
    -- The columns of the orientation matrix are the directions parallel to RD, TD, and ND
    -- (m = [RD, TD, ND]) and the rows represent the directions of the rotated frame (
    -- crystal frame ) axis.
    m = inverse . rotMat . fromQuaternion

colorAdjust :: RGBDoubleColor -> RGBDoubleColor
colorAdjust (RGBDoubleColor (r,g,b)) = RGBDoubleColor (adj r, adj g, adj b)
  where adj x = 0.15 + 0.85 * x

getRGBColor :: RGBDoubleColor -> RGBColor
getRGBColor (RGBDoubleColor (r,g,b)) = let
  ru = round $ r * 255
  gu = round $ g * 255
  bu = round $ b * 255
  in RGBColor (ru, gu, bu)

-- | Calculates the IPF color of a pole direction in a certain symmetry group.
ipfColor :: AngularDist -> RGBDoubleColor
ipfColor (AngularDist (angleB, angleR, angleG))
  | angleR > angleG && angleR > angleB = RGBDoubleColor (1, angleG / angleR, angleB / angleR)
  | angleG > angleR && angleG > angleB = RGBDoubleColor (angleR / angleG, 1, angleB / angleG)
  | otherwise                          = RGBDoubleColor (angleR / angleB, angleG / angleB, 1)

-- | Calculates the IPF color in one reference frame directions for a given orientation
--  represented in quaternion and its symmetry group.
getIPFColor :: Symm -> RefFrame -> Quaternion -> (Normal3D, RGBDoubleColor)
getIPFColor symm ref = let
  foo (fzn, as) = (fzn, colorAdjust $ ipfColor as)
  base          = getIPFSchema symm
  symOps        = getSymmOps symm
  in foo . findMinAngleBase symOps base . invPole ref

-- | Finds the set of angles between a given direction and the planes of an 'IPFBase'. The
-- returned angles between the direction @v@ and the planes are in radians and in the
-- following order: @(v <-> plane RG, v <-> plane GB, v <-> plane BR)@
findAnglesBase :: Normal3D -> IPFBase -> AngularDist
findAnglesBase v IPFBase{..} = let
  func plane x
    | alpha > 0.5 * pi = alpha - 0.5 * pi
    | otherwise        = 0.5 * pi - alpha
    where alpha = angle' plane x
  in AngularDist (func planeRG v, func planeGB v, func planeBR v)

-- | Finds the minimum set of angles between a given direction and all symmetric planes of
-- 'IPFBase'. The returned angles between the direction @v@ and the planes are in radians
-- and in the following order: @(v <-> plane RG, v <-> plane GB, v <-> plane BR)@
findMinAngleBase :: Vector SymmOp -> IPFBase -> Normal3D -> (Normal3D, AngularDist)
findMinAngleBase symOps base x = let
  symmX  = U.map mkNormal $ getAllSymmVec symOps $ fromNormal x
  ibase  = U.minIndex $ U.map (angDist . flip findAnglesBase base) symmX
  finalX = symmX U.! ibase
  alphas = findAnglesBase finalX base
  angDist (AngularDist (a1, a2, a3)) = a1*a1 + a2*a2 + a3*a3
  in (finalX , alphas)

genIPFLegend :: FilePath -> Int -> IO ()
genIPFLegend name n = do
  qs <- V.replicateM n randomIO
  let
    (ns, cs) = V.unzip $ V.map (getIPFColor Cubic ND) qs
    sps      = getBothProj $ V.map lambertSO3Proj ns
    out      = V.zipWith func sps cs
    func (Vec2 x y) (RGBDoubleColor (r, g, b)) =
      show x ++ " " ++ show y ++ " " ++ show r ++ " " ++ show g ++ " " ++ show b
  writeFile name (unlines $ V.toList out)
