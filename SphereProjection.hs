{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}

-- |
-- Module      : Hammer.Texture.SphereProjection
-- Copyright   : (c) 2013 Edgar Gomes
-- License     : Privete-style (see LICENSE)
--
-- Maintainer  : Edgar Gomes <talktoedgar@gmail.com>
-- Stability   : experimental
-- Portability : tested on GHC only
--
-- This module defines SO3 sphere projections.
--
module Hammer.Texture.SphereProjection
       ( lambertSO3Proj
       , steroSO3Proj
       , getLowerProj
       , getUpperProj
       , getBothProj
       ) where

import qualified Data.Vector                as V

import           Data.Vector                  (Vector)

import           Hammer.Math.Algebra

-- =======================================================================================

-- | Define symmetry types. Reflection through the origin (v = -v) or no symmetry at all.
data SO3Proj =
    UpperSO3 Vec2
  | LowerSO3 Vec2
  deriving (Show)

-- | Lambert transformation of SO3 in R2. In spherical coordenates @(theta, omega) =>
-- (theta, 2*sin(omega/2))@ where theta is the azimuthal angle to @RD@ axis and omega is
-- the angle to @ND@ axis.
lambertSO3Proj :: Normal3 -> SO3Proj
lambertSO3Proj n
  | z < 0     = LowerSO3 $ Vec2 (-k * x) (-k * y)
  | otherwise = UpperSO3 $ Vec2 (k  * x) (k  * y)
  where
    Vec3 x y z = fromNormal n
    k          = sqrt(2 / (1 + abs z))

-- | Sterographic transformation from SO3 to R2. In spherical coordenates @(theta, omega)
-- => (theta, 2*tan(omega/2))@ where theta is the azimuthal angle to @RD@ axis and omega
-- is the angle to @ND@ axis.
steroSO3Proj :: Normal3 -> SO3Proj
steroSO3Proj n
  | z < 0     = LowerSO3 $ Vec2 (-k * x) (-k * y)
  | otherwise = UpperSO3 $ Vec2 (k  * x) (k  * y)
  where
    Vec3 x y z = fromNormal n
    k          = (1 / (1 + abs z))

getBothProj :: Vector SO3Proj -> Vector Vec2
getBothProj = let
  both (UpperSO3 x) = x
  both (LowerSO3 x) = x
  in V.map both

getLowerProj :: Vector SO3Proj -> Vector Vec2
getLowerProj = let
  isUpper (UpperSO3 _) = True
  isUpper _            = False
  in V.map (\(UpperSO3 x) -> x) . V.filter isUpper

getUpperProj :: Vector SO3Proj -> Vector Vec2
getUpperProj = let
  isUpper (UpperSO3 _) = True
  isUpper _            = False
  in V.map (\(UpperSO3 x) -> x) . V.filter isUpper

-- Get the maximum distance form origin (max. radius) of a given projection.
--getSO3ProjMaxDist :: SphereProjection -> Double
--getSO3ProjMaxDist projType = case projType of
--  Lambert      _ -> sqrt (2)
--  Sterographic _ -> 1
