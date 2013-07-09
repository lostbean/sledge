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
       ( SO3Proj
       , lambertSO3Proj
       , steroSO3Proj
       , so3ProjCoord
       , isUpperSO3
       , isLowerSO3
       , getLowerProj
       , getUpperProj
       , getBothProj
       ) where

import qualified Data.Vector as V

import           Data.Vector (Vector)

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
  | z < 0     = LowerSO3 $ Vec2  (k * x)  (k * y)
  | otherwise = UpperSO3 $ Vec2 (-k * x) (-k * y)
  where
    Vec3 x y z = fromNormal n
    -- max. radius corrected to 1
    k          = sqrt(1 / (1 + abs z))

-- | Sterographic transformation from SO3 to R2. In spherical coordenates @(theta, omega)
-- => (theta, 2*tan(omega/2))@ where theta is the azimuthal angle to @RD@ axis and omega
-- is the angle to @ND@ axis.
steroSO3Proj :: Normal3 -> SO3Proj
steroSO3Proj n
  | z < 0     = LowerSO3 $ Vec2  (k * x)  (k * y)
  | otherwise = UpperSO3 $ Vec2 (-k * x) (-k * y)
  where
    Vec3 x y z = fromNormal n
    k          = (1 / (1 + abs z))

so3ProjCoord :: SO3Proj -> Vec2
so3ProjCoord (UpperSO3 x) = x
so3ProjCoord (LowerSO3 x) = x

isUpperSO3 :: SO3Proj -> Bool
isUpperSO3 (UpperSO3 _) = True
isUpperSO3 _            = False

isLowerSO3 :: SO3Proj -> Bool
isLowerSO3 = not . isUpperSO3

getBothProj :: Vector SO3Proj -> Vector Vec2
getBothProj = V.map so3ProjCoord

getLowerProj :: Vector SO3Proj -> Vector Vec2
getLowerProj = V.map so3ProjCoord . V.filter isLowerSO3

getUpperProj :: Vector SO3Proj -> Vector Vec2
getUpperProj = V.map so3ProjCoord . V.filter isUpperSO3
