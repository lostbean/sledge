{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE FlexibleContexts           #-}
{-# LANGUAGE MultiParamTypeClasses      #-}
{-# LANGUAGE TypeFamilies               #-}

module Hammer.Texture.SH.HyperSphere
       ( -- * Unity Sphere (SO2)
         SO2 (..)
       , so2ToCart
       , getSO2Grid
       , so2ToSO3
         -- * Unity Hyper-Sphere (SO3)
       , SO3 (..)
       , so3ToQuaternion
       , quaternionToSO3
       , so3ToCart
       , getSO3Grid
       , so3ToSO2
         -- * VTK Ploting
       , mkSO2
       , renderSO2VTK
       , renderSO2PointsVTK
       , mkSO3
       , renderSO3ShellVTK
       , renderSO3SolidVTK
       , renderSO3PointsVTK
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Data.Vector         (Vector)

import           Hammer.Math.Algebra
import           Hammer.Texture.Orientation
import           Hammer.Render.VTK.VTKRender

import           Debug.Trace
dbg s x = trace (s L.++ show x) x


data SO2 =
  SO2
  { so2Theta :: Double
  , so2Phi   :: Double
  } deriving (Show)

data SO3 =
  SO3
  { so3Omega :: Double
  , so3Theta :: Double
  , so3Phi   :: Double
  } deriving (Show)

cartToSO2 :: Double -> Double -> SO2
cartToSO2 x1 x2 = let
  x = 2 * x1 * sqrt (1 - x1*x1 - x2*x2)
  y = 2 * x2 * sqrt (1 - x1*x1 - x2*x2)
  z = 1 - 2 * (x1*x1 + x2*x2)
  in SO2 {so2Theta = acos z, so2Phi = atan (y / x)}

-- ========================================= SO2 =========================================

so2ToCart :: SO2 -> Vec3
so2ToCart SO2{..} = let
  x = sin so2Theta * cos so2Phi
  y = sin so2Theta * sin so2Phi
  z = cos so2Theta
  in Vec3 x y z

so2ToSO3 :: Double -> SO2 -> SO3
so2ToSO3 x SO2{..} = SO3 {so3Omega = x, so3Theta = so2Theta, so3Phi = so2Phi}

getSO2IsoGrid :: Int -> Int -> Vector SO2
getSO2IsoGrid nTheta nPhi = let
  step_theta =     1  / (fromIntegral nTheta)
  step_phi   = 2 * pi / (fromIntegral nPhi)
  funcSO2 ip it = let
    t = acos $ 1 - 2 * step_theta * (fromIntegral it)
    p = step_phi * (fromIntegral ip)
    in SO2 {so2Theta = t, so2Phi = p}
  qlm = V.fromList [ (funcSO2 np nt) | nt <- [1..nTheta-1], np <- [0..nPhi-1] ]
  q0  = SO2 0 0
  qn  = SO2 {so2Theta = pi, so2Phi = 2*pi}
  in q0 `V.cons` qlm `V.snoc` qn

getSO2Grid :: Int -> Int -> Vector SO2
getSO2Grid nTheta nPhi = let
  step_theta =     pi / (fromIntegral nTheta)
  step_phi   = 2 * pi / (fromIntegral nPhi)
  funcSO2 ip it = let
    t = step_theta * (fromIntegral it)
    p = step_phi   * (fromIntegral ip)
    in SO2 {so2Theta = t, so2Phi = p}
  qlm = V.fromList [ (funcSO2 np nt) | nt <- [1..nTheta-1], np <- [0..nPhi-1] ]
  q0  = SO2 0 0
  qn  = SO2 {so2Theta = pi, so2Phi = 2*pi}
  in q0 `V.cons` qlm `V.snoc` qn

getSO2Cells :: Int -> Int -> Vector SO2Cell
getSO2Cells nTheta nPhi = let
  ll0 = V.replicate nPhi 0
  lln = V.replicate nPhi (nPhi * (nTheta-1) + 1)
  llm = V.generate (nTheta-1) (\n -> V.enumFromN (n * nPhi + 1) nPhi)
  ll  = ll0 `V.cons` llm `V.snoc` lln
  mkStrip l1 l2 = let
    func i a b = SO2Cell (b, l2 V.! (i+1), l1 V.! (i+1), a)
    merge      = SO2Cell (V.head l2, V.head l1, V.last l1, V.last l2)
    in merge `V.cons` V.izipWith func (V.init l1) (V.init l2)
  foo acc i l = acc V.++ mkStrip l (ll V.! (i+1))
  in V.ifoldl' foo V.empty (V.init ll)

-- local instance to avoid conflict when exported.
newtype SO2Cell = SO2Cell (Int, Int, Int, Int)

instance RenderCell SO2Cell where
  makeCell (SO2Cell (a, b, c, d)) = U.fromList [a, b, c, d]
  getType _                       = VTK_QUAD

-- ========================================= SO3 =========================================

so3ToQuaternion :: SO3 -> Quaternion
so3ToQuaternion SO3{..} = let
  q0 = cos so3Omega
  sT = sin so3Theta
  sO = sin so3Omega
  q1 = sO * sT * cos so3Phi
  q2 = sO * sT * sin so3Phi
  q3 = sO * cos so3Theta
  in mkQuaternion $ Vec4 q0 q1 q2 q3

quaternionToSO3 :: Quaternion -> SO3
quaternionToSO3 q = let
  (q0, Vec3 q1 _ q3) = splitQuaternion q
  omega = acos q0
  theta = acos $ q3 / sin omega
  phi   = acos $ q1 / (sin omega * sin theta)
  in SO3 {so3Omega = omega, so3Theta = theta, so3Phi = phi}

so3ToCart :: SO3 -> Vec3
so3ToCart x = (so3Omega x) *& (so2ToCart $ so3ToSO2 x)

so3ToSO2 :: SO3 -> SO2
so3ToSO2 SO3{..} = SO2 {so2Theta = so3Theta, so2Phi = so3Phi}

getSO3Grid :: Int -> Int -> Int -> Vector SO3
getSO3Grid nOmega nTheta nPhi = let
  step_w = pi / (fromIntegral nOmega)
  sphere = getSO2Grid nTheta nPhi
  getW   = (step_w *) . fromIntegral
  ws     = V.generate (nOmega+1) getW
  func w = V.map (so2ToSO3 w) sphere
  in V.concatMap func ws

getSO3IsoGrid :: Int -> Int -> Int -> Vector SO3
getSO3IsoGrid nOmega nTheta nPhi = let
  step_w = 2 * pi / (fromIntegral nOmega)
  sphere = getSO2IsoGrid nTheta nPhi
  getW   = (step_w *) . fromIntegral
  foo i  = let
    w = getW i
    in (3 / 4 * (w - sin w)) ** (1 / 3)
  ws     = V.generate (nOmega+1) foo
  func w = V.map (so2ToSO3 w) sphere
  in V.concatMap func ws

getSO3Cells :: Int -> Int -> Int -> Vector SO3Cell
getSO3Cells nOmega nTheta nPhi = let
  so2cells = getSO2Cells nTheta nPhi
  layers   = V.enumFromN 0 nOmega
  -- number of points according to getSO2Grid
  psSize   = (nTheta - 1) * nPhi + 2
  func i (SO2Cell (a,b,c,d)) = let
    k1 = psSize * i
    k2 = psSize * (i + 1)
    in SO3Cell (a+k1, b+k1, c+k1, d+k1, a+k2, b+k2, c+k2, d+k2)
  in V.concatMap (\i -> V.map (func i) so2cells) layers

newtype SO3Cell = SO3Cell (Int, Int, Int, Int, Int, Int, Int, Int)

instance RenderCell SO3Cell where
  makeCell (SO3Cell (a, b, c, d, e, f, g, h)) =
    U.fromList [a, b, c, d, e, f, g, h]
  getType _ = VTK_HEXAHEDRON

-- ==================================== Render tools =====================================

mkSO2 :: Int -> Int -> (Vector SO2, VTK Vec3)
mkSO2 nTheta nPhi = let
  cells = getSO2Cells   nTheta nPhi
  grid  = getSO2IsoGrid nTheta nPhi
  ps    = V.map so2ToCart grid
  vtk   = mkUGVTK "Sphere" (V.convert ps) cells
  in (grid, vtk)

mkSO3 :: Int -> Int -> Int -> (Vector SO3, VTK Vec3)
mkSO3 nOmega nTheta nPhi = let
  cells = getSO3Cells   nOmega nTheta nPhi
  grid  = getSO3IsoGrid nOmega nTheta nPhi
  ps    = V.map so3ToCart grid
  vtk   = mkUGVTK "Hyper-sphere" (V.convert ps) cells
  in (grid, vtk)

renderSO2VTK :: (SO2 -> Double) -> VTK Vec3
renderSO2VTK feval = let
  (grid, vtk) = mkSO2 60 60
  func i _    = feval (grid V.!i)
  attr        = mkPointAttr ("Intensity") func
  in addDataPoints vtk attr

renderSO3SolidVTK :: (SO3 -> Double) -> VTK Vec3
renderSO3SolidVTK feval = let
  (grid, vtk) = mkSO3 30 30 30
  func i _    = feval (grid V.!i)
  attr        = mkPointAttr ("Intensity") func
  in addDataPoints vtk attr

renderSO3ShellVTK :: (SO3 -> Double) -> VTK Vec3
renderSO3ShellVTK feval = let
  (grid, vtk) = mkSO2 30 30
  ws          = [0, 2*pi/30 .. 2*pi]
  func w i _  = feval $ so2ToSO3 w (grid V.!i)
  foo acc w   = let
    attr = mkPointAttr ("Intensity - " ++ show w) (func w)
    in addDataPoints acc attr
  in L.foldl' foo vtk ws

renderSO2PointsVTK :: Vector SO2 -> VTK Vec3
renderSO2PointsVTK lq = let
  pos  = V.map so2ToCart lq
  pids = V.enumFromN (0 :: Int) (V.length pos)
  in mkUGVTK "SO2 points" (V.convert pos) pids

renderSO3PointsVTK :: Vector SO3 -> VTK Vec3
renderSO3PointsVTK lq = let
  pos  = V.map so3ToCart lq
  pids = V.enumFromN (0 :: Int) (V.length pos)
  in mkUGVTK "SO3 points" (V.convert pos) pids
