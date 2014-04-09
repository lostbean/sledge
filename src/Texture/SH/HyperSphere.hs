{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE FlexibleContexts           #-}
{-# LANGUAGE MultiParamTypeClasses      #-}
{-# LANGUAGE TypeFamilies               #-}

module Texture.SH.HyperSphere
       ( -- * Unity Sphere (SO2)
         SO2 (..)
       , so2ToCart
       , cartToSO2
       , getSO2Grid
       , getSO2IsoGrid
       , so2ToSO3
         -- * Unity Hyper-Sphere (SO3)
       , SO3 (..)
       , so3ToQuaternion
       , quaternionToSO3
       , so3ToCart
       , getSO3Grid
       , getSO3IsoGrid
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

import qualified Data.List                   as L
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Generic.Base    as G
import qualified Data.Vector.Generic.Mutable as M


import           Data.Vector.Unboxed         (Vector)
import           Control.Monad               (liftM)

import           Hammer.Math.Algebra
import           Texture.Orientation
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

cartToSO2 :: Vec3 -> SO2
cartToSO2 n = let
  s = len n
  in SO2 {so2Theta = acos (_3 n / s), so2Phi = atan2 (_2 n) (_1 n)}

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
  qlm = U.fromList [ (funcSO2 np nt) | nt <- [1..nTheta-1], np <- [0..nPhi-1] ]
  q0  = SO2 0 0
  qn  = SO2 {so2Theta = pi, so2Phi = 2*pi}
  in q0 `U.cons` qlm `U.snoc` qn

getSO2Grid :: Int -> Int -> Vector SO2
getSO2Grid nTheta nPhi = let
  step_theta =     pi / (fromIntegral nTheta)
  step_phi   = 2 * pi / (fromIntegral nPhi)
  funcSO2 ip it = let
    t = step_theta * (fromIntegral it)
    p = step_phi   * (fromIntegral ip)
    in SO2 {so2Theta = t, so2Phi = p}
  qlm = U.fromList [ (funcSO2 np nt) | nt <- [1..nTheta-1], np <- [0..nPhi-1] ]
  q0  = SO2 0 0
  qn  = SO2 {so2Theta = pi, so2Phi = 2*pi}
  in q0 `U.cons` qlm `U.snoc` qn

getSO2Cells :: Int -> Int -> Vector SO2Cell
getSO2Cells nTheta nPhi = let
  ll0 = U.replicate nPhi 0
  lln = U.replicate nPhi (nPhi * (nTheta-1) + 1)
  llm = V.generate (nTheta-1) (\n -> U.enumFromN (n * nPhi + 1) nPhi)
  ll  = ll0 `V.cons` llm `V.snoc` lln
  mkStrip l1 l2 = let
    func i a b = SO2Cell (b, l2 U.! (i+1), l1 U.! (i+1), a)
    merge      = SO2Cell (U.head l2, U.head l1, U.last l1, U.last l2)
    in merge `U.cons` U.izipWith func (U.init l1) (U.init l2)
  foo acc i l = acc U.++ mkStrip l (ll V.! (i+1))
  in V.ifoldl' foo U.empty (V.init ll)

-- local instance to avoid conflict when exported.
newtype SO2Cell = SO2Cell (Int, Int, Int, Int) deriving (G.Vector U.Vector, M.MVector U.MVector, U.Unbox)

instance RenderCell SO2Cell where
  makeCell (SO2Cell (a, b, c, d)) = U.fromList [a, b, c, d]
  getType _                       = VTK_QUAD

-- ========================================= SO3 =========================================

so3ToQuaternion :: SO3 -> Quaternion
so3ToQuaternion SO3{..} = let
  q0 = cos (so3Omega / 2)
  sO = sin (so3Omega / 2)
  sT = sin so3Theta
  q1 = sO * sT * cos so3Phi
  q2 = sO * sT * sin so3Phi
  q3 = sO * cos so3Theta
  in mkQuaternion $ Vec4 q0 q1 q2 q3

quaternionToSO3 :: Quaternion -> SO3
quaternionToSO3 q = let
  (q0, qv) = splitQuaternion q
  so2   = cartToSO2 qv
  omega = 2 * acos q0
  in so2ToSO3 omega so2

so3ToCart :: SO3 -> Vec3
so3ToCart x = (so3Omega x) *& (so2ToCart $ so3ToSO2 x)

so3ToSO2 :: SO3 -> SO2
so3ToSO2 SO3{..} = SO2 {so2Theta = so3Theta, so2Phi = so3Phi}

getSO3Grid :: Int -> Int -> Int -> Vector SO3
getSO3Grid nOmega nTheta nPhi = let
  step_w = pi / (fromIntegral nOmega)
  sphere = getSO2Grid nTheta nPhi
  getW   = (step_w *) . fromIntegral
  ws     = U.generate (nOmega+1) getW
  func w = U.map (so2ToSO3 w) sphere
  in U.concatMap func ws

getSO3IsoGrid :: Int -> Int -> Int -> Vector SO3
getSO3IsoGrid nOmega nTheta nPhi = let
  step_w = 2 * pi / (fromIntegral nOmega)
  sphere = getSO2IsoGrid nTheta nPhi
  getW   = (step_w *) . fromIntegral
  foo i  = let
    w = getW i
    in (3 / 4 * (w - sin w)) ** (1 / 3)
  ws     = U.generate (nOmega+1) foo
  func w = U.map (so2ToSO3 w) sphere
  in U.concatMap func ws

getSO3Cells :: Int -> Int -> Int -> Vector SO3Cell
getSO3Cells nOmega nTheta nPhi = let
  so2cells = getSO2Cells nTheta nPhi
  layers   = U.enumFromN 0 nOmega
  -- number of points according to getSO2Grid
  psSize   = (nTheta - 1) * nPhi + 2
  func i (SO2Cell (a,b,c,d)) = let
    k1 = psSize * i
    k2 = psSize * (i + 1)
    in SO3Cell ((a+k1, b+k1, c+k1, d+k1), (a+k2, b+k2, c+k2, d+k2))
  in U.concatMap (\i -> U.map (func i) so2cells) layers

newtype SO3Cell = SO3Cell ((Int, Int, Int, Int), (Int, Int, Int, Int))
                deriving (G.Vector U.Vector, M.MVector U.MVector, U.Unbox)

instance RenderCell SO3Cell where
  makeCell (SO3Cell ((a, b, c, d), (e, f, g, h))) =
    U.fromList [a, b, c, d, e, f, g, h]
  getType _ = VTK_HEXAHEDRON

-- ==================================== Render tools =====================================

mkSO2 :: Int -> Int -> (Vector SO2, VTK Vec3)
mkSO2 nTheta nPhi = let
  cells = getSO2Cells   nTheta nPhi
  grid  = getSO2Grid nTheta nPhi
  ps    = U.map so2ToCart grid
  vtk   = mkUGVTK "Sphere" (U.convert ps) cells
  in (grid, vtk)

mkSO3 :: Int -> Int -> Int -> (Vector SO3, VTK Vec3)
mkSO3 nOmega nTheta nPhi = let
  cells = getSO3Cells   nOmega nTheta nPhi
  grid  = getSO3Grid nOmega nTheta nPhi
  ps    = U.map so3ToCart grid
  vtk   = mkUGVTK "Hyper-sphere" (U.convert ps) cells
  in (grid, vtk)

renderSO2VTK :: (SO2 -> Double) -> VTK Vec3
renderSO2VTK feval = let
  (grid, vtk) = mkSO2 60 60
  func i _    = feval (grid U.!i)
  attr        = mkPointAttr ("Intensity") func
  in addDataPoints vtk attr

renderSO3SolidVTK :: (SO3 -> Double) -> VTK Vec3
renderSO3SolidVTK feval = let
  (grid, vtk) = mkSO3 30 30 30
  func i _    = feval (grid U.!i)
  attr        = mkPointAttr ("Intensity") func
  in addDataPoints vtk attr

renderSO3ShellVTK :: (SO3 -> Double) -> VTK Vec3
renderSO3ShellVTK feval = let
  (grid, vtk) = mkSO2 30 30
  ws          = [0, 2*pi/30 .. 2*pi]
  func w i _  = feval $ so2ToSO3 w (grid U.!i)
  foo acc w   = let
    attr = mkPointAttr ("Intensity - " ++ show w) (func w)
    in addDataPoints acc attr
  in L.foldl' foo vtk ws

renderSO2PointsVTK :: Vector SO2 -> VTK Vec3
renderSO2PointsVTK lq = let
  pos  = U.map so2ToCart lq
  pids = U.enumFromN (0 :: Int) (U.length pos)
  in mkUGVTK "SO2 points" (U.convert pos) pids

renderSO3PointsVTK :: Vector SO3 -> VTK Vec3
renderSO3PointsVTK lq = let
  pos  = U.map so3ToCart lq
  pids = U.enumFromN (0 :: Int) (U.length pos)
  in mkUGVTK "SO3 points" (U.convert pos) pids

-- -------------------------------------------- Unbox SO2 ----------------------------------------------------

newtype instance U.MVector s SO2 = MV_SO2 (U.MVector s (Double, Double))
newtype instance U.Vector    SO2 = V_SO2  (U.Vector    (Double, Double))

instance U.Unbox SO2

instance M.MVector U.MVector SO2 where
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
  basicLength (MV_SO2 v)                  = M.basicLength v
  basicUnsafeSlice i n (MV_SO2 v)         = MV_SO2 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_SO2 v1) (MV_SO2 v2)   = M.basicOverlaps v1 v2
  basicUnsafeNew n                        = MV_SO2 `liftM` M.basicUnsafeNew n
  basicUnsafeReplicate n (SO2 x y)        = MV_SO2 `liftM` M.basicUnsafeReplicate n (x, y)
  basicUnsafeRead (MV_SO2 v) i            = M.basicUnsafeRead v i >>= (\(x, y) -> return $ SO2 x y)
  basicUnsafeWrite (MV_SO2 v) i (SO2 x y) = M.basicUnsafeWrite v i (x, y)
  basicClear (MV_SO2 v)                   = M.basicClear v
  basicSet (MV_SO2 v) (SO2 x y)           = M.basicSet v (x, y)
  basicUnsafeCopy (MV_SO2 v1) (MV_SO2 v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeGrow (MV_SO2 v) n            = MV_SO2 `liftM` M.basicUnsafeGrow v n

instance G.Vector U.Vector SO2 where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_SO2 v)           = V_SO2 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_SO2 v)              = MV_SO2 `liftM` G.basicUnsafeThaw v
  basicLength (V_SO2 v)                  = G.basicLength v
  basicUnsafeSlice i n (V_SO2 v)         = V_SO2 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_SO2 v) i          = G.basicUnsafeIndexM v i >>= (\(x, y) -> return $ SO2 x y)
  basicUnsafeCopy (MV_SO2 mv) (V_SO2 v ) = G.basicUnsafeCopy mv v
  elemseq _ (SO2 x y) t                  = G.elemseq (undefined :: Vector a) x $
                                           G.elemseq (undefined :: Vector a) y t

-- -------------------------------------------- Unbox SO3 ----------------------------------------------------

newtype instance U.MVector s SO3 = MV_SO3 (U.MVector s (Double, Double, Double))
newtype instance U.Vector    SO3 = V_SO3  (U.Vector    (Double, Double, Double))

instance U.Unbox SO3

instance M.MVector U.MVector SO3 where
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
  basicLength (MV_SO3 v)                    = M.basicLength v
  basicUnsafeSlice i n (MV_SO3 v)           = MV_SO3 $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_SO3 v1) (MV_SO3 v2)     = M.basicOverlaps v1 v2
  basicUnsafeNew n                          = MV_SO3 `liftM` M.basicUnsafeNew n
  basicUnsafeReplicate n (SO3 x y z)        = MV_SO3 `liftM` M.basicUnsafeReplicate n (x, y, z)
  basicUnsafeRead (MV_SO3 v) i              = M.basicUnsafeRead v i >>= (\(x,y,z) -> return $ SO3 x y z)
  basicUnsafeWrite (MV_SO3 v) i (SO3 x y z) = M.basicUnsafeWrite v i (x, y, z)
  basicClear (MV_SO3 v)                     = M.basicClear v
  basicSet (MV_SO3 v) (SO3 x y z)           = M.basicSet v (x, y, z)
  basicUnsafeCopy (MV_SO3 v1) (MV_SO3 v2)   = M.basicUnsafeCopy v1 v2
  basicUnsafeGrow (MV_SO3 v) n              = MV_SO3 `liftM` M.basicUnsafeGrow v n

instance G.Vector U.Vector SO3 where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_SO3 v)          = V_SO3 `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_SO3 v)             = MV_SO3 `liftM` G.basicUnsafeThaw v
  basicLength (V_SO3 v)                 = G.basicLength v
  basicUnsafeSlice i n (V_SO3 v)        = V_SO3 $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_SO3 v) i         = G.basicUnsafeIndexM v i >>= (\(x,y,z) -> return $ SO3 x y z)
  basicUnsafeCopy (MV_SO3 mv) (V_SO3 v) = G.basicUnsafeCopy mv v
  elemseq _ (SO3 x y z) t               = G.elemseq (undefined :: Vector a) x $
                                          G.elemseq (undefined :: Vector a) y $
                                          G.elemseq (undefined :: Vector a) z t
