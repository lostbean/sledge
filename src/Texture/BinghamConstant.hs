{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}

{- |
Solves the multivariate confluent hypergeometric function for
F = 2 * 1F1 (1/2; (d + 1) / 2; z1 .. zn) using the truncate sum. The
two variable case is one of the Humbert series.
-}

module Texture.BinghamConstant
       ( computeF
       , computedFdz1
       , computedFdz2
       , computedFdz3
       , computedFdz4
       , computeAllF
       , surface_area_sphere
       ) where

import qualified Data.Vector.Unboxed as U

import Data.Vector.Unboxed          (Vector)
import Data.Word

-- import Debug.Trace
-- dbg s x = trace (s L.++ show x) x

-- | Computation of normalization constant from the concentration values of a Bingham
-- distribution in 4 dimensions using saddle point approximations (Kume and Wood, 2005).
computeF :: Double -> Double -> Double -> Double -> Double
computeF z1 z2 z3 z4 = saddleC3 (U.fromList [-z1, -z2, -z3, -z4])

-- | Numerical calculation of 1st partial derivative of the normalization constant.
computedFdz1 :: Double -> Double -> Double -> Double -> Double
computedFdz1 z1 z2 z3 z4 = (-1) * dC (U.fromList [1,0,0,0]) (U.fromList [-z1, -z2, -z3, -z4])

-- | Numerical calculation of 2nd partial derivative of the normalization constant.
computedFdz2 :: Double -> Double -> Double -> Double -> Double
computedFdz2 z1 z2 z3 z4 = (-1) * dC (U.fromList [0,1,0,0]) (U.fromList [-z1, -z2, -z3, -z4])

-- | Numerical calculation of 3rd partial derivative of the normalization constant.
computedFdz3 :: Double -> Double -> Double -> Double -> Double
computedFdz3 z1 z2 z3 z4 = (-1) * dC (U.fromList [0,0,1,0]) (U.fromList [-z1, -z2, -z3, -z4])

-- | Numerical calculation of 4th partial derivative of the normalization constant.
computedFdz4 :: Double -> Double -> Double -> Double -> Double
computedFdz4 z1 z2 z3 z4 = (-1) * dC (U.fromList [0,0,0,1]) (U.fromList [-z1, -z2, -z3, -z4])

-- | Fast computation of normalization constant and its partial derivatives.
computeAllF :: Double -> Double -> Double -> Double -> (Double, (Double, Double, Double, Double))
computeAllF z1 z2 z3 z4 = let
  delta = 0.0001
  cs  = U.fromList [-z1, -z2, -z3, -z4]
  cs2 = U.fromList [-z1, delta - z2, -z3, -z4]
  cs3 = U.fromList [-z1, -z2, delta - z3, -z4]
  cs4 = U.fromList [-z1, -z2, -z3, delta - z4]
  x  = saddleC3 cs
  x2 = saddleC3 cs2
  x3 = saddleC3 cs3
  x4 = saddleC3 cs4
  d1 = x - d2 - d3 - d4
  d2 = -(x2-x) / delta
  d3 = -(x3-x) / delta
  d4 = -(x4-x) / delta
  in (x, (d1, d2, d3, d4))

-- ==================================== tools ============================================

-- | Computes the surface area of a unit sphere with dimension d (e.g. circle ~> d=1)
surface_area_sphere :: Int -> Double
surface_area_sphere d
  | d == 0 = 2
  | d == 1 = 2 * pi
  | d == 2 = 4 * pi
  | d == 3 = 2 * pi * pi
  | otherwise = (2 * pi / (fromIntegral $ d-1)) * surface_area_sphere (d-2)

-- ========================== Saddle point approximation =================================

factTable :: Vector Double
factTable = U.scanl' (\acc i -> acc * i) 1 (U.enumFromN 1 10000)

fact :: Int -> Double
fact = (factTable U.!)

-- | Calculates the Bingham normalization constant using saddle point approximation. This
-- method work for an arbitrary number of dimensions and is faster than truncated series
-- approximation. Further information on: Kume and Wood 2005 and Kume and Wood 2007
saddleC3 :: Vector Double -> Double
saddleC3 cs = a / b
  where
    t = findRoot cs
    n = U.length cs
    rho j = dKdt j cs t / sqrt ((dKdt 2 cs t)^j)
    prod = let
      func i = 1 / sqrt (cs U.! i - t)
      in U.product $ U.generate n func
    w = (1 / 8) * (rho 4) - (5 / 24) * (rho 3)^(2::Int)
    a = sqrt (2 * pi^(n-1)) * exp (w - t) * prod
    b = sqrt (dKdt 2 cs t)

dKdt :: Int -> Vector Double -> Double -> Double
dKdt j cs t = let
  n = U.length cs
  func i = fact (j - 1) / (cs U.! i - t)^j
  in 0.5 * (U.sum $ U.generate n func)

findRoot :: Vector Double -> Double
findRoot cs = go (lmin + delta)
  where
    delta = -0.001
    lmin  = U.minimum cs
    f  x  = (dKdt 1 cs x) - 1
    df x  = (f (x + delta) - f x) / delta
    go x
      | abs (x1 - x) < 1e-10 = x1
      | otherwise            = go x1
      where x1 = x - (f x / df x)

-- ============================== Partial derivatives ====================================

-- | Numerical calculation of partial derivatives of the normalization constant.
dC :: Vector Word -> Vector Double -> Double
dC js cs = let
  delta = 0.0001
  a = saddleC3 cs
  b = saddleC3 (U.zipWith func cs js)
  func x j = x + delta * fromIntegral j
  in (b-a)/delta

-- ================================ Testing functions ====================================

checkDerivatives :: Vector Double -> Double
checkDerivatives xs = let
  gends n = U.generate size (\i -> if i == n then 1 else 0)
  size = U.length xs
  c1 = sum [dC (gends i) xs | i <- [0 .. size-1]]
  c2 = saddleC3 xs
  in (c2 + c1) / c1
