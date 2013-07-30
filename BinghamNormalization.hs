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

module Hammer.Texture.BinghamNormalization
       ( computeF
       , surface_area_sphere
       , computedFdz1
       , computedFdz2
       , computedFdz3
       ) where

import qualified Data.List           as L
import qualified Data.Vector.Unboxed as U

import Data.Vector.Unboxed          (Vector)

import Math.Gamma

--import Debug.Trace
--dbg s x = trace (s L.++ show x) x

epsilon :: Double
epsilon = 1e-8

iteration_factor :: Int
iteration_factor = 10

min_iterations :: Int
min_iterations = 10

-- ===================================== 1F1 =============================================

-- | Computes the hypergeometric function 1F1(a;b;z1,z2,z3) with a = 1/2, b = (dim+1)/2
computeF :: Int -> Double -> Double -> Double -> Double
computeF dime z1 z2 z3 = case L.sort $ L.filter ((>epsilon).abs) [z1, z2, z3] of
  [s1]
    | s1 > 0    -> rp * rp * compute_1F1_1d a b s1 iter
    | otherwise -> exp s1 * computeF dime (-s1) 0 0
  [s1, s2]
    | s1 > 0    -> rp * compute_1F1_2d a b s1 s2 iter
    | otherwise -> exp s1 * computeF dime (-s1) (s2-s1) 0
  [s1, s2, s3]
    | s1 > 0    -> compute_1F1_3d a b s1 s2 s3 iter
    | otherwise -> exp s1 * computeF dime (-s1) (s2-s1) (s3-s1)
  _ -> surface_area_sphere dime   -- Uniform
  where
    rp = sqrt pi -- equal to gamma(1/2)
    a  = 0.5
    b  = 0.5 * (fromIntegral $ dime + 1)
    -- number of iteractions
    zmax = floor $ max (max (abs z1) (abs z2)) (abs z3)
    iter = max (zmax * iteration_factor) min_iterations

-- | Computes the hypergeometric function 1F1(a;b;z1,z2,z3) in canonical form (z1 > z2 > z3 > 0)
compute_1F1_3d :: Double -> Double -> Double -> Double -> Double -> Int -> Double
compute_1F1_3d a b z1 z2 z3 iter = let
  func :: Double -> Int -> Int -> Int -> Double
  func !acc !i !j !k
    | i > iter  = acc
    | j > iter  = func acc (i+1) 0 0
    | k > iter  = func acc i (j+1) 0
    | isToSmall = func acc i (j+1) 0
    | otherwise = func (acc + exp g) i j (k+1)
    where
      isToSmall = (di > z1 || dj > z2 || dk > z3) && (exp g < epsilon * acc)
      di = fromIntegral i
      dj = fromIntegral j
      dk = fromIntegral k
      g = lnGamma (di + a) + lnGamma (dj + a) + lnGamma (dk + a) - lnGamma (di + dj + dk + b)
          + di * log z1 + dj * log z2 + dk * log z3
          - logFact i - logFact j - logFact k
  in 2 * sqrt(pi) * (func 0 0 0 0)

-- | Computes the hypergeometric function 1F1(a;b;z1,z2) in canonical form (z1 > z2 > 0)
compute_1F1_2d :: Double -> Double -> Double -> Double -> Int -> Double
compute_1F1_2d a b z1 z2 iter = let
  func :: Double -> Int -> Int -> Double
  func !acc !i !j
    | i > iter  = acc
    | j > iter  = func acc (i+1) 0
    | isToSmall = func acc (i+1) 0
    | otherwise = func (acc + exp g) i (j+1)
    where
      isToSmall = (di > z1 || dj > z2) && (exp g < epsilon * acc)
      di = fromIntegral i
      dj = fromIntegral j
      g = lnGamma (di + a) + lnGamma (dj + a) - lnGamma (di + dj + b) +
          di * log z1 + dj * log z2 - logFact i - logFact j
  in 2 * sqrt(pi) * (func 0 0 0)

-- | Computes the hypergeometric function 1F1(a;b;z1) in canonical form (z1 > 0)
compute_1F1_1d :: Double -> Double -> Double -> Int -> Double
compute_1F1_1d a b z1 iter = let
  func :: Double -> [Int] -> Double
  func acc [] = acc
  func acc (i:s)
    | (di > z1) && (exp g < epsilon * acc) = acc
    | otherwise = func (acc + exp g) s
    where
      di = fromIntegral i
      g = lnGamma (di + a) - lnGamma (di + b) + di * (log z1) - logFact i
  in 2 * sqrt(pi) * (func 0 [0..iter])

-- ==================================== d1F1/dz ==========================================

-- | Computes the hypergeometric function 1F1(a;b;z1,z2,z3) with a = 1/2, b = (dim+1)/2
computedFdz1 :: Int -> Double -> Double -> Double -> Double
computedFdz1 dime z1 z2 z3 = case L.sort $ L.filter ((>epsilon).abs) [z1, z2, z3] of
  [s1]
    | s1 > 0    -> rp * rp * compute_d1F1_dz_1d a b s1 iter
    | otherwise -> computeF dime s1 0 0 - exp s1 * computedFdz1 dime (-s1) 0 0
  [s1, s2]
    | s1 > 0    -> rp * compute_d1F1_dz1_2d a b s1 s2 iter
    | otherwise -> computeF dime s1 s2 0 - exp s1 *
                   ( computedFdz1 dime (-s1) (s2-s1) 0 +
                     computedFdz2 dime (-s1) (s2-s1) 0 )
  [s1, s2, s3]
    | s1 > 0    -> compute_d1F1_dz1_3d a b s1 s2 s3 iter
    | otherwise -> computeF dime s1 s2 s3 -
                   exp s1 * ( computedFdz1 dime (-s1) (s2-s1) (s3-s1) +
                              computedFdz2 dime (-s1) (s2-s1) (s3-s1) +
                              computedFdz3 dime (-s1) (s2-s1) (s3-s1) )
  _ -> (surface_area_sphere dime) / (fromIntegral $ dime + 1)   -- Uniform
  where
    rp = sqrt pi
    a  = 0.5
    b  = 0.5 * (fromIntegral $ dime + 1)
    -- number of iteractions
    zmax = floor $ max (max (abs z1) (abs z2)) (abs z3)
    iter = max (zmax * iteration_factor) min_iterations

-- | Computes the hypergeometric function 1F1(a;b;z1,z2,z3) with a = 1/2, b = (dim+1)/2
computedFdz2 :: Int -> Double -> Double -> Double -> Double
computedFdz2 dime z1 z2 z3 = case L.sort $ L.filter ((>epsilon).abs) [z1, z2, z3] of
  [s1]
    | s1 > 0    -> 0.5 * rp * rp * compute_1F1_1d a b2 s1 iter  -- (dim+2)
    | otherwise -> 0.5 * exp s1 * computeF (dime+2) (-s1) 0 0   -- (dim+2)
  [s1, s2]
    | s1 > 0    -> rp * compute_d1F1_dz2_2d a b s1 s2 iter
    | otherwise -> exp s1 * computedFdz2 dime (-s1) (s2-s1) 0
  [s1, s2, s3]
    | s1 > 0    -> compute_d1F1_dz2_3d a b s1 s2 s3 iter
    | otherwise -> exp s1 * computedFdz3 dime (-s1) (s2-s1) (s3-s1)
  _ -> (surface_area_sphere dime) / (fromIntegral $ dime + 1)   -- Uniform
  where
    rp = sqrt pi
    a  = 0.5
    b  = 0.5 * (fromIntegral $ dime + 1)
    b2 = 0.5 * (fromIntegral $ dime + 3)
    -- number of iteractions
    zmax = floor $ max (max (abs z1) (abs z2)) (abs z3)
    iter = max (zmax * iteration_factor) min_iterations

-- | Computes the hypergeometric function 1F1(a;b;z1,z2,z3) with a = 1/2, b = (dim+1)/2
computedFdz3 :: Int -> Double -> Double -> Double -> Double
computedFdz3 dime z1 z2 z3 = case L.sort $ L.filter ((>epsilon).abs) [z1, z2, z3] of
  [s1]
    | s1 > 0    -> 0.5 * rp * rp * compute_1F1_1d a b2 s1 iter      -- (dim+2)
    | otherwise -> 0.5 * exp s1 * computeF (dime+2) (-s1) 0 0       -- (dim+2)
  [s1, s2]
    | s1 > 0    -> 0.5 * rp * compute_1F1_2d a b2 s1 s2 iter        -- (dim+2)
    | otherwise -> 0.5 * exp s1 * computeF (dime+2) (-s1) (s2-s1) 0 -- (dim+2)
  [s1, s2, s3]
    | s1 > 0    -> compute_d1F1_dz3_3d a b s1 s2 s3 iter
    | otherwise -> exp s1 * computedFdz2 dime (-s1) (s2-s1) (s3-s1)
  _ -> (surface_area_sphere dime) / (fromIntegral $ dime + 1)       -- Uniform
  where
    rp = sqrt pi
    a  = 0.5
    b  = 0.5 * (fromIntegral $ dime + 1)
    b2 = 0.5 * (fromIntegral $ dime + 3)
    -- number of iteractions
    zmax = floor $ max (max (abs z1) (abs z2)) (abs z3)
    iter = max (zmax * iteration_factor) min_iterations

-- ----------------------------------- d1F1/dzx 3D ---------------------------------------

-- | Computes the hypergeometric function 1F1(a;b;z1,z2,z3) in canonical form (z1 > z2 > z3 > 0)
compute_d1F1_dz1_3d :: Double -> Double -> Double -> Double -> Double -> Int -> Double
compute_d1F1_dz1_3d a b z1 z2 z3 iter = let
  func :: Double -> Int -> Int -> Int -> Double
  func !acc !i !j !k
    | i > iter  = acc
    | j > iter  = func acc (i+1) 0 0
    | k > iter  = func acc i (j+1) 0
    | isToSmall = func acc i (j+1) 0
    | otherwise = func (acc + exp g) i j (k+1)
    where
      isToSmall = (di > z1 || dj > z2 || dk > z3) && (exp g < epsilon * acc)
      di = fromIntegral i
      dj = fromIntegral j
      dk = fromIntegral k
      g = lnGamma (di + a) + lnGamma (dj + a) + lnGamma (dk + a) - lnGamma (di + dj + dk + b)
          + (di-1) * log z1 + dj * log z2 + dk * log z3
          - logFact (i-1) - logFact j - logFact k
  in 2 * sqrt(pi) * (func 0 1 0 0)

-- | Computes the hypergeometric function 1F1(a;b;z1,z2,z3) in canonical form (z1 > z2 > z3 > 0)
compute_d1F1_dz2_3d :: Double -> Double -> Double -> Double -> Double -> Int -> Double
compute_d1F1_dz2_3d a b z1 z2 z3 iter = let
  func :: Double -> Int -> Int -> Int -> Double
  func !acc !i !j !k
    | i > iter  = acc
    | j > iter  = func acc (i+1) 1 0
    | k > iter  = func acc i (j+1) 0
    | isToSmall = func acc i (j+1) 0
    | otherwise = func (acc + exp g) i j (k+1)
    where
      isToSmall = (di > z1 || dj > z2 || dk > z3) && (exp g < epsilon * acc)
      di = fromIntegral i
      dj = fromIntegral j
      dk = fromIntegral k
      g = lnGamma (di + a) + lnGamma (dj + a) + lnGamma (dk + a) - lnGamma (di + dj + dk + b)
          + di*log z1 + (dj-1)*log z2 + dk*log z3
          - logFact i - logFact (j-1) - logFact k
  in 2 * sqrt(pi) * (func 0 0 1 0)

-- | Computes the hypergeometric function 1F1(a;b;z1,z2,z3) in canonical form (z1 > z2 > z3 > 0)
compute_d1F1_dz3_3d :: Double -> Double -> Double -> Double -> Double -> Int -> Double
compute_d1F1_dz3_3d a b z1 z2 z3 iter = let
  func :: Double -> Int -> Int -> Int -> Double
  func !acc !i !j !k
    | i > iter  = acc
    | j > iter  = func acc (i+1) 0 1
    | k > iter  = func acc i (j+1) 1
    | isToSmall = func acc i (j+1) 1
    | otherwise = func (acc + exp g) i j (k+1)
    where
      isToSmall = (di > z1 || dj > z2 || dk > z3) && (exp g < epsilon * acc)
      di = fromIntegral i
      dj = fromIntegral j
      dk = fromIntegral k
      g = lnGamma (di + a) + lnGamma (dj + a) + lnGamma (dk + a) - lnGamma (di + dj + dk + b)
          + di * log z1 + dj * log z2 + (dk-1) * log z3
          - logFact i - logFact j - logFact (k-1)
  in 2 * sqrt(pi) * (func 0 0 0 1)

-- ---------------------------------- d1F1/dz 2D -----------------------------------------

-- | Computes the hypergeometric function 1F1(a;b;z1,z2) in canonical form (z1 > z2 > 0)
compute_d1F1_dz1_2d :: Double -> Double -> Double -> Double -> Int -> Double
compute_d1F1_dz1_2d a b z1 z2 iter = let
  func :: Double -> Int -> Int -> Double
  func !acc !i !j
    | i > iter  = acc
    | j > iter  = func acc (i+1) 0
    | isToSmall = func acc (i+1) 0
    | otherwise = func (acc + exp g) i (j+1)
    where
      isToSmall = (di > z1 || dj > z2) && (exp g < epsilon * acc)
      di = fromIntegral i
      dj = fromIntegral j
      g = lnGamma (di + a) + lnGamma (dj + a) - lnGamma (di + dj + b) +
          (di-1) * log z1 + dj * log z2 - logFact (i-1) - logFact j
  in 2 * sqrt(pi) * (func 0 1 0)

-- | Computes the hypergeometric function 1F1(a;b;z1,z2) in canonical form (z1 > z2 > 0)
compute_d1F1_dz2_2d :: Double -> Double -> Double -> Double -> Int -> Double
compute_d1F1_dz2_2d a b z1 z2 iter = let
  func :: Double -> Int -> Int -> Double
  func !acc !i !j
    | i > iter  = acc
    | j > iter  = func acc (i+1) 1
    | isToSmall = func acc (i+1) 1
    | otherwise = func (acc + exp g) i (j+1)
    where
      isToSmall = (di > z1 || dj > z2) && (exp g < epsilon * acc)
      di = fromIntegral i
      dj = fromIntegral j
      g = lnGamma (di + a) + lnGamma (dj + a) - lnGamma (di + dj + b) +
          di * log z1 + (dj-1) * log z2 - logFact i - logFact (j-1)
  in 2 * sqrt(pi) * (func 0 0 1)

-- ------------------------------------ d1F1/dz 1D ---------------------------------------

-- | Computes the hypergeometric function 1F1(a;b;z1) in canonical form (z1 > 0)
compute_d1F1_dz_1d :: Double -> Double -> Double -> Int -> Double
compute_d1F1_dz_1d a b z1 iter = let
  func :: Double -> [Int] -> Double
  func acc [] = acc
  func acc (i:s)
    | (di > z1) && (exp g < epsilon * acc) = acc
    | otherwise = func (acc + exp g) s
    where
      di = fromIntegral i
      g = lnGamma (di + a) - lnGamma (di + b) + (di-1) * (log z1) - logFact (i-1)
  in 2 * sqrt(pi) * (func 0 [1..iter])

-- ==================================== tools ============================================

-- | Computes the surface area of a unit sphere with dimension d
surface_area_sphere :: Int -> Double
surface_area_sphere d
  | d == 0 = 2
  | d == 1 = 2 * pi
  | d == 2 = 4 * pi
  | d == 3 = 2 * pi * pi
  | otherwise = (2 * pi / (fromIntegral $ d-1)) * surface_area_sphere (d-2)

logFactTable :: Vector Double
logFactTable = U.scanl' (\acc i -> acc + log i) 0 (U.enumFromN 1 10000)

logFact :: Int -> Double
logFact = (logFactTable U.!)
