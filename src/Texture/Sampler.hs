{-# LANGUAGE FlexibleContexts #-}

module Texture.Sampler
       ( hitAndRunSlice
       , hitAndRun
       ) where

import qualified Data.Vector.Unboxed as U

import System.Random
import System.IO

import Hammer.VTK
import Hammer.Math.Algebra

import Texture.Bingham
import Texture.Orientation
import Texture.HyperSphere

hitAndRunSlice :: (HasRandomDir a)=>
                  (a -> Double) -> a -> Double -> Int -> IO [a]
hitAndRunSlice p x0 max n = go x0 0
  where
    go xi i
      | i >= n    = return []
      | otherwise = do
        u      <- randomRIO (0, p xi)
        shoter <- linearShoter xi
        let scan = do
              k <- randomRIO (0, abs max)
              let
                xj = shoter k
                fj = p xj
              if fj >= u
                then go xj (i+1) >>= \xs -> return (xj : xs)
                else go xi i
        scan

hitAndRun :: (HasRandomDir a)=>
             (a -> Double) -> a -> Double -> Int -> IO [a]
hitAndRun p x0 max n = go x0 0
  where
    go xi i
      | i >= n    = return []
      | otherwise = do
        shoter <- linearShoter xi
        k      <- randomRIO (0, max)
        let
          xj = shoter k
          fi = p xi
          fj = p xj
        t <- randomRIO (0, 1)
        if min (fj/fi) 1 >= t
          then go xj (i+1) >>= \xs -> return (xj : xs)
          else go xi i

class HasRandomDir a where
  linearShoter :: a -> IO (Double -> a)

instance HasRandomDir Double where
  linearShoter x = do
    t <- randomRIO (0,1) >>= return . func
    return (\k -> x + k * t)
    where
      func :: Double -> Double
      func x = if x >= 0.5 then 1 else (-1)

instance HasRandomDir Vec2 where
  linearShoter x = do
    t <- sampleOne >>= func
    return (\k -> x &+ k *& t)
    where
      func v
        | l > 1     = sampleOne >>= func
        | l == 0    = return v
        | otherwise = return (v &* (1/l))
        where l = norm v
      sampleOne = do
        x <- randomRIO (-1,1)
        y <- randomRIO (-1,1)
        return (Vec2 x y)

instance HasRandomDir Quaternion where
  linearShoter x = do
    t <- sampleOne >>= func
    return (\k -> x #<= (toQuaternion $ mkAxisPair t (Rad k)))
    where
      func v
        | l > 1     = sampleOne >>= func
        | l == 0    = return v
        | otherwise = return (v &* (1/l))
        where l = norm v
      sampleOne = do
        x <- randomRIO (-1,1)
        y <- randomRIO (-1,1)
        z <- randomRIO (-1,1)
        return (Vec3 x y z)

-- ===================================== Testing =========================================

testFunc = let
  f1 = multiNormalPDF (Mat2 (Vec2 0.1 0.1) (Vec2 0.1 0.2)) (Vec2 3 3)
  f2 = multiNormalPDF (Mat2 (Vec2 1 0.5) (Vec2 0.5 1)) (Vec2 (-5) (-5))
  in \x -> f1 x + f2 x

viewFunc = do
  h <- openFile "data-func" WriteMode
  let xs = [(Vec2 x y) | x <- [-20, -19.5 .. 20], y <- [-20, -19.5 .. 20]]
  mapM_ (\v@(Vec2 x y) -> hPutStrLn h $ show x ++ "  " ++ show y ++ " " ++ show (testFunc v)) xs
  hClose h

testSampler n = do
  h <- openFile "data-samp" WriteMode
  xs <- hitAndRunSlice testFunc (Vec2 0 0) 20 n
  mapM_ (\(Vec2 x y) -> hPutStrLn h $ show x ++ "  " ++ show y) xs
  hClose h

testSamplerQuality n = let
  da1 = (1, mkQuaternion  (Vec4 0 0 1 0))
  da2 = (10, mkQuaternion  (Vec4 0 1 0 0))
  da3 = (10, mkQuaternion (Vec4 1 0 0 1))
  din = mkBingham da1 da2 da3
  in do
    xs <- hitAndRunSlice (binghamPDF din) (zerorot) (2*pi) n
    writeQuater "Bing-PDF-In-testSamplerQuality" $ renderBingham din
    writeQuater "Bing-Samples-testSamplerQuality" $ renderPoints xs
    let dout = fitBingham (U.fromList xs)
    writeQuater "Bing-PDF-Out-testSamplerQuality" $ renderBingham dout

testSamplerMultiModal n = let
  da1 = (1, mkQuaternion  (Vec4 0 0 1 0))
  da2 = (1, mkQuaternion  (Vec4 0 1 0 0))
  da3 = (20, mkQuaternion (Vec4 1 0 0 1))
  da  = mkBingham da1 da2 da3
  db1 = (1, mkQuaternion  (Vec4 1 0 0 (-1)))
  db2 = (20, mkQuaternion (Vec4 0 1 0 0))
  db3 = (1, mkQuaternion  (Vec4 1 0 0 1))
  db  = mkBingham db1 db2 db3
  in do
    xs <- hitAndRunSlice (\q -> binghamPDF da q + binghamPDF db q) (zerorot) (2*pi) n
    writeQuater "Bing-PDF-A-testSamplerMultiModal" $ renderBingham da
    writeQuater "Bing-PDF-B-testSamplerMultiModal" $ renderBingham db
    writeQuater "Bing-Samples-testSamplerMultiModal" $ renderPoints xs


renderPoints :: [Quaternion] -> VTK Vec3
renderPoints = renderSO3PointsVTK . U.map (quaternionToSO3) . U.fromList

writeQuater :: (RenderElemVTK a)=> String -> VTK a -> IO ()
writeQuater name = writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vtu") True

-- ==================================== n-Dim Normal =====================================

-- | Inverse of error function.
invERF :: Double -> Double
invERF x
  | x < 0     = -invERF (-x)
  | otherwise = let
    a = 0.147
    y1 = (2 / (pi * a) + log (1 - x * x) / 2)
    y2 = sqrt (y1 * y1 - ( 1 / a) * log (1 - x * x))
    in sqrt (y2 - y1)

-- | Compute the normal PDF.
normalPDF :: Double -> Double -> Double -> Double
normalPDF x mu sigma = let
  dx = x - mu
  s2 = 2 * sigma * sigma
  in exp(-(dx * dx) / s2) / (sqrt (2 * pi) * sigma)

-- | Generate a random sample from a univariate normal distribution.
-- The random input must range [0, 1]
normalSample :: Double -> Double -> Double -> Double
normalSample mu sigma rnd = mu + sigma * sqrt 2 * invERF (2 * rnd - 1)

-- | Compute a multivariate normal pdf in principal components form.
-- Uses the eigenvectors and eigenvalues to calculate the inverse of
-- the covariance matrix.
multiNormalPDF :: Mat2 -> Vec2 -> Vec2 -> Double
multiNormalPDF cv mu x = let
  d  = 4
  dx = x &- mu
  k  = ((2 * pi) ** (d/2))
  -- Is the same as: SUM ( <dx, V[i]> / z[i] )^2
  e  = ((dx .* (transpose cv)) &. dx)
  in (exp $ (-0.5) * e) / k
