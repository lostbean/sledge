module Main where

import qualified Data.Vector.Unboxed as U

import System.IO

import Hammer.VTK
import Linear.Vect

import Texture.Bingham
import Texture.Orientation
import Texture.HyperSphere
import Texture.Sampler

main :: IO ()
main = do
  testSampFit   1000
  testSampMulti 1000

  -- Simple test in R2
  viewFuncR2          -- plot multi-modal target function
  testSamplerR2 1000  -- plot samples

-- ==================================== n-Dim Normal =====================================

testSampFit :: Int -> IO ()
testSampFit n = let
  da1 = (1,  mkQuaternion (Vec4 0 0 1 0))
  da2 = (10, mkQuaternion (Vec4 0 1 0 0))
  da3 = (10, mkQuaternion (Vec4 1 0 0 1))
  din = mkBingham da1 da2 da3
  in do
    xs <- hitAndRunSlice defaultCfg (binghamPDF din) (zerorot) n
    let dout = fitBingham (U.fromList xs)
    writeQuater "Bing-PDF-In-testSamplerQuality"  $ renderBingham din
    writeQuater "Bing-Samples-testSamplerQuality" $ renderPoints  xs
    writeQuater "Bing-PDF-Out-testSamplerQuality" $ renderBingham dout

testSampMulti :: Int -> IO ()
testSampMulti n = let
  da1 = (1,  mkQuaternion (Vec4 0 0 1 0))
  da2 = (1,  mkQuaternion (Vec4 0 1 0 0))
  da3 = (20, mkQuaternion (Vec4 1 0 0 1))
  da  = mkBingham da1 da2 da3
  db1 = (1,  mkQuaternion (Vec4 1 0 0 (-1)))
  db2 = (20, mkQuaternion (Vec4 0 1 0 0))
  db3 = (1,  mkQuaternion (Vec4 1 0 0 1))
  db  = mkBingham db1 db2 db3
  in do
    xs <- hitAndRunSlice defaultCfg (\q -> binghamPDF da q + binghamPDF db q) (zerorot) n
    writeQuater "Bing-PDF-A-testSamplerMultiModal"   $ renderBingham da
    writeQuater "Bing-PDF-B-testSamplerMultiModal"   $ renderBingham db
    writeQuater "Bing-Samples-testSamplerMultiModal" $ renderPoints  xs

renderPoints :: [Quaternion] -> VTK Vec3
renderPoints = renderSO3PointsVTK . U.map (quaternionToSO3) . U.fromList

writeQuater :: (RenderElemVTK a)=> String -> VTK a -> IO ()
writeQuater n = writeUniVTKfile ("/home/edgar/Desktop/" ++ n ++ ".vtu") True

-- ============================= Testing R2 (Euclidian) ==================================

testFunc :: Vec2 -> Double
testFunc = let
  f1 = multiNormalPDF (Mat2 (Vec2 0.1 0.1) (Vec2 0.1 0.2)) (Vec2 3      3 )
  f2 = multiNormalPDF (Mat2 (Vec2 1   0.5) (Vec2 0.5   1)) (Vec2 (-5) (-5))
  in \x -> f1 x + f2 x

viewFuncR2 :: IO ()
viewFuncR2 = do
  h <- openFile "data-func" WriteMode
  let xs = [(Vec2 x y) | x <- [-20, -19.5 .. 20], y <- [-20, -19.5 .. 20]]
  mapM_ (\v@(Vec2 x y) -> hPutStrLn h $ show x ++ "  " ++ show y ++ " " ++ show (testFunc v)) xs
  hClose h

testSamplerR2 :: Int -> IO ()
testSamplerR2 n = do
  h <- openFile "data-samp" WriteMode
  xs <- hitAndRunSlice defaultCfg testFunc (Vec2 0 0) n
  mapM_ (\(Vec2 x y) -> hPutStrLn h $ show x ++ "  " ++ show y) xs
  hClose h

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

-- ============================= Testing R1 (Euclidian) ==================================

-- -- | Inverse of error function.
-- invERF :: Double -> Double
-- invERF x
--   | x < 0     = -invERF (-x)
--   | otherwise = let
--     a = 0.147
--     y1 = (2 / (pi * a) + log (1 - x * x) / 2)
--     y2 = sqrt (y1 * y1 - ( 1 / a) * log (1 - x * x))
--     in sqrt (y2 - y1)

-- -- | Compute the normal PDF.
-- normalPDF :: Double -> Double -> Double -> Double
-- normalPDF x mu sigma = let
--   dx = x - mu
--   s2 = 2 * sigma * sigma
--   in exp(-(dx * dx) / s2) / (sqrt (2 * pi) * sigma)

-- -- | Generate a random sample from a univariate normal distribution.
-- -- The random input must range [0, 1]
-- normalSample :: Double -> Double -> Double -> Double
-- normalSample mu sigma rnd = mu + sigma * sqrt 2 * invERF (2 * rnd - 1)
