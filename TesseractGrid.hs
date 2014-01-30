{-# LANGUAGE DeriveGeneric #-}

module Hammer.Texture.TesseractGrid
       ( TesseractPoint (..)
       , TesseractGrid (..)
       , tesseractToQuaternion
       , quaternionToTesseract
       , volumeGridCell
       , genTesseractGrid2
       , genTesseractGrid
       , getTesseractIx
       , binningTesseract
       , tesseractZipWith
       , tesseractMap

       , testNormalization
       , printTesseract
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Data.Vector  (Vector)
import           GHC.Generics (Generic)

import           Data.Vector.Binary
import           Data.Binary
import           System.Random

import           Hammer.Math.Algebra
import           Hammer.Texture.Orientation
import           Hammer.Render.VTK.VTKRender

import           Debug.Trace
dbg s x = trace (s L.++ show x) x

data TesseractPoint
  = Cell1 Vec3
  | Cell2 Vec3
  | Cell3 Vec3
  | Cell4 Vec3
  deriving (Show)

data TesseractGrid a
  = TesseractGrid
    { gridSize :: Int
    , cell1    :: Vector a
    , cell2    :: Vector a
    , cell3    :: Vector a
    , cell4    :: Vector a
    } deriving (Show, Generic)

instance (Binary a)=> Binary (TesseractGrid a)

tesseractToQuaternion :: TesseractPoint -> Quaternion
tesseractToQuaternion tp = case tp of
  Cell1 v -> mkQuaternion $ Vec4 1 (_1 v) (_2 v) (_3 v)
  Cell2 v -> mkQuaternion $ Vec4 (_1 v) 1 (_2 v) (_3 v)
  Cell3 v -> mkQuaternion $ Vec4 (_1 v) (_2 v) 1 (_3 v)
  Cell4 v -> mkQuaternion $ Vec4 (_1 v) (_2 v) (_3 v) 1

quaternionToTesseract :: Quaternion -> TesseractPoint
quaternionToTesseract q
  | qa1 == qmax = Cell1 $ Vec3 (q2/q1) (q3/q1) (q4/q1)
  | qa2 == qmax = Cell2 $ Vec3 (q1/q2) (q3/q2) (q4/q2)
  | qa3 == qmax = Cell3 $ Vec3 (q1/q3) (q2/q3) (q4/q3)
  | otherwise   = Cell4 $ Vec3 (q1/q4) (q2/q4) (q3/q4)
  where
    (q1, qv) = splitQuaternion q
    q2  = _1 qv
    q3  = _2 qv
    q4  = _3 qv
    qa1 = abs q1
    qa2 = abs q2
    qa3 = abs q3
    qa4 = abs q4
    qmax = max (max qa1 qa2) (max qa3 qa4)

volumeGridCell :: TesseractPoint -> Double
volumeGridCell tp = let
  foo v = let v2 = 1 + v &. v in v2*v2
  in case tp of
    Cell1 v -> foo v
    Cell2 v -> foo v
    Cell3 v -> foo v
    Cell4 v -> foo v

genTesseractGrid2 :: Int -> (Quaternion -> a) -> TesseractGrid a
genTesseractGrid2 m func = genTesseractGrid m (func . tesseractToQuaternion)

genTesseractGrid :: Int -> (TesseractPoint -> a) -> TesseractGrid a
genTesseractGrid m func = let
  step = 2 / (fromIntegral m)
  x0   = step/2 - 1
  foo  = (x0 +) . (step *) . fromIntegral
  linX = [foo i | i <- [0 .. (m - 1)]]
  cell = V.fromList [ Vec3 x y z
                    | x <- linX
                    , y <- linX
                    , z <- linX ]
  in TesseractGrid { gridSize = m
                   , cell1    = V.map (func . Cell1) cell
                   , cell2    = V.map (func . Cell2) cell
                   , cell3    = V.map (func . Cell3) cell
                   , cell4    = V.map (func . Cell4) cell
                   }

getTesseractIx :: Int -> TesseractPoint -> (Int, Int, Int, Int)
getTesseractIx m tp = let
  step = 2 / (fromIntegral m)
  func = min (fromIntegral $ m - 1) . floor . (/ step) . (+ 1)
  foo c (Vec3 x y z)
    | otherwise       = (c, func x, func y, func z)
  in case tp of
    Cell1 v -> foo 0 v
    Cell2 v -> foo 1 v
    Cell3 v -> foo 2 v
    Cell4 v -> foo 3 v

binningTesseract :: Int -> Vector Quaternion -> TesseractGrid Double
binningTesseract m qs = let
  pos2ix (c, x, y, z) = m*m*m*c + m*m*x + m*y + z
  us = V.map (pos2ix . getTesseractIx m . quaternionToTesseract) qs
  v0 = V.replicate (4*m*m*m) 0
  v1 = V.replicate (V.length qs) 1
  ts = V.accumulate_ (+) v0 us v1
  getRange i = V.slice (pos2ix (i, 0, 0, 0)) (m*m*m) ts
  in TesseractGrid { gridSize = m
                   , cell1    = getRange 0
                   , cell2    = getRange 1
                   , cell3    = getRange 2
                   , cell4    = getRange 3
                   }

tesseractZipWith :: (a -> b -> c) -> TesseractGrid a -> TesseractGrid b -> TesseractGrid c
tesseractZipWith func ta tb
  | gridSize ta /= gridSize tb = error "[TesseractGrid] Can't zip different grid sizes."
  | otherwise = TesseractGrid { gridSize = gridSize ta
                              , cell1    = V.zipWith func (cell1 ta) (cell1 tb)
                              , cell2    = V.zipWith func (cell2 ta) (cell2 tb)
                              , cell3    = V.zipWith func (cell3 ta) (cell3 tb)
                              , cell4    = V.zipWith func (cell4 ta) (cell4 tb)
                              }

tesseractMap :: (a -> b) -> TesseractGrid a -> TesseractGrid b
tesseractMap func ta = ta { cell1 = V.map func (cell1 ta)
                          , cell2 = V.map func (cell2 ta)
                          , cell3 = V.map func (cell3 ta)
                          , cell4 = V.map func (cell4 ta)
                          }

-- ============================= Test Function =========================================== 

testBinning :: IO (TesseractGrid Double)
testBinning = let
  n  = 10000000
  m  = 10
  genqs = V.fromList . take n . randoms
  in do
    gen <- newStdGen
    return $ binningTesseract m (genqs gen)

testNormalization :: IO (TesseractGrid Double)
testNormalization = let
  n  = 10000000 :: Int
  m  = 10 :: Int
  ns = genTesseractGrid m volumeGridCell
  k  = (fromIntegral n) / (fromIntegral $ 4*m*m*m)
  in do
    ts <- testBinning
    return $ tesseractMap (/k) $ tesseractZipWith (*) ns ts

printTesseract :: (RenderElemVTK a)=> TesseractGrid a -> String -> IO ()
printTesseract t name = let
  m     = gridSize t
  step  = 2 / (fromIntegral m)
  vtk   = mkSPVTK "Tesseract" (m, m, m) (-1, -1, -1) (step, step, step)
  attr1 = mkPointAttr "Cell-1" (\i _ -> (cell1 t) V.! i)
  attr2 = mkPointAttr "Cell-2" (\i _ -> (cell2 t) V.! i)
  attr3 = mkPointAttr "Cell-3" (\i _ -> (cell3 t) V.! i)
  attr4 = mkPointAttr "Cell-4" (\i _ -> (cell4 t) V.! i)
  vtk2  = L.foldl' (\acc attr -> addDataPoints acc attr) vtk [attr1, attr2, attr3, attr4]
  in writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vti") True vtk2
