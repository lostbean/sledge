{-# LANGUAGE DeriveGeneric   #-}
{-# LANGUAGE RecordWildCards #-}

module Texture.TesseractGrid
       ( TesseractPoint (..)
       , TesseractGrid (..)
       , tesseractToQuaternion
       , quaternionToTesseract
       , volumeGridCell
       , genQuaternionGrid
       , genTesseractGrid2
       , genTesseractGrid
       , getTesseractPos
       , emptyTesseract
       , binningTesseract
       , accuTesseract
       , updateTesseract
       , tesseractZipWith
       , tesseractMap
       , maxTesseractPoint
       , plotTesseract
       , plotTesseractPoints

       , testNormalization
       , printTesseract
       ) where

import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Data.Vector  (Vector)
import           GHC.Generics (Generic)

import           Data.Binary
import           System.Random

import           Hammer.Math.Algebra
import           Hammer.VTK
import           Texture.Orientation

--import           Debug.Trace
--dbg s x = trace (s ++ show x) x

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

-- ============================= Basic Functions ========================================= 

tesseractToQuaternion :: TesseractPoint -> Quaternion
tesseractToQuaternion tp = case tp of
  Cell1 v -> mkQuaternion $ Vec4 1 (_1 v) (_2 v) (_3 v)
  Cell2 v -> mkQuaternion $ Vec4 (_1 v) 1 (_2 v) (_3 v)
  Cell3 v -> mkQuaternion $ Vec4 (_1 v) (_2 v) 1 (_3 v)
  Cell4 v -> mkQuaternion $ Vec4 (_1 v) (_2 v) (_3 v) 1

quaternionToTesseract :: Quaternion -> TesseractPoint
quaternionToTesseract q
  | i == 0    = Cell1 $ Vec3 (q2/qmax) (q3/qmax) (q4/qmax)
  | i == 1    = Cell2 $ Vec3 (q1/qmax) (q3/qmax) (q4/qmax)
  | i == 2    = Cell3 $ Vec3 (q1/qmax) (q2/qmax) (q4/qmax)
  | otherwise = Cell4 $ Vec3 (q1/qmax) (q2/qmax) (q3/qmax)
  where
    (q1, qv) = splitQuaternion q
    Vec3 q2 q3 q4 = qv
    qs   = U.fromList [q1, q2, q3, q4]
    i    = U.maxIndex $ U.map abs qs
    qmax = qs U.! i

volumeGridCell :: TesseractPoint -> Double
volumeGridCell tp = let
  foo v = let v2 = 1 + v &. v in v2*v2
  in case tp of
    Cell1 v -> foo v
    Cell2 v -> foo v
    Cell3 v -> foo v
    Cell4 v -> foo v

tessPosToIx :: Int -> (Int, Int, Int) -> Int
tessPosToIx m (x, y, z) = m * m * z + m * y + x

ixToTessPos :: Int -> Int -> (Int, Int, Int)
ixToTessPos m ix = (x, y, z)
  where
    --isOutRange = ix < 0 || ix >= (m * m * m)
    (z, a) = quotRem ix (m * m)
    (y, x) = quotRem a m

getTesseractPos :: Int -> TesseractPoint -> (Int, (Int, Int, Int))
getTesseractPos m tp = let
  step = 2 / (fromIntegral m)
  func = min (m-1) . floor . (\x -> (x + 1) / step)
  foo c (Vec3 x y z) = (c, (func x, func y, func z))
  in case tp of
    Cell1 v -> foo 1 v
    Cell2 v -> foo 2 v
    Cell3 v -> foo 3 v
    Cell4 v -> foo 4 v

getTesseractPoint :: Int -> (Int, (Int, Int, Int)) -> TesseractPoint
getTesseractPoint m (cell, (ix, iy, iz))
  | cell <= 1 = Cell1 v
  | cell == 2 = Cell2 v
  | cell == 3 = Cell3 v
  | otherwise = Cell4 v
  where
    step = 2 / (fromIntegral m)
    foo = (\x -> x - 1 + 0.5 * step) . (* step) . fromIntegral
    v = Vec3 (foo ix) (foo iy) (foo iz)

-- ============================== Grid functions ========================================= 

genQuaternionGrid :: Int -> U.Vector Quaternion
genQuaternionGrid m = let
  step = 2 / (fromIntegral m)
  x0   = step/2 - 1
  foo  = (x0 +) . (step *) . fromIntegral
  lin  = [foo i | i <- [0 .. (m - 1)]]
  -- Use X -> Y -> Z sequence (X is fast increment)
  cell = U.fromList [ Vec3 x y z
                    | z <- lin
                    , y <- lin
                    , x <- lin ]
  in (U.map (tesseractToQuaternion . Cell1) cell) U.++
     (U.map (tesseractToQuaternion . Cell2) cell) U.++
     (U.map (tesseractToQuaternion . Cell3) cell) U.++
     (U.map (tesseractToQuaternion . Cell4) cell)

genTesseractGrid2 :: Int -> (Quaternion -> a) -> TesseractGrid a
genTesseractGrid2 m func = genTesseractGrid m (func . tesseractToQuaternion)

genTesseractGrid :: Int -> (TesseractPoint -> a) -> TesseractGrid a
genTesseractGrid m func = let
  step = 2 / (fromIntegral m)
  x0   = step/2 - 1
  foo  = (x0 +) . (step *) . fromIntegral
  lin  = [foo i | i <- [0 .. (m - 1)]]
  -- Use X -> Y -> Z sequence (X is fast increment)
  cell = V.fromList [ Vec3 x y z
                    | z <- lin
                    , y <- lin
                    , x <- lin ]
  in TesseractGrid { gridSize = m
                   , cell1    = V.map (func . Cell1) cell
                   , cell2    = V.map (func . Cell2) cell
                   , cell3    = V.map (func . Cell3) cell
                   , cell4    = V.map (func . Cell4) cell
                   }

binningTesseract :: (Num a)=> Vector Quaternion -> TesseractGrid a -> TesseractGrid a
binningTesseract qs = accuTesseract qs (V.replicate (V.length qs) 1)

accuTesseract :: (Num a)=> Vector Quaternion -> Vector a -> TesseractGrid a -> TesseractGrid a
accuTesseract qs vs tess@TesseractGrid{..} = let
  us = V.map (getTesseractPos gridSize . quaternionToTesseract) $ V.take (V.length vs) qs
  vp = V.map (tessPosToIx gridSize . snd) us
  func cell i = let
    is = V.findIndices ((== i) . fst) us
    vi = V.map (vp V.!) is
    vx = V.map (vs V.!) is
    in V.accumulate_ (+) cell vi vx
  in tess { cell1 = func cell1 1
          , cell2 = func cell2 2
          , cell3 = func cell3 3
          , cell4 = func cell4 4
          }

updateTesseract :: (Num a)=> Vector Quaternion -> Vector a -> TesseractGrid a -> TesseractGrid a
updateTesseract qs vs tess@TesseractGrid{..} = let
  us = V.map (getTesseractPos gridSize . quaternionToTesseract) $ V.take (V.length vs) qs
  vp = V.map (tessPosToIx gridSize . snd) us
  func cell i = let
    is = V.findIndices ((== i) . fst) us
    vi = V.map (vp V.!) is
    vx = V.map (vs V.!) is
    in V.update_ cell vi vx
  in tess { cell1 = func cell1 1
          , cell2 = func cell2 2
          , cell3 = func cell3 3
          , cell4 = func cell4 4
          }

emptyTesseract :: Int -> a -> TesseractGrid a
emptyTesseract m x = let
  s = abs $ m * m * m
  in TesseractGrid
     { gridSize = abs m
     , cell1 = V.replicate s x
     , cell2 = V.replicate s x
     , cell3 = V.replicate s x
     , cell4 = V.replicate s x
     }

-- ============================= Traversing functions ====================================

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

maxTesseractPoint :: (Ord a)=> TesseractGrid a -> TesseractPoint
maxTesseractPoint TesseractGrid{..} = let
  m1 = V.maxIndex cell1
  m2 = V.maxIndex cell2
  m3 = V.maxIndex cell3
  m4 = V.maxIndex cell4
  c  = V.maxIndex $ V.fromList [ cell1 V.! m1, cell2 V.! m2
                               , cell3 V.! m3, cell4 V.! m4]
  ms = V.fromList [m1, m2, m3, m4]
  in (getTesseractPoint gridSize (c+1, ixToTessPos gridSize $ ms V.! c))

-- ============================= Plotting function ======================================= 

plotTesseract :: (RenderElemVTK a)=> TesseractGrid a -> VTK a
plotTesseract t = let
  m     = gridSize t
  step  = 2 / (fromIntegral m)
  orig  = 0.5 * step - 1
  attr1 = mkPointAttr "Cell-1" ((cell1 t) V.!)
  attr2 = mkPointAttr "Cell-2" ((cell2 t) V.!)
  attr3 = mkPointAttr "Cell-3" ((cell3 t) V.!)
  attr4 = mkPointAttr "Cell-4" ((cell4 t) V.!)
  attrs = [attr1, attr2, attr3, attr4]
  in mkSPVTK "Tesseract" (m+1, m+1, m+1) (orig, orig, orig) (step, step, step) attrs

plotTesseractPoints :: U.Vector Quaternion -> VTK Vec3
plotTesseractPoints qs = let
  (cs, ps) = U.unzip $ U.map (foo . quaternionToTesseract) qs
  attr  = mkPointValueAttr  "Cell" (\i _ -> cs U.! i)
  foo x = case x of
    Cell1 v -> (1 :: Int, v)
    Cell2 v -> (2, v)
    Cell3 v -> (3, v)
    Cell4 v -> (4, v)
  in mkUGVTK "Tesseract" ps (V.generate (U.length ps) id) [attr] []

-- ============================= Test Function =========================================== 

testBinning :: IO (TesseractGrid Double)
testBinning = let
  n  = 10000000
  m  = 10
  genqs = V.fromList . take n . randoms
  t0 = emptyTesseract m 0
  in do
    gen <- newStdGen
    return $ binningTesseract (genqs gen) t0

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
printTesseract t fname = let
  vtk = plotTesseract t
  in writeUniVTKfile ("/home/edgar/Desktop/" ++ fname ++ ".vti") True vtk

test :: IO (Bool)
test = do
  let m = 10
  c <- randomRIO (1, 4)
  n <- randomRIO (0, m*m*m - 1)
  let
    t0 = emptyTesseract m 0
    p = getTesseractPoint m (c, ixToTessPos m n)
    q = tesseractToQuaternion p
    t1 :: TesseractGrid Double
    t1 = binningTesseract (V.singleton q) t0
    cs = V.fromList [cell1 t1, cell2 t1, cell3 t1, cell4 t1]
    cm = V.maxIndex (V.map V.maximum cs)
    nm = V.maxIndex (cs V.! cm)
    qm = tesseractToQuaternion $ maxTesseractPoint t1
    vtk = plotTesseractPoints (U.fromList [q, qm])
  printTesseract t1 "tessTest"
  writeUniVTKfile ("/home/edgar/Desktop/tessTest.vtu") True vtk
  print (q, qm)
  putStr "(cell, linear index, pos) = " >> print (c, n, ixToTessPos m n)
  putStr "TessPoint = "            >> print p
  putStr "TessPoint(recalc.) = "   >> print (quaternionToTesseract q)
  putStr "linear index(recalc) = " >> print (getTesseractPos m p)
  putStr "(cell, linear index)(recalc) = " >> print (cm+1, nm)
  return $ c == (cm+1) && n == nm && q == qm
