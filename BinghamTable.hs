{-# LANGUAGE BangPatterns      #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}

module Hammer.Texture.BinghamTable
       ( tableIndex

       , interpolateF3D
       , interpolatedFdZ13D
       , interpolatedFdZ23D
       , interpolatedFdZ33D

       , interpolateF2D
       , interpolatedFdZ12D
       , interpolatedFdZ22D

       , interpolateF1D
       , interpolatedFdZ1D

       , findNearestNaive

       , writeAllTables
       ) where

import qualified Data.Vector.Unboxed as VU

import           Data.Vector.Unboxed (Vector)
import           Data.Bits           (shiftR)

import           Data.Binary
import           Data.Array.Unboxed
import           System.IO.Unsafe

import           Hammer.Math.Algebra
import           Hammer.Texture.BinghamNormalization

import           Paths_hammer
--getDataFileName = undefined

-- ============================= Load table from file ====================================

readTable :: (Binary a)=> String -> a
readTable tableName = unsafePerformIO . unsafeInterleaveIO $ do
  file <- getDataFileName $ "datatables/" ++ tableName ++ ".data"
  decodeFile file

-- ================================= Lookup tables =======================================

tableF3D :: UArray (Int, Int, Int) Double
tableF3D = readTable "tableF3D"

tabledFdZ13D :: UArray (Int, Int, Int) Double
tabledFdZ13D = readTable "tabledFdZ1_3D"

tabledFdZ23D :: UArray (Int, Int, Int) Double
tabledFdZ23D =  readTable "tabledFdZ2_3D"

tabledFdZ33D :: UArray (Int, Int, Int) Double
tabledFdZ33D =  readTable "tabledFdZ3_3D"


tableF2D :: UArray (Int, Int) Double
tableF2D = readTable "tableF2D"

tabledFdZ12D :: UArray (Int, Int) Double
tabledFdZ12D = readTable "tabledFdZ1_2D"

tabledFdZ22D :: UArray (Int, Int) Double
tabledFdZ22D =  readTable "tabledFdZ2_2D"


tableF1D :: UArray Int Double
tableF1D = readTable "tableF1D"

tabledFdZ1D :: UArray Int Double
tabledFdZ1D = readTable "tabledFdZ_1D"

findNearestNaive :: Vec3 -> ((Int, Int, Int), Vec3)
findNearestNaive x = let
  arr = uvarray tabledY
  p   = VU.minIndex $ VU.map (normsqr . (x &-)) arr
  in (ix (uvbounds tabledY) p, arr VU.! p)

tabledY :: MultiVector (Int, Int, Int) Vec3
tabledY = let
  bd = bounds tableF3D
  func addr
    | f == 0    = zero
    | otherwise = Vec3 (x/f) (y/f) (z/f)
    where
      f = tableF3D ! addr
      x = tabledFdZ13D ! addr
      y = tabledFdZ23D ! addr
      z = tabledFdZ33D ! addr
  vec = VU.generate (rangeSize bd) (func . ix bd)
  in mkMultiVector bd vec

-- ======================================================================================= 

data (IxVector i)=> MultiVector i e = MultiVector
                   { uvbounds :: (i, i)
                   , uvarray  :: Vector e
                   } deriving (Show)

class (Ord a)=> IxVector a where
  ix        :: (a, a) -> Int -> a
  pos       :: (a, a) -> a -> Int
  rangeIx   :: (a, a) -> [a]
  rangePos  :: (a, a) -> [Int]
  ixSize    :: (a, a) -> Int

instance IxVector (Int, Int, Int) where
  ix ((xi, yi, zi), (xf, yf, _)) = let
    dx = 1 + xf - xi
    dy = 1 + yf - yi
    dxy = dx * dy
    in \p -> let
      (z, xy) = quotRem p  dxy
      (y, x)  = quotRem xy dy
      in (xi + x, yi + y, zi + z)
  pos ((xi, yi, zi), (xf, yf, _)) = let
    dx = 1 + xf - xi
    dy = 1 + yf - yi
    dxy = dx * dy
    in \(x,y,z) -> (x - xi) + dx*(y - yi) + dxy*(z - zi)
  rangeIx ((xi, yi, zi), (xf, yf, zf)) = [(i,j,k) | i <- [xi..xf]
                                                  , j <- [yi..yf]
                                                  , k <- [zi..zf]]
  rangePos bd = [0 .. (ixSize bd) - 1]
  ixSize ((xi, yi, zi), (xf, yf, zf)) = let
    dx = 1 + xf - xi
    dy = 1 + yf - yi
    dz = 1 + zf - zi
    in dx * dy * dz

mkMultiVector :: (IxVector i, VU.Unbox e)=> (i, i) -> Vector e -> MultiVector i e
mkMultiVector bd vec
  | sb == sv  = MultiVector bd vec
  | otherwise = error "[MultiVector] Array size and input vector don't match!"
  where
    sb = ixSize bd
    sv = VU.length vec

-- ======================================================================================= 

tableIndex :: Vector Double
tableIndex = VU.fromList [ 0.000000e+00, 1.000000e-01, 2.000000e-01, 3.000000e-01
                         , 4.000000e-01, 5.000000e-01, 6.000000e-01, 7.000000e-01
                         , 8.000000e-01, 9.000000e-01, 1.000000e+00, 1.100000e+00
                         , 1.200000e+00, 1.300000e+00, 1.400000e+00, 1.500000e+00
                         , 1.600000e+00, 1.700000e+00, 1.800000e+00, 1.900000e+00
                         , 2.000000e+00, 2.200000e+00, 2.400000e+00, 2.600000e+00
                         , 2.800000e+00, 3.000000e+00, 3.200000e+00, 3.400000e+00
                         , 3.600000e+00, 3.800000e+00, 4.000000e+00, 4.500000e+00
                         , 5.000000e+00, 5.500000e+00, 6.000000e+00, 6.500000e+00
                         , 7.000000e+00, 7.500000e+00, 8.000000e+00, 8.500000e+00
                         , 9.000000e+00, 9.500000e+00, 1.000000e+01, 1.050000e+01
                         , 1.100000e+01, 1.150000e+01, 1.200000e+01, 1.250000e+01
                         , 1.300000e+01, 1.350000e+01, 1.400000e+01, 1.450000e+01
                         , 1.500000e+01, 1.550000e+01, 1.600000e+01, 1.650000e+01
                         , 1.700000e+01, 1.750000e+01, 1.800000e+01, 1.850000e+01
                         , 1.900000e+01, 1.950000e+01, 2.000000e+01, 2.100000e+01
                         , 2.200000e+01, 2.300000e+01, 2.400000e+01, 2.500000e+01
                         , 2.600000e+01, 2.700000e+01, 2.800000e+01, 2.900000e+01
                         , 3.000000e+01 ]

-- =========================== Interpolate value from table ==============================

binarySearch :: (VU.Unbox a, Ord a) => Vector a -> a -> Int
binarySearch vec x = loop 0 (VU.length vec - 1)
  where
    loop !l !u
      | u <= l    = l
      | otherwise = case compare x' x of
        LT -> loop (k+1) u
        GT -> loop l     k
        EQ -> k
      where
        k  = (u + l) `shiftR` 1
        x' = VU.unsafeIndex vec k

getBounds :: Double -> (Int, Int)
getBounds y
  | i >= maxi = (maxi-1, maxi)
  | i <= 0    = (0     , 1   )
  | otherwise = (i-1   , i   )
  where
    i    = binarySearch tableIndex y
    maxi = VU.length tableIndex - 1

getX :: Int -> Double
getX = (tableIndex VU.!)

interpolateTable3D :: ((Int, Int, Int) -> Double) -> Double -> Double -> Double -> Double
interpolateTable3D table z1 z2 z3 = let
  (il, iu) = getBounds z1
  (jl, ju) = getBounds z2
  (kl, ku) = getBounds z3

  di = getX iu - getX il
  dj = getX ju - getX jl
  dk = getX ku - getX kl

  dz1 = z1 - getX il
  dz2 = z2 - getX jl
  dz3 = z3 - getX kl

  ylll = table (il, jl, kl)
  yllu = table (il, jl, ku)
  ylul = table (il, ju, kl)
  yluu = table (il, ju, ku)
  yull = table (iu, jl, kl)
  yulu = table (iu, jl, ku)
  yuul = table (iu, ju, kl)
  yuuu = table (iu, ju, ku)
  -- interpolate over k
  yll = ylll + dz3 * (yllu - ylll) / dk
  ylu = ylul + dz3 * (yluu - ylul) / dk
  yul = yull + dz3 * (yulu - yull) / dk
  yuu = yuul + dz3 * (yuuu - yuul) / dk
  -- interpolate over j
  yl = yll + dz2 * (ylu - yll) / dj
  yu = yul + dz2 * (yuu - yul) / dj
  -- interpolate over i
  in yl + dz1 * (yu - yl) / di

interpolateTable2D :: ((Int, Int) -> Double) -> Double -> Double -> Double
interpolateTable2D table z1 z2 = let
  (il, iu) = getBounds z1
  (jl, ju) = getBounds z2

  di = getX iu - getX il
  dj = getX ju - getX jl

  dz1 = z1 - getX il
  dz2 = z2 - getX jl

  yll = table (il, jl)
  ylu = table (il, ju)
  yul = table (iu, jl)
  yuu = table (iu, ju)
  -- interpolate over j
  yl = yll + dz2 * (ylu - yll) / dj
  yu = yul + dz2 * (yuu - yul) / dj
  -- interpolate over i
  in yl + dz1 * (yu - yl) / di

interpolateTable1D :: (Int -> Double) -> Double -> Double
interpolateTable1D table z1 = let
  (il, iu) = getBounds z1
  di  = getX iu - getX il
  dz1 = z1 - getX il
  yl = table il
  yu = table iu
  -- interpolate over i
  in yl + dz1 * (yu - yl) / di

{-# INLINE fast3Sort #-}
fast3Sort :: (Ord a)=> (a,a,a) -> (a,a,a)
fast3Sort v@(a, b, c)
  | (a >= b) && (b >= c) = v
  | otherwise            = (a', b', c')
  where
    minab = min a b
    maxab = max a b
    a'    = max maxab c
    b'    = max (min maxab c) minab
    c'    = min minab c

{-# INLINE fast2Sort #-}
fast2Sort :: (Ord a)=> (a,a) -> (a,a)
fast2Sort v@(a, b)
  | a >= b    = v
  | otherwise = (b, a)


interpolateF3D :: Double -> Double -> Double -> Double
interpolateF3D z1 z2 z3 = interpolateTable3D ((tableF3D!) . fast3Sort) z1 z2 z3

interpolatedFdZ13D :: Double -> Double -> Double -> Double
interpolatedFdZ13D z1 z2 z3 = interpolateTable3D ((tabledFdZ13D!) . fast3Sort) z1 z2 z3

interpolatedFdZ23D :: Double -> Double -> Double -> Double
interpolatedFdZ23D z1 z2 z3 = interpolateTable3D ((tabledFdZ23D!) . fast3Sort) z1 z2 z3

interpolatedFdZ33D :: Double -> Double -> Double -> Double
interpolatedFdZ33D z1 z2 z3 = interpolateTable3D ((tabledFdZ33D!) . fast3Sort) z1 z2 z3


interpolateF2D :: Double -> Double -> Double
interpolateF2D z1 z2 = interpolateTable2D ((tableF2D!) . fast2Sort) z1 z2

interpolatedFdZ12D :: Double -> Double -> Double
interpolatedFdZ12D z1 z2 = interpolateTable2D ((tabledFdZ12D!) . fast2Sort) z1 z2

interpolatedFdZ22D :: Double -> Double -> Double
interpolatedFdZ22D z1 z2 = interpolateTable2D ((tabledFdZ22D!) . fast2Sort) z1 z2


interpolateF1D :: Double -> Double
interpolateF1D z1 = interpolateTable1D (tableF1D!) z1

interpolatedFdZ1D :: Double -> Double
interpolatedFdZ1D z1 = interpolateTable1D (tabledFdZ1D!) z1

-- =========================== Table generation ==============================

genTable3D :: (Double -> Double -> Double -> Double) -> UArray (Int, Int, Int) Double
genTable3D foo = let
  n = VU.length tableIndex - 1
  f = (tableIndex VU.!)
  in array ((0, 0, 0), (n, n, n)) [ ((i, j, k), foo (f i) (f j) (f k))
                                  | i <- [0..n]
                                  , j <- [0..n]
                                  , k <- [0..n] ]

genTable2D :: (Double -> Double -> Double) -> UArray (Int, Int) Double
genTable2D foo = let
  n = VU.length tableIndex - 1
  f = (tableIndex VU.!)
  in array ((0, 0), (n, n)) [ ((i, j), foo (f i) (f j))
                            | i <- [0..n]
                            , j <- [0..n] ]

genTable1D :: (Double -> Double) -> UArray Int Double
genTable1D foo = let
  n = VU.length tableIndex - 1
  f = (tableIndex VU.!)
  in array (0, n) [ (i, foo (f i)) | i <- [0..n] ]


writeAllTables  :: IO ()
writeAllTables = do
  writeTableF3D
  writeTabledF13D
  writeTabledF23D
  writeTabledF33D

  writeTableF2D
  writeTabledF12D
  writeTabledF22D

  writeTableF1D
  writeTabledF1D

writeTableF3D :: IO ()
writeTableF3D = encodeFile "tableF3D.data" (genTable3D (computeF 3))

writeTabledF13D :: IO ()
writeTabledF13D = encodeFile "tabledFdZ1_3D.data" (genTable3D (computedFdz1 3))

writeTabledF23D :: IO ()
writeTabledF23D = encodeFile "tabledFdZ2_3D.data" (genTable3D (computedFdz2 3))

writeTabledF33D :: IO ()
writeTabledF33D = encodeFile "tabledFdZ3_3D.data" (genTable3D (computedFdz3 3))


writeTableF2D :: IO ()
writeTableF2D = encodeFile "tableF2D.data" (genTable2D (computeF 3 0.0))

writeTabledF12D :: IO ()
writeTabledF12D = encodeFile "tabledFdZ1_2D.data" (genTable2D (computedFdz1 3 0.0))

writeTabledF22D :: IO ()
writeTabledF22D = encodeFile "tabledFdZ2_2D.data" (genTable2D (computedFdz2 3 0.0))


writeTableF1D :: IO ()
writeTableF1D = encodeFile "tableF1D.data" (genTable1D (computeF 3 0.0 0.0))

writeTabledF1D :: IO ()
writeTabledF1D = encodeFile "tabledFdZ_1D.data" (genTable1D (computedFdz1 3 0.0 0.0))
