module TestKernel where

import qualified Data.Vector.Unboxed as U
import qualified Data.Vector         as V
import qualified Data.List           as L

import           Data.Function       (on)

import Hammer.Math.Algebra
import Hammer.VTK
import System.Random

import qualified Data.KDtree  as KD
import qualified Data.VPtree  as VP
import qualified Data.MVPtree as MVP

import Texture.Orientation
import Texture.Symmetry
import Texture.IsoSphere
import Texture.TesseractGrid
import Texture.Kernel

import Data.Time.Clock


instance VP.Metric Quaternion where
  dist = getMisoAngle Cubic

instance MVP.Metric Quaternion where
  dist = getMisoAngle Cubic

-- | Generate random values of Rodrigues-Frank C4 symmetry (90 <100>, 90 <010> 90 <001>)
-- with FZ planes at (+-45 <100>,+-45 <010>,+-45 <001>)
getFRFZ :: IO [Quaternion]
getFRFZ = newStdGen >>= \g -> do
  return $ map (toFZ Cubic) (randoms g)

-- | Test VP tree on Rodrigues-Frank space
testKernel :: IO ()
testKernel = do
  --rs  <- fmap (U.fromList . take 10000) getFRFZ
  --let rs = genIsoSphereSO3Grid Cubic (Deg 5)
  let rs = U.filter (isInRodriFZ Cubic) $ genQuaternionGrid 80
  let
    vp = VP.fromVector rs
    t = toQuaternion $ mkAxisPair (Vec3 1 0 0) (Deg 44)
    d = fromAngle (Deg 7.5)

    vpo = VP.nearNeighbors vp d t
    vplist = L.sort $ map (\(i,_,_) -> i) vpo

    bfos = U.filter ((< d) . VP.dist t . snd) (U.imap (,) rs)
    nbfs@(ibfs,_) = U.minimumBy (compare `on` (VP.dist t . snd)) bfos
    bflists = L.sort $ map fst $ U.toList bfos

  putStrLn "=========VP tree==========="
  putStrLn $ "ix: " ++ show vplist
  putStrLn $ "#: " ++ show (length vplist)

  putStrLn "==========Brutal Force=========="
  putStrLn $ "ix: " ++ show bflists
  putStrLn $ "#: " ++ show (length bflists)

  putStrLn "========== Checking =========="
  putStrLn $ "check nears: " ++ show (bflists == vplist)

  {--
  putStrLn "========== Performance =========="
  let nsample =  10000
  ps <- fmap (U.fromList . take nsample) getFRFZ

  ta0 <- getCurrentTime
  putStr $ "Sampling with brutal force " ++ show nsample ++ " rotations. Total points: "
  let func p = U.filter ((< d) . VP.dist p . snd) (U.imap (,) rs)
  putStrLn $ show $ U.sum $ U.map (U.length . func) ps
  ta1 <- getCurrentTime
  putStrLn $ "Calculation time: " ++ show (diffUTCTime ta1 ta0)

  tc0 <- getCurrentTime
  putStr $ "Searching nearest point in VP tree "
  putStrLn $ show $ U.sum $ U.map ((\(Just (i,_,_)) -> i) . VP.nearestNeighbor vp) ps
  tc1 <- getCurrentTime
  putStrLn $ "Calculation time: " ++ show (diffUTCTime tc1 tc0)

  td0 <- getCurrentTime
  putStr $ "Sampling with VP tree " ++ show nsample ++ " rotations. Total points: "
  putStrLn $ show $ U.sum $ U.map (length . VP.nearNeighbors vp d) ps
  td1 <- getCurrentTime
  putStrLn $ "Calculation time: " ++ show (diffUTCTime td1 td0)
  --}

  -- View Gaussian kernel
  --qs <- fmap (U.fromList . take 10) getFRFZ
  let qs = U.map (toQuaternion . mkAxisPair (Vec3 1 1 1) . Deg) $ U.enumFromStepN (-45) 2 45
  let
    vs   = addManyKernels (Deg 3) vp qs (U.replicate (U.length rs) 0)
    attr = mkPointValueAttr "Intensity" (\i _ -> vs U.! i)
    vtk  = renderQuaternions rs []
  writeUniVTKfile "/home/edgar/kernel.vtu" True (addPointValueAttr vtk attr)
