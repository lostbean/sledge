module TestKernel where

import qualified Data.Vector.Unboxed as U
import qualified Data.Vector         as V
import qualified Data.List           as L

import           Data.Function       (on)

import Hammer.Math.Algebra
import Hammer.VTK
import System.Random

import qualified Data.BlazeVPtree as VP

import Texture.Orientation
import Texture.Symmetry
import Texture.IsoSphere
import Texture.TesseractGrid
import Texture.Kernel

import Data.Time.Clock

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

  putStrLn "========== Performance =========="
  let nsample =  10000
  ps <- fmap (U.fromList . take nsample) getFRFZ

  printTime "brutal force" $ let
    func p = U.filter ((< d) . VP.dist p . snd) (U.imap (,) rs)
    in U.map (U.length . func) ps

  printTime "VP tree" $ U.map (length . VP.nearNeighbors vp d) ps

  printTime "VP tree" $ U.map ((\(Just (i,_,_)) -> i) . VP.nearestThanNeighbor vp (4/80)) ps

  -- View Gaussian kernel
  --qs <- fmap (U.fromList . take 10) getFRFZ
  let qs = U.map (toQuaternion . mkAxisPair (Vec3 1 1 1) . Deg) $ U.enumFromStepN (-45) 2 45
  let
    vs   = addManyKernels (Deg 3) vp qs (U.replicate (U.length rs) 0)
    attr = mkPointValueAttr "Intensity" (\i _ -> vs U.! i)
    vtk  = renderQuaternions rs []
  writeUniVTKfile "/home/edgar/kernel.vtu" True (addPointValueAttr vtk attr)

printTime :: (U.Unbox s, Show s, Num s)=> [Char] -> U.Vector s -> IO ()
printTime name ps = do
  ta0 <- getCurrentTime
  putStr $ "Sampling with " ++ name ++ ". Total points: "
  putStrLn $ show $ U.sum ps
  ta1 <- getCurrentTime
  putStrLn $ "Calculation time: " ++ show (diffUTCTime ta1 ta0)
