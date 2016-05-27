module Main where

import qualified Data.Vector.Unboxed as U

import Linear.Vect
import Hammer.VTK

import qualified Data.BlazeVPtree as VP

import Texture.Orientation
import Texture.Symmetry
import Texture.IsoSphere
import Texture.TesseractGrid
import Texture.Kernel
import Texture.Sampler
import Texture.ODF ()

main :: IO ()
main = do
  testKernelSampling 1000

-- | Test VP tree on Rodrigues-Frank space
testKernelSampling :: Int -> IO ()
testKernelSampling n = do
  -- Initial ODF
  let
    vs = U.filter (isInRodriFZ Cubic) $ genQuaternionGrid 80
    vp = VP.fromVector vs

  let
    qs1 = U.map (toQuaternion . mkAxisPair (Vec3   1    1  1) . Deg) $ U.enumFromStepN (-45) 2 45
    qs2 = U.map (toQuaternion . mkAxisPair (Vec3 (-1) (-1) 1) . Deg) $ U.enumFromStepN (-45) 2 45
    qs = qs1 U.++ qs2
    ks = addManyKernels (Deg 3) vp qs (U.replicate (U.length vs) 0)

  renderODF (vs, ks) "/home/edgar/InputODF"

  -- Sampling ODF
  let pdf x = maybe 0 (\((i,_,_)) -> ks U.! i) (VP.nearestThanNeighbor vp (4/80) x)
  xs <- hitAndRunSlice defaultCfg pdf zerorot n
  let ts = addManyKernels (Deg 1.5) vp (U.fromList xs) (U.replicate (U.length vs) 0)

  renderODF (vs, ts) "/home/edgar/OutputODF"

-- | Render ODF
renderODF :: (U.Vector Quaternion, U.Vector Double) -> String -> IO ()
renderODF (qs, vs) file = let
  attr = mkPointValueAttr "Intensity" (\i _ -> vs U.! i)
  vtk  = renderQuaternions qs []
  in writeUniVTKfile (file ++ ".vtu") True (addPointValueAttr vtk attr)
