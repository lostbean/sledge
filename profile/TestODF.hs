module TestODF where

import qualified Data.Vector.Unboxed as U

import Hammer.Math.Algebra
import Hammer.VTK

import Texture.Orientation
import Texture.Symmetry
import Texture.IsoSphere
import Texture.TesseractGrid
import Texture.Kernel
import Texture.Sampler
import Texture.ODF

testODF :: Int -> IO ()
testODF n = do
  -- Initial ODF
  let
    odf0 = buildEmptyODF (Deg 2.5) Cubic (Deg 3)
    qs1  = U.map (toQuaternion . mkAxisPair (Vec3 1 1 1) . Deg) $ U.enumFromStepN (-45) 2 45
    qs2  = U.map (toQuaternion . mkAxisPair (Vec3 (-1) (-1) 1) . Deg) $ U.enumFromStepN (-45) 2 45
    odf1 = addPoints (qs1 U.++ qs2) odf0
    vtkodf1 = renderODFVTK odf1

  writeUniVTKfile ("/home/edgar/InputODF" ++ ".vtu") True vtkodf1

  -- Sampling ODF
  xs <- hitAndRunSlice defaultCfg (getODFeval odf1) zerorot n
  let vtkodf2 = renderODFVTK $ addPoints (U.fromList xs) odf0

  writeUniVTKfile ("/home/edgar/OutputODF" ++ ".vtu") True vtkodf2
