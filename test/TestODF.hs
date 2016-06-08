module Main where

import qualified Data.Vector.Unboxed as U

import Linear.Vect
import Hammer.VTK

import Texture.Orientation
import Texture.Symmetry
import Texture.IsoSphere
import Texture.Sampler
import Texture.ODF

main :: IO ()
main = do
  testODF 1000

testODF :: Int -> IO ()
testODF n = do
  -- Initial ODF
  let
    odf0 = buildEmptyODF (Deg 2.5) Cubic (Deg 3)
    qs1  = U.map (toQuaternion . mkAxisPair (Vec3   1    1   1 ) . Deg) $ U.enumFromStepN (-45) 2    45
    qs2  = U.map (toQuaternion . mkAxisPair (Vec3 (-1) (-1)  1 ) . Deg) $ U.enumFromStepN (-45) 2    45
    qs3  = U.map (toQuaternion . mkAxisPair (Vec3   1  (-1)  1 ) . Deg) $ U.enumFromStepN 30    0.75 33
    qs4  = U.map (toQuaternion . mkAxisPair (Vec3 (-1)   1 (-1)) . Deg) $ U.enumFromStepN 30    0.75 33
    odf1 = addPoints (qs1 U.++ qs2 U.++ qs3 U.++ qs4) odf0
    vtkodf1 = renderODFVTK odf1

  writeUniVTKfile ("InputODF" ++ ".vtu") True vtkodf1

  -- Sampling ODF
  xs <- hitAndRunSlice defaultCfg (getODFeval odf1) mempty n
  let
    vtkodf2 = renderODFVTK $ addPoints (U.fromList xs) odf0
    vtkxs   = renderQuaternions (U.map (toFZ Cubic) $ U.fromList xs) []

  writeUniVTKfile ("SampledPoints" ++ ".vtu") True vtkxs
  writeUniVTKfile ("OutputODF"     ++ ".vtu") True vtkodf2
