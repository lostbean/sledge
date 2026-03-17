module BinghamTest (
    testDist,
    testSample,
    testNormalization,
) where

import qualified Data.Vector.Unboxed as U
import Hammer.VTK
import Linear.Vect
import System.Random (newStdGen, randoms)

import Texture.Bingham
import Texture.Bingham.Constant (surface_area_sphere)
import Texture.HyperSphere
import Texture.Orientation

writeQuater :: (RenderElemVTK a) => String -> VTK a -> IO ()
writeQuater name = writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vtu") True

renderPoints :: [Quaternion] -> VTK Vec3D
renderPoints = renderSO3PointsVTK . U.map quaternionToSO3 . U.fromList

testDist :: Bingham
testDist =
    let
        d1 = (30, mkQuaternion (Vec4 0 0 1 0))
        d2 = (2, mkQuaternion (Vec4 0 1 0 0))
        d3 = (1, mkQuaternion (Vec4 1 0 0 0))
        dist = mkBingham d1 d2 d3
     in
        dist

testSample :: Int -> IO ()
testSample n =
    let
        d1 = (5, mkQuaternion (Vec4 0 0 1 0))
        d2 = (2, mkQuaternion (Vec4 0 0 0 (-1)))
        d3 = (1, mkQuaternion (Vec4 0 1 0 0))
        dist = mkBingham d1 d2 d3
        prod = bingProduct dist dist
     in
        do
            a <- sampleBingham dist n
            putStrLn $ show dist
            putStrLn $ showPretty $ scatter dist
            putStrLn $ showPretty $ getScatterMatrix $ U.fromList a
            writeQuater "Bing-PDF-testSample" $ renderBingham dist
            writeQuater "Bing-PDFProduct-testSample" $ renderBingham prod
            writeQuater "Bing-Samples-testSample" $ renderPoints a
            writeQuater "Euler-PDF-testSample" $ renderBinghamToEuler (Deg 5) (72, 36, 72) dist

-- | Use Monte Carlo integration to check to normalization of the distribution.
testNormalization :: Double -> Double -> Double -> IO Double
testNormalization z1 z2 z3 =
    let
        d1 = (z1, mkQuaternion (Vec4 0 0 1 0))
        d2 = (z2, mkQuaternion (Vec4 0 1 0 0))
        d3 = (z3, mkQuaternion (Vec4 0 0 0 1))
        dist = mkBingham d1 d2 d3
     in
        do
            gen <- newStdGen
            let
                n = 1000000
                qs = take n $ randoms gen
                s = sum $ map (binghamPDF dist) qs
                v = surface_area_sphere 3
            return $ (v / fromIntegral n) * s
