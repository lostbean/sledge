module BinghamTest (
    test,
) where

import Linear.Vect
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck
import TestOrphans ()
import Texture.Bingham
import Texture.Orientation

test :: TestTree
test =
    testGroup
        "Bingham"
        [ testCase "Bingham PDF normalization (rough check)" testNormalization
        , testProperty "Bingham PDF is positive" prop_pdfPositive
        , testProperty "Bingham mode has highest PDF" prop_modeIsMax
        ]

{- | Rough check of normalization using a few points.
A better check would be Monte Carlo integration, but that's slow for unit tests.
-}
testNormalization :: Assertion
testNormalization = do
    let
        d1 = (10, mkQuaternion (Vec4 1 0 0 0))
        d2 = (5, mkQuaternion (Vec4 0 1 0 0))
        d3 = (2, mkQuaternion (Vec4 0 0 1 0))
        dist = mkBingham d1 d2 d3
    -- Just check it doesn't crash and returns a reasonable value
    assertBool "Normalization is positive" (normalization dist > 0)

prop_pdfPositive :: Double -> Double -> Double -> Quaternion -> Property
prop_pdfPositive z1 z2 z3 q =
    let
        -- mkBingham expects (concentration, direction)
        -- concentrations are relative to the mode which will have 0.
        dist =
            mkBingham
                (z1, mkQuaternion (Vec4 1 0 0 0))
                (z2, mkQuaternion (Vec4 0 1 0 0))
                (z3, mkQuaternion (Vec4 0 0 1 0))
        pdf = binghamPDF dist q
     in
        counterexample ("PDF: " ++ show pdf) $ pdf >= 0

prop_modeIsMax :: Double -> Double -> Double -> Quaternion -> Property
prop_modeIsMax z1 z2 z3 q =
    let
        dist =
            mkBingham
                (z1, mkQuaternion (Vec4 1 0 0 0))
                (z2, mkQuaternion (Vec4 0 1 0 0))
                (z3, mkQuaternion (Vec4 0 0 1 0))
        pdfMode = binghamPDF dist (bingMode dist)
        pdfQ = binghamPDF dist q
     in
        counterexample ("PDF Mode: " ++ show pdfMode ++ ", PDF Q: " ++ show pdfQ) $
            pdfMode >= pdfQ - 1e-9
