module BinghamTest (
    test,
) where

import qualified Data.Vector.Unboxed as U
import Linear.Vect
import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck
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

-- We need Arbitrary instance for Quaternion, but it's already in TestTexture.
-- To avoid duplication, we might want to move these instances to a shared module.
-- For now, I'll just use the one from TestTexture if I can, or re-define here if needed.
-- Since they are orphans in TestTexture, I can't easily import them unless I import TestTexture.

instance Arbitrary Quaternion where
    arbitrary = do
        w <- choose (-1, 1)
        x <- choose (-1, 1)
        y <- choose (-1, 1)
        z <- choose (-1, 1)
        return $ mkQuaternion (Vec4 w x y z)
