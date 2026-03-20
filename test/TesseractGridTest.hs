module TesseractGridTest (
    test,
) where

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck
import TestOrphans ()
import Texture.Orientation
import Texture.TesseractGrid

test :: TestTree
test =
    testGroup
        "TesseractGrid"
        [ testProperty "quaternion-tesseract roundtrip" prop_roundtrip
        , testCase "Tesseract grid point indexing" testIndexing
        ]

prop_roundtrip :: Quaternion -> Property
prop_roundtrip q =
    let
        p = quaternionToTesseract q
        q' = tesseractToQuaternion p
        -- We expect q and q' to represent the same orientation.
        -- Due to discretization in TesseractGrid, there might be some error?
        -- Wait, quaternionToTesseract is a mapping to the 4-cell projection.
        -- Let's check if it's really a roundtrip.
        omega = getOmega (q -@- q')
     in
        counterexample ("Original Q: " ++ show q ++ "\nTessPoint: " ++ show p ++ "\nRecalculated Q: " ++ show q' ++ "\nOmega: " ++ show omega) $
            omega < 1e-9

testIndexing :: Assertion
testIndexing = do
    let m = 10
        c = 2
        pos = (5, 5, 5)
        p = getTesseractPoint m (c, pos)
        (c', pos') = getTesseractPos m p
    assertEqual "Cell index" c c'
    assertEqual "Position" pos pos'
