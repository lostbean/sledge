module IPFTest (
    test,
) where

import Linear.Vect
import Test.Tasty
import Test.Tasty.HUnit
import Texture.IPF
import Texture.Orientation
import Texture.Symmetry

test :: TestTree
test =
    testGroup
        "IPF"
        [ testCase "IPF color for cubic e3" testIPFCubicE3
        ]

testIPFCubicE3 :: Assertion
testIPFCubicE3 = do
    let
        symm = Cubic
        ref = ND
        q = zerorot -- Identity rotation
        (_, RGBDoubleColor (r, g, b)) = getIPFColor symm ref q
    -- For Cubic and ND with identity rotation, the inverse pole is e3.
    -- e3 in our IPF schema for Cubic corresponds to [0 0 1], which is the Red direction.
    -- So we expect Red to be dominant.
    assertBool "Red is dominant" (r >= g && r >= b)
