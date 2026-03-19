module ODFTest (
    test,
) where

import qualified Data.Vector.Unboxed as U
import Linear.Vect
import Test.Tasty
import Test.Tasty.HUnit
import Texture.ODF
import Texture.Orientation
import Texture.Symmetry

test :: TestTree
test =
    testGroup
        "ODF"
        [ testCase "Build and add points to ODF" testODF
        ]

testODF :: Assertion
testODF = do
    let
        kw = Deg 5
        symm = Cubic
        step = Deg 10
        odf = buildEmptyODF kw symm step
        -- Create some points near identity
        pts = U.fromList [mkQuaternion (Vec4 1 0 0 0), mkQuaternion (Vec4 0.99 0.01 0 0)]
        odf' = addPoints pts odf
        eval = getODFeval odf'
        val = eval (mkQuaternion (Vec4 1 0 0 0))
    assertBool "Intensity at identity is positive" (val > 0)
    assertEqual "Grid size is positive" True (odfGridSize odf' > 0)
