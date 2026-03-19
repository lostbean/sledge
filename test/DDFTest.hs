module DDFTest (
    test,
) where

import qualified Data.Vector.Unboxed as U
import Linear.Vect
import Test.Tasty
import Test.Tasty.HUnit
import Texture.DDF
import Texture.Orientation
import Texture.Symmetry

test :: TestTree
test =
    testGroup
        "DDF"
        [ testCase "Build and add points to DDF" testDDF
        ]

testDDF :: Assertion
testDDF = do
    let
        kw = Deg 5
        symm = Cubic
        step = Deg 10
        ddf = buildEmptyDDF kw symm step
        -- Create some points near e1
        pts = U.fromList [mkNormal (Vec3 1 0 0), mkNormal (Vec3 1 0.1 0)]
        ddf' = addPoints pts ddf
        eval = getDDFeval ddf'
        val = eval (mkNormal (Vec3 1 0 0))
    assertBool "Intensity at e1 is positive" (val > 0)
    assertEqual "Grid size is positive" True (ddfGridSize ddf' > 0)
