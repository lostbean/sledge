module SphereProjTest (
    test,
) where

import Linear.Vect
import Test.Tasty
import Test.Tasty.HUnit
import Texture.SphereProjection

test :: TestTree
test =
    testGroup
        "SphereProjection"
        [ testCase "Lambert projection" testLambert
        , testCase "Stereographic projection" testStero
        ]

testLambert :: Assertion
testLambert = do
    let n = mkNormal (Vec3 0 0 1)
        p = lambertSO3Proj n
    assertBool "Upper projection" (isUpperSO3 p)
    assertEqual "Center" (Vec2 0 0) (so3ProjCoord p)

testStero :: Assertion
testStero = do
    let n = mkNormal (Vec3 0 0 1)
        p = steroSO3Proj n
    assertBool "Upper projection" (isUpperSO3 p)
    assertEqual "Center" (Vec2 0 0) (so3ProjCoord p)
