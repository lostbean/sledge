module SamplerTest (
    test,
) where

import Linear.Vect
import Test.Tasty
import Test.Tasty.HUnit
import Texture.Bingham
import Texture.Orientation
import Texture.Sampler

test :: TestTree
test =
    testGroup
        "Sampler"
        [ testCase "Hit and Run Slice Sampler (Double)" testHRS_Double
        , testCase "Hit and Run Slice Sampler (Quaternion)" testHRS_Quaternion
        ]

testHRS_Double :: Assertion
testHRS_Double = do
    let
        -- Normal distribution (not normalized, but it shouldn't matter for the sampler)
        p x = exp (-(x * x) / 2)
        cfg = defaultCfg{initShotDist = 1.0, maxShotDist = 10.0}
    samples <- hitAndRunSlice cfg p 0 100
    assertBool "Got samples" (length samples == 100)

testHRS_Quaternion :: Assertion
testHRS_Quaternion = do
    let
        d1 = (10, mkQuaternion (Vec4 1 0 0 0))
        d2 = (5, mkQuaternion (Vec4 0 1 0 0))
        d3 = (2, mkQuaternion (Vec4 0 0 1 0))
        dist = mkBingham d1 d2 d3
        p = binghamPDF dist
        cfg = defaultCfg
    samples <- hitAndRunSlice cfg p (bingMode dist) 100
    assertBool "Got samples" (length samples == 100)
