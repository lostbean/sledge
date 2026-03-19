{-# LANGUAGE RecordWildCards #-}

module Main where

import Test.Tasty

import BinghamTest
import DDFTest
import IPFTest
import ODFTest
import SamplerTest
import SphereProjTest
import TesseractGridTest
import TestKernel
import TestTexture

main :: IO ()
main =
    defaultMain $
        testGroup
            "Tests"
            [ TestKernel.test
            , TestTexture.testTexture
            , TestTexture.testOrientation
            , BinghamTest.test
            , TesseractGridTest.test
            , SamplerTest.test
            , DDFTest.test
            , ODFTest.test
            , SphereProjTest.test
            , IPFTest.test
            ]
