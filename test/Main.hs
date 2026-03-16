{-# LANGUAGE RecordWildCards #-}

module Main where

import Test.Tasty

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
            ]
