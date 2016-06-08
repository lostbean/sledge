{-# LANGUAGE RecordWildCards #-}
module Main where

import Test.Tasty

import TestTexture
import TestKernel

main :: IO ()
main = defaultMain
     $ testGroup "Tests"
     [ TestKernel.test
     , TestTexture.testTexture
     , TestTexture.testOrientation
     , TestTexture.testAverageQuaternion
     , TestTexture.testWeightedAverageQuaternion
     ]
