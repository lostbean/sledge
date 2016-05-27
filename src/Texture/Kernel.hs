{-# LANGUAGE FlexibleContexts #-}
module Texture.Kernel
  ( gaussianKernel
  , addOneKernel
  , addManyKernels
  , addOneKernelWithConst
  , addManyKernelsWithConst
  ) where

import Control.Monad.ST
import Data.BlazeVPtree
import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Texture.Orientation

gaussianKernel :: (Angle a)=> a -> Double -> Double
gaussianKernel a x = exp (k * x * x)
  where
    k = -0.5 / (w * w)
    w = fromAngle a

addKernelM :: (Metric p, U.Unbox p, Angle a, GM.MVector v Double) => a -> VPtree p -> p -> v s Double -> ST s ()
addKernelM a tree p acc = mapM_ func fs
  where
    fs = nearNeighbors tree (3 * fromAngle a) p
    func (i, _, d) = do
      old <- GM.read acc i
      GM.write acc i (old + gaussianKernel a d)

addOneKernel :: (Metric p, U.Unbox p, Angle a, G.Vector v Double) => a -> VPtree p -> p -> v Double -> v Double
addOneKernel a tree p = G.modify (addKernelM a tree p)

addManyKernels :: (Metric p, U.Unbox p, Angle a, G.Vector v p, G.Vector v Double) => a -> VPtree p -> v p -> v Double -> v Double
addManyKernels a tree ps = G.modify (\v -> G.mapM_ (\p -> addKernelM a tree p v) ps)


addKernelWithConstM :: (Metric p, U.Unbox p, Angle a, GM.MVector v Double) => a -> Double -> VPtree p -> p -> v s Double -> ST s ()
addKernelWithConstM a k tree p acc = mapM_ func fs
  where
    fs = nearNeighbors tree (3 * fromAngle a) p
    func (i, _, d) = do
      old <- GM.read acc i
      GM.write acc i (old + k * gaussianKernel a d)


addOneKernelWithConst :: ( Metric p
                         , U.Unbox p
                         , Angle a
                         , G.Vector v Double
                         ) => a -> Double -> VPtree p -> p -> v Double -> v Double
addOneKernelWithConst  a k tree p = G.modify (addKernelWithConstM a k tree p)


addManyKernelsWithConst :: ( Metric p
                           , U.Unbox p
                           , Angle a
                           , G.Vector v p
                           , G.Vector v Double
                           ) => a -> Double -> VPtree p -> v p -> v Double -> v Double
addManyKernelsWithConst a k tree ps = G.modify (\v -> G.mapM_ (\p -> addKernelWithConstM a k tree p v) ps)
