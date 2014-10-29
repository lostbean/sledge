{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}

module Texture.Kernel
       ( gaussianKernel
       , addOneKernel
       , addManyKernels
       ) where

import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Generic as G

import Control.Monad.ST
import Control.Monad.Primitive

import Data.BlazeVPtree

import Texture.Orientation

gaussianKernel :: (Angle a)=> a -> Double -> Double
gaussianKernel a x = exp (k * x * x)
  where
    k = -0.5 / (w * w)
    w = fromAngle a

addOneKernel :: (Metric Quaternion, Angle a, G.Vector v Double)
             => a -> VPtree Quaternion -> Quaternion -> v Double -> v Double
addOneKernel a tree p = G.modify (addKernelM a tree p)

addManyKernels :: (Metric Quaternion, Angle a, G.Vector v Quaternion, G.Vector v Double)
               => a -> VPtree Quaternion -> v Quaternion -> v Double -> v Double
addManyKernels a tree ps = G.modify (\v -> G.mapM_ (\p -> addKernelM a tree p v) ps)

addKernelM :: (Metric Quaternion, Angle a, GM.MVector v Double, PrimMonad m)
           => a -> VPtree Quaternion -> Quaternion -> v (PrimState m) Double -> m ()
addKernelM a tree p acc = mapM_ func fs
  where
    fs = nearNeighbors tree (3 * fromAngle a) p
    fuu x = acos . composeQ0 x . invert
    func (i, q, d) = do
      old <- GM.read acc i
      GM.write acc i (old + gaussianKernel a d)
