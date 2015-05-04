{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE FlexibleContexts #-}

module Texture.Kernel
       ( gaussianKernel
       , addOneKernel
       , addManyKernels
       , removeOneKernel
       , removeManyKernels
       ) where

import qualified Data.Vector.Generic.Mutable as GM
import qualified Data.Vector.Generic as G
import qualified Data.Vector.Unboxed as U

import Control.Monad.ST

import Data.BlazeVPtree

import Texture.Orientation

gaussianKernel :: (Angle a)=> a -> Double -> Double
gaussianKernel a x = exp (k * x * x)
  where
    k = -0.5 / (w * w)
    w = fromAngle a

addKernelM :: (Metric p, U.Unbox p, Angle a, GM.MVector v Double)
           => a -> VPtree p -> p -> v s Double -> ST s ()
addKernelM a tree p acc = mapM_ func fs
  where
    fs = nearNeighbors tree (3 * fromAngle a) p
    func (i, q, d) = do
      old <- GM.read acc i
      GM.write acc i (old + gaussianKernel a d)

addOneKernel :: (Metric p, U.Unbox p, Angle a, G.Vector v Double)
             => a -> VPtree p -> p -> v Double -> v Double
addOneKernel a tree p = G.modify (addKernelM a tree p)

addManyKernels :: (Metric p, U.Unbox p, Angle a, G.Vector v p, G.Vector v Double)
               => a -> VPtree p -> v p -> v Double -> v Double
addManyKernels a tree ps = G.modify (\v -> G.mapM_ (\p -> addKernelM a tree p v) ps)


removeKernelM :: (Metric p, U.Unbox p, Angle a, GM.MVector v Double)
              => a -> VPtree p -> p -> v s Double -> ST s ()
removeKernelM a tree p acc = mapM_ func fs
  where
    fs = nearNeighbors tree (3 * fromAngle a) p
    func (i, q, d) = do
      old <- GM.read acc i
      GM.write acc i (old - gaussianKernel a d)

removeOneKernel :: (Metric p, U.Unbox p, Angle a, G.Vector v p, G.Vector v Double)
                => a -> VPtree p -> p -> v Double -> v Double
removeOneKernel a tree p = G.modify (removeKernelM a tree p)

removeManyKernels :: (Metric p, U.Unbox p, Angle a, G.Vector v p, G.Vector v Double)
                  => a -> VPtree p -> v p -> v Double -> v Double
removeManyKernels a tree ps = G.modify (\v -> G.mapM_ (\p -> removeKernelM a tree p v) ps)
