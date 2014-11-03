{-# LANGUAGE RecordWildCards #-}

module Texture.DDF
       ( DDF ( ddfIntensity
             , ddfGrid
             , ddfGridSize
             , ddfGridStep
             , ddfTree
             , ddfSymm
             , ddfKernelWidth
             )
       , buildEmptyDDF
       , resetDDF
       , addPoints
       , getDDFeval
       , renderDDFVTK
       ) where

import qualified Data.Vector.Unboxed as U

import Hammer.Math.Algebra
import Hammer.VTK

import qualified Data.BlazeVPtree as VP

import Texture.Orientation
import Texture.Symmetry
import Texture.IsoSphere
import Texture.Kernel

instance VP.Metric Normal3 where
  dist v1 v2 = acosSafe (v1 &. v2)

data DDF
  = DDF
  { ddfIntensity   :: U.Vector Double
  , ddfGrid        :: U.Vector Normal3
  , ddfGridSize    :: Int
  , ddfGridStep    :: Rad
  , ddfTree        :: VP.VPtree Normal3
  , ddfSymm        :: Symm
  , ddfKernelWidth :: Rad
  } deriving (Show)


buildEmptyDDF :: Deg -> Symm -> Deg -> DDF
buildEmptyDDF kw symm step
  = DDF
 { ddfIntensity   = U.replicate n 0
 , ddfGrid        = qs
 , ddfGridSize    = n
 , ddfGridStep    = w
 , ddfTree        = VP.fromVector qs
 , ddfSymm        = symm
 , ddfKernelWidth = toAngle (fromAngle kw)
 }
 where
   (s, w) = getOptSubDivisionLevel step
   n  = U.length qs
   qs = U.map mkNormal $ vertices $ isoSphere s

resetDDF :: DDF -> DDF
resetDDF ddf = ddf { ddfIntensity = U.replicate (ddfGridSize ddf) 0}

addPoints :: U.Vector Normal3 -> DDF -> DDF
addPoints qs ddf@DDF{..} = ddf { ddfIntensity = is }
  where is = addManyKernels ddfKernelWidth ddfTree qs ddfIntensity

getDDFeval :: DDF -> (Normal3 -> Double)
getDDFeval DDF{..} = maybe 0 (\((i,_,_)) -> ddfIntensity U.! i) . func
  where
    step = fromAngle ddfGridStep
    func = VP.nearestThanNeighbor ddfTree step

-- | Render DDF
renderDDFVTK :: DDF -> VTK Normal3
renderDDFVTK DDF{..} = let
  attr = mkPointValueAttr "Intensity" (\i _ -> ddfIntensity U.! i)
  in mkUGVTK "DDF" ddfGrid (U.generate ddfGridSize id) [attr] []
