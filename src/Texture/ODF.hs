{-# LANGUAGE RecordWildCards #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
module Texture.ODF
  ( ODF ( odfIntensity
        , odfGrid
        , odfGridSize
        , odfGridStep
        , odfTree
        , odfSymm
        , odfKernelWidth
        )
  , buildEmptyODF
  , resetODF
  , addPoints
  , addPointsParallel
  , addPointsWithConst
  , integrateODFwith
  , getODFeval
  , getMaxOrientation
  , renderODFVTK
  ) where

import Data.Maybe
import Hammer.Parallel.Vector (autoParFilter, runEval)
import Hammer.VTK
import Linear.Vect
import qualified Data.BlazeVPtree as VP
import qualified Data.Vector.Unboxed as U

import Texture.Orientation
import Texture.Symmetry
import Texture.IsoSphere
import Texture.TesseractGrid
import Texture.Kernel

instance VP.Metric Quaternion where
  dist = getMisoAngle Cubic

data ODF
  = ODF
  { odfIntensity   :: U.Vector Double
  , odfGrid        :: U.Vector Quaternion
  , odfGridSize    :: Int
  , odfGridStep    :: Int
  , odfTree        :: VP.VPtree Quaternion
  , odfSymm        :: Symm
  , odfKernelWidth :: Rad
  } deriving (Show)


buildEmptyODF :: Deg -> Symm -> Deg -> ODF
buildEmptyODF kw symm step = ODF
  { odfIntensity   = U.replicate n 0
  , odfGrid        = qs
  , odfGridSize    = n
  , odfGridStep    = s
  , odfTree        = VP.fromVector qs
  , odfSymm        = symm
  , odfKernelWidth = toAngle (fromAngle kw)
  }
  where
    n  = U.length qs
    s  = abs $ round (4 / fromAngle step)
    qs = runEval $ autoParFilter (isInRodriFZ symm) (genQuaternionGrid s)

resetODF :: ODF -> ODF
resetODF odf = odf {odfIntensity = U.replicate (odfGridSize odf) 0}

addPointsWithConst :: U.Vector Quaternion -> Double -> Maybe Rad -> ODF -> ODF
addPointsWithConst qs k customWidth odf@ODF{..} = odf { odfIntensity = is }
  where is    = addManyKernelsWithConst width k odfTree qs odfIntensity
        width = fromMaybe odfKernelWidth customWidth

addPoints :: U.Vector Quaternion -> ODF -> ODF
addPoints qs odf@ODF{..} = odf { odfIntensity = is }
  where is = addManyKernels odfKernelWidth odfTree qs odfIntensity

addPointsParallel :: U.Vector Quaternion -> ODF -> ODF
addPointsParallel qs odf@ODF{..} = odf { odfIntensity = is }
  where is = addManyKernelsParallel odfKernelWidth odfTree qs odfIntensity

getODFeval :: ODF -> (Quaternion -> Double)
getODFeval ODF{..} = maybe 0 (\(i, _, _) -> odfIntensity U.! i) . func
  where
    step = 4 / fromIntegral odfGridStep
    func = VP.nearestThanNeighbor odfTree step

getMaxOrientation :: ODF -> (Quaternion, Double)
getMaxOrientation ODF{..} = (odfGrid U.! i, odfIntensity U.! i)
  where i = U.maxIndex odfIntensity

integrateODFwith :: (U.Unbox a, Num a)=> (Quaternion -> Double -> a) -> ODF -> a
integrateODFwith func ODF{..} = U.foldl' (+) 0 $ U.zipWith func odfGrid odfIntensity

-- | Render ODF
renderODFVTK :: ODF -> VTK Vec3D
renderODFVTK ODF{..} = let
  attr = mkPointValueAttr "Intensity" (\i _ -> odfIntensity U.! i)
  vtk  = renderQuaternions odfGrid []
  in addPointValueAttr vtk attr
