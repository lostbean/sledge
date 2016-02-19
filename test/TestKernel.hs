module TestKernel
  ( test
  ) where

import Data.Function (on)
import System.Random
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector.Unboxed as U
import qualified Data.List           as L

import Hammer.Math.Algebra
import Texture.Orientation
import Texture.Symmetry
import Texture.TesseractGrid
import Texture.Kernel
import Texture.ODF
import qualified Data.BlazeVPtree as VP

test :: TestTree
test = testGroup "Kernel" [testKernel]

-- | Generate random values of Rodrigues-Frank C4 symmetry (90 <100>, 90 <010> 90 <001>)
-- with FZ planes at (+-45 <100>,+-45 <010>,+-45 <001>)
getFRFZ :: IO [Quaternion]
getFRFZ = newStdGen >>= \g -> do
  return $ map (toFZ Cubic) (randoms g)

-- | Test VP tree on Rodrigues-Frank space
testKernel = let
  rs = U.filter (isInRodriFZ Cubic) $ genQuaternionGrid 80
  vp = VP.fromVector rs
  t = toQuaternion $ mkAxisPair (Vec3 1 0 0) (Deg 44)
  d = fromAngle (Deg 7.5)

  vpo = VP.nearNeighbors vp d t
  vplist = L.sort $ map (\(i,_,_) -> i) vpo

  bfos = U.filter ((< d) . VP.dist t . snd) (U.imap (,) rs)
  nbfs@(ibfs,_) = U.minimumBy (compare `on` (VP.dist t . snd)) bfos
  bflists = L.sort $ map fst $ U.toList bfos

  in testCase "check nears" $ bflists @=? vplist
