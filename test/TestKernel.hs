module TestKernel (
    test,
) where

import qualified Data.List as L
import qualified Data.Vector.Unboxed as U
import Test.Tasty
import Test.Tasty.HUnit

import qualified Data.BlazeVPtree as VP
import Linear.Vect
import Texture.ODF ()
import Texture.Orientation
import Texture.Symmetry
import Texture.TesseractGrid

test :: TestTree
test = testGroup "Kernel" [testKernel]

-- | Test VP tree on Rodrigues-Frank space
testKernel :: TestTree
testKernel =
    let
        rs = U.filter (isInRodriFZ Cubic) $ genQuaternionGrid 80
        vp = VP.fromVector rs
        t = toQuaternion $ mkAxisPair (Vec3 1 0 0) (Deg 44)
        d = fromAngle (Deg 7.5)

        vpo = VP.nearNeighbors vp d t
        vplist = L.sort $ map (\(i, _, _) -> i) vpo

        bfos = U.filter ((< d) . VP.dist t . snd) (U.imap (,) rs)
        bflists = L.sort $ map fst $ U.toList bfos
     in
        testCase "check nears" $ bflists @=? vplist
