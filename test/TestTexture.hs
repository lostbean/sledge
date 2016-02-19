{-# LANGUAGE RecordWildCards #-}

module TestTexture
  ( testTexture
  , testOrientation
  ) where


import Control.Applicative
import Data.Vector (Vector)
import Test.QuickCheck
import Test.Tasty
import Test.Tasty.HUnit      as HU
import Test.Tasty.QuickCheck as QC

import qualified Data.Vector as V
import Hammer.Math.Algebra
import Texture.Orientation
import Texture.Symmetry

-- ================================================================================

testTexture :: TestTree
testTexture = testGroup "Fundamental Zone"
  [ QC.testProperty "fundamental zone"         testFundamentalZone
  , QC.testProperty "fundamental zone(matrix)" testFundamentalZoneMatrix
  ]

instance Arbitrary Euler where
  arbitrary = liftA3 mkEuler x2 x1 x2
    where
      x1 = Deg <$> choose (0, 180)
      x2 = Deg <$> choose (0, 360)

instance Arbitrary Vec3 where
  arbitrary = normalize <$> liftA3 Vec3 p p p
    where p = choose (0,1)

instance Arbitrary Quaternion where
  arbitrary = toQuaternion <$> (arbitrary :: Gen Euler)

instance Arbitrary Rodrigues where
  arbitrary = fromQuaternion <$> arbitrary

msgFail :: (Show a, Testable prop)=> a -> prop -> Property
msgFail text = printTestCase ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")

errLimit :: Double
errLimit = 1e-7

testFundamentalZone :: Quaternion -> Property
testFundamentalZone q = let
  fz    = getSymmOps Cubic
  r     = fromQuaternion q :: Rodrigues
  inFZ1 = (toFZ Cubic q)   :: Quaternion
  inFZ2 = toFZDirect fz_r r
  fz_r  = V.map (fromQuaternion . symmOp) $ V.convert fz :: Vector Rodrigues
  delta = getOmega $ inFZ1 -#- (toQuaternion inFZ2)
  err   = abs delta < errLimit || abs (delta - 2*pi) < errLimit
  in msgFail ("No match in Fundamental Zone! " ++ show (delta)) err


testFundamentalZoneMatrix :: Quaternion -> Property
testFundamentalZoneMatrix q = let
  m     = fromQuaternion q :: RotMatrix
  inFZ1 = (toFZ Cubic q) :: Quaternion
  inFZ2 = toFZDirect mSymm m
  delta = getOmega $ inFZ1 -#- (toQuaternion inFZ2)
  err   = abs delta < errLimit || abs (delta - 2*pi) < errLimit

  e  = fromQuaternion q         :: Euler
  a1 = toAngle $ getOmega inFZ2 :: Deg
  a2 = toAngle $ getOmega inFZ1 :: Deg

  in msgFail ("No match in Fundamental Zone!" ++ show (e, a1, a2, delta)) err

mSymm :: Vector RotMatrix
mSymm = V.map func $ V.fromList ls
  where
    func (va, vb, vc) = mkUnsafeRotMatrix $
                        Mat3 (mkVec3 va) (mkVec3 vb) (mkVec3 vc)
    ls = [(( 1, 0, 0),
           ( 0, 1, 0),
           ( 0, 0, 1)),


          (( 1, 0, 0),
           ( 0, 0, 1),
           ( 0,-1, 0)),

          (( 1, 0, 0),
           ( 0,-1, 0),
           ( 0, 0,-1)),

          (( 1, 0, 0),
           ( 0, 0,-1),
           ( 0, 1, 0)),

          (( 0, 0,-1),
           ( 0, 1, 0),
           ( 1, 0, 0)),

          ((-1, 0, 0),
           ( 0, 1, 0),
           ( 0, 0,-1)),

          (( 0, 0, 1),
           ( 0, 1, 0),
           (-1, 0, 0)),

          (( 0, 1, 0),
           (-1, 0, 0),
           ( 0, 0, 1)),

          ((-1, 0, 0),
           ( 0,-1, 0),
           ( 0, 0, 1)),

          (( 0,-1, 0),
           ( 1, 0, 0),
           ( 0, 0, 1)),



          (( 0, 1, 0),
           ( 1, 0, 0),
           ( 0, 0,-1)),

          ((-1, 0, 0),
           ( 0, 0, 1),
           ( 0, 1, 0)),

          (( 0, 0, 1),
           ( 0,-1, 0),
           ( 1, 0, 0)),

          (( 0,-1, 0),
           (-1, 0, 0),
           ( 0, 0,-1)),

          ((-1, 0, 0),
           ( 0, 0,-1),
           ( 0,-1, 0)),

          (( 0, 0,-1),
           ( 0,-1, 0),
           (-1, 0, 0)),



          (( 0, 1, 0),
           ( 0, 0, 1),
           ( 1, 0, 0)),

          (( 0, 0, 1),
           ( 1, 0, 0),
           ( 0, 1, 0)),

          (( 0, 0,-1),
           (-1, 0, 0),
           ( 0, 1, 0)),

          (( 0,-1, 0),
           ( 0, 0, 1),
           (-1, 0, 0)),

          (( 0, 0, 1),
           (-1, 0, 0),
           ( 0,-1, 0)),

          (( 0,-1, 0),
           ( 0, 0,-1),
           ( 1, 0, 0)),

          (( 0, 0,-1),
           ( 1, 0, 0),
           ( 0,-1, 0)),

          (( 0, 1, 0),
           ( 0, 0,-1),
           (-1, 0, 0))

         ]

-- ================================== Test ======================================

testOrientation :: TestTree
testOrientation = testGroup "Orientations"
  [ QC.testProperty "deg/rad conversion"       testAngleConv
  , QC.testProperty "fundamental zone(matrix)" testFundamentalZoneMatrix
  , testOrientationModule
  , testConv
  ]

instance Arbitrary Deg where
  arbitrary = Deg <$> arbitrary

(=@=) :: Deg -> Deg -> Bool
a =@= b = abs (a - b) < Deg 0.0001

testAngleConv :: Deg -> Bool
testAngleConv an = an =@= (toAngle $ fromAngle an)

-- | Verify the correctness of this module
testOrientationModule :: TestTree
testOrientationModule = let
  testAngle =
    [ toQuaternion $ mkEuler (Deg   0.0)  (Deg (-190.0)) (Deg 0.0)
    , toQuaternion $ mkEuler (Deg 380.0)  (Deg      0.0) (Deg 0.0)
    ]
  in testCaseSteps "orientation" $ \step -> do
    step "testgetAngle id"
    mapM_ (testgetAngle id) testAngle
    step "testgetAngle getShortAngle"
    mapM_ (testgetAngle getShortAngle) testAngle
    step "testgetAngle getAbsShortAngle"
    mapM_ (testgetAngle getAbsShortAngle) testAngle
    step "composition"
    testComposition
    step "misorientation"
    testMisso

testConv :: TestTree
testConv = testGroup "Rotation angles"
  [ QC.testProperty "Matrix" $ \q -> let
    m   = (fromQuaternion q) :: RotMatrix
    qm  = toQuaternion m
    in Deg 0 =@= (toAngle $ getOmega (q -@- qm ))
  , QC.testProperty "Axis-pair" $ \q -> let
    ap  = (fromQuaternion q) :: AxisPair
    qap = toQuaternion ap
    in Deg 0 =@= (toAngle $ getOmega (q -@- qap))
  , QC.testProperty "Euler" $ \q -> let
    e   = (fromQuaternion q) :: Euler
    qe  = toQuaternion e
    in Deg 0 =@= (toAngle $ getOmega (q -@- qe ))
  , QC.testProperty "Rodrigues" $ \q -> let
    r   = (fromQuaternion q) :: Rodrigues
    qr  = toQuaternion r
    in Deg 0 =@= (toAngle $ getOmega (q -@- qr ))
  ]

testgetAngle :: (Double -> Double) -> Quaternion -> Assertion
testgetAngle foo q0 = do
  let
    toQ q = q                  :: Quaternion
    toA q = (fromQuaternion q) :: AxisPair
    toE q = (fromQuaternion q) :: Euler
    toR q = (fromQuaternion q) :: Rodrigues
    toM q = (fromQuaternion q) :: RotMatrix
    isOk a = a >= Deg 0 && a < Deg 360
    in do
      HU.assertBool "Quaternion" (isOk . toAngle . foo . getOmega . toQ $ q0)
      HU.assertBool "Matrix"     (isOk . toAngle . foo . getOmega . toM $ q0)
      HU.assertBool "Axis-pair"  (isOk . toAngle . foo . getOmega . toA $ q0)
      HU.assertBool "Euler"      (isOk . toAngle . foo . getOmega . toE $ q0)
      HU.assertBool "Rodrigues"  (isOk . toAngle . foo . getOmega . toR $ q0)

testComposition :: IO ()
testComposition = let
  e1 = mkEuler (Deg 30) (Deg 0 ) (Deg 0)
  e2 = mkEuler (Deg 0 ) (Deg 30) (Deg 0)
  e3 = mkEuler (Deg 30) (Deg 30) (Deg 0)
  q1 = toQuaternion e1
  q2 = toQuaternion e2
  m1 = fromQuaternion q1 :: RotMatrix
  m2 = fromQuaternion q2
  a1 = fromQuaternion q1 :: AxisPair
  a2 = fromQuaternion q2
  r1 = fromQuaternion q1 :: Rodrigues
  r2 = fromQuaternion q2
  ce = e1 #<= e2
  cq = q1 #<= q2
  cm = m1 #<= m2
  ca = a1 #<= a2
  cr = r1 #<= r2
  func = toAngle . getOmega . (e3 -@-) . fromQuaternion
  in do
    -- Rotation e1 followed by e2
    -- where e1 #<= e2 = Euler(30.0 deg, 30.0 deg, 0.0 deg)
    HU.assertBool "Euler"      (Deg 0 =@= (func $ toQuaternion ce))
    HU.assertBool "Quaternion" (Deg 0 =@= (func $ toQuaternion cq))
    HU.assertBool "RotMatrix"  (Deg 0 =@= (func $ toQuaternion cm))
    HU.assertBool "AxisPair"   (Deg 0 =@= (func $ toQuaternion ca))
    HU.assertBool "Rodrigues"  (Deg 0 =@= (func $ toQuaternion cr))

testMisso :: IO ()
testMisso = let
  e1 = mkEuler (Deg 90) (Deg 0)  (Deg 0)
  e2 = mkEuler (Deg 0)  (Deg 90) (Deg 0)
  q1 = toQuaternion e1
  q2 = toQuaternion e2
  m   = q1 -@- q2
  mi  = q2 -@- q1
  cm  = q1 -#- q2
  cmi = q2 -#- q1
  testm =  q1 #<= (q2 -#- q1)
  func a b = toAngle $ getOmega (a -@- b)
  in do
    HU.assertBool " e1 -@- e2 == e1 -#- e2"   (Deg 0 =@= (func (invert  m)  mi))
    HU.assertBool " e2 -@- e1 == e2 -#- e1"   (Deg 0 =@= (func (invert cm) cmi))
    HU.assertBool " e1 #<= (e2 -#- e1) == e2" (Deg 0 =@= (func testm q2))
