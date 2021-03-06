{-# LANGUAGE
    FlexibleInstances
  , ScopedTypeVariables
  , TypeSynonymInstances
  #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}
module TestTexture
  ( testTexture
  , testOrientation
  , testAverageQuaternion
  , testWeightedAverageQuaternion
  , testSymmAverageQuaternion
  , testRotProperties
  ) where

import Control.Applicative
import Data.Monoid ((<>), mempty)
import Data.Vector (Vector)
import Test.QuickCheck
import Test.Tasty
import Test.Tasty.HUnit      as HU
import Test.Tasty.QuickCheck as QC

import qualified Data.Vector as V
import Linear.Vect
import Linear.Mat
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

instance Arbitrary Vec3D where
  arbitrary = normalize <$> liftA3 Vec3 p p p
    where p = choose (0,1)

instance Arbitrary Quaternion where
  arbitrary = toQuaternion <$> (arbitrary :: Gen Euler)

instance Arbitrary Rodrigues where
  arbitrary = fromQuaternion <$> arbitrary

instance Arbitrary RotMatrix where
  arbitrary = fromQuaternion <$> arbitrary

instance Arbitrary AxisPair where
  arbitrary = fromQuaternion <$> arbitrary

instance Arbitrary Deg where
  arbitrary = Deg <$> arbitrary

(~=) :: (NearZero a, Ord a, Num a) => a -> a -> Bool
a ~= b = isMainlyZero $ a - b

msgFail :: (Show a, Testable prop)=> a -> prop -> Property
msgFail text = counterexample ("\x1b[7m Fail: " ++ show text ++ "! \x1b[0m")

omegaDeg :: Rot q => q -> q -> Deg
omegaDeg a b = toAngle . getOmega $ a -@- b :: Deg

omegaSymmDeg :: Rot q => Symm -> q -> q -> Deg
omegaSymmDeg symm a b = toAngle $ getMisoAngle symm (toQuaternion a) (toQuaternion b)

binaryTest :: (Show a, Show b) => String -> (a -> b -> Bool) -> a -> b -> Property
binaryTest msg test a b = counterexample (msg ++ ": " ++ show a ++ " " ++ show b) $ a `test` b

-- ================================================================================

testRotProperties :: TestTree
testRotProperties = testGroup "Orientations properties"
  [ testMonoidInstances
  , testGroupInstances
  , testRotInstances
  ]

testMonoidInstances :: TestTree
testMonoidInstances = testGroup "Monoid instances"
  [ QC.testProperty "Quaternion" (\a -> testMonoidClass (a :: Quaternion))
  , QC.testProperty "Euler"      (\a -> testMonoidClass (a :: Euler)     )
  , QC.testProperty "AxisPair"   (\a -> testMonoidClass (a :: AxisPair)  )
  , QC.testProperty "Rodrigues"  (\a -> testMonoidClass (a :: Rodrigues) )
  , QC.testProperty "RotMatrix"  (\a -> testMonoidClass (a :: RotMatrix) )
  ]

testMonoidClass :: Rot a => a -> a -> a -> Property
testMonoidClass q1 q2 q3 = let
  qa = (q1 <> q2) <> q3
  qb = q1 <> (q2 <> q3)
  qc = q1 <> mempty
  qd = mempty <> q1
  in binaryTest "associative" (~=) (omegaDeg qa qb) (Deg 0) .&&.
     binaryTest "right id"    (~=) (omegaDeg q1 qc) (Deg 0) .&&.
     binaryTest "left id"     (~=) (omegaDeg q1 qd) (Deg 0)

testGroupInstances :: TestTree
testGroupInstances = testGroup "Group instances"
  [ QC.testProperty "Quaternion" (\a -> testGroupClass (a :: Quaternion))
  , QC.testProperty "Euler"      (\a -> testGroupClass (a :: Euler)     )
  , QC.testProperty "AxisPair"   (\a -> testGroupClass (a :: AxisPair)  )
  , QC.testProperty "Rodrigues"  (\a -> testGroupClass (a :: Rodrigues) )
  , QC.testProperty "RotMatrix"  (\a -> testGroupClass (a :: RotMatrix) )
  ]

testGroupClass :: Rot a => a -> Property
testGroupClass q1 = let
  qa = q1 <> invert q1
  qb = invert q1 <> q1
  in binaryTest "right invert" (~=) (omegaDeg qa mempty) (Deg 0) .&&.
     binaryTest "left invert"  (~=) (omegaDeg qb mempty) (Deg 0)

testRotInstances :: TestTree
testRotInstances = testGroup "Rot instances"
  [ QC.testProperty "Quaternion" (\a -> testRotClass (a :: Quaternion))
  , QC.testProperty "Euler"      (\a -> testRotClass (a :: Euler)     )
  , QC.testProperty "AxisPair"   (\a -> testRotClass (a :: AxisPair)  )
  , QC.testProperty "Rodrigues"  (\a -> testRotClass (a :: Rodrigues) )
  , QC.testProperty "RotMatrix"  (\a -> testRotClass (a :: RotMatrix) )
  ]

testRotClass :: forall a . Rot a => a -> Property
testRotClass r1 = let
  q  = toQuaternion r1
  r2 = fromQuaternion q :: a
  in binaryTest "angle to quaternion"      (~=) (getOmega q )    (getOmega r1) .&&.
     binaryTest "angle between conversion" (~=) (getOmega r1)    (getOmega r2) .&&.
     binaryTest "conversion"               (~=) (omegaDeg r1 r2) (Deg 0)

-- ================================================================================

testFundamentalZone :: Quaternion -> Property
testFundamentalZone q = let
  fz    = getSymmOps Cubic
  r     = fromQuaternion q :: Rodrigues
  inFZ1 = toFZ Cubic     q :: Quaternion
  inFZ2 = toFZDirect fz_r r
  fz_r  = V.map (fromQuaternion . symmOp) $ V.convert fz :: Vector Rodrigues
  delta = getOmega $ inFZ1 -#- toQuaternion inFZ2
  err   = abs delta < epsilon || abs (delta - 2*pi) < epsilon
  in msgFail ("No match in Fundamental Zone! " ++ show delta) err

testFundamentalZoneMatrix :: Quaternion -> Property
testFundamentalZoneMatrix q = let
  m     = fromQuaternion q :: RotMatrix
  inFZ1 = toFZ Cubic     q :: Quaternion
  inFZ2 = toFZDirect mSymm m
  delta = getOmega $ inFZ1 -#- toQuaternion inFZ2
  err   = abs delta < epsilon || abs (delta - 2*pi) < epsilon

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
  [ testAngleConv
  , testOrientationModule
  ]

testAngleConv :: TestTree
testAngleConv = QC.testProperty "deg/rad conversion" $
  \an -> (an :: Deg) ~= toAngle (fromAngle an)

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

testgetAngle :: (Double -> Double) -> Quaternion -> Assertion
testgetAngle foo q0 = let
  toQ q = q                :: Quaternion
  toA q = fromQuaternion q :: AxisPair
  toE q = fromQuaternion q :: Euler
  toR q = fromQuaternion q :: Rodrigues
  toM q = fromQuaternion q :: RotMatrix
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
    HU.assertBool "Euler"      (Deg 0 ~= func (toQuaternion ce))
    HU.assertBool "Quaternion" (Deg 0 ~= func (toQuaternion cq))
    HU.assertBool "RotMatrix"  (Deg 0 ~= func (toQuaternion cm))
    HU.assertBool "AxisPair"   (Deg 0 ~= func (toQuaternion ca))
    HU.assertBool "Rodrigues"  (Deg 0 ~= func (toQuaternion cr))

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
    HU.assertBool " e1 -@- e2 == e1 -#- e2"   (Deg 0 ~= func (invert  m)  mi)
    HU.assertBool " e2 -@- e1 == e2 -#- e1"   (Deg 0 ~= func (invert cm) cmi)
    HU.assertBool " e1 #<= (e2 -#- e1) == e2" (Deg 0 ~= func testm q2)

-- ================================================================================

testAverageQuaternion :: TestTree
testAverageQuaternion = testGroup "Simple average of two quaternions"
  [ QC.testProperty "averageTwoQuaternion"      (testQAvg (\a b -> averageTwoQuaternion (1, a) (1, b)))
  , QC.testProperty "averageQuaternion"         (testQAvg (\a b -> averageQuaternion [a, b]))
  , QC.testProperty "averageWeightedQuaternion" (testQAvg (\a b -> averageWeightedQuaternion [(1, a), (1, b)]))
  ]

testQAvg :: (Quaternion -> Quaternion -> Quaternion) -> Quaternion -> Quaternion -> Property
testQAvg func q1 q2 = let
  avg = func q1 q2
  in binaryTest "same angular dist" (~=) (omegaDeg avg q1) (omegaDeg avg q2)

-- ================================================================================

testWeightedAverageQuaternion :: TestTree
testWeightedAverageQuaternion = testGroup "Weighted average of two quaternions"
  [ QC.testProperty "averageTwoQuaternion"      (testWQAvg averageTwoQuaternion)
  , QC.testProperty "averageWeightedQuaternion" (testWQAvg (\a b -> averageWeightedQuaternion [a, b]))
  ]

testWQAvg :: ((Double, Quaternion) -> (Double, Quaternion) -> Quaternion) -> Quaternion -> Quaternion -> Property
testWQAvg func q1 q2 = let
  avg1 = func (1,q1) (0,q2)
  avg2 = func (0,q1) (1,q2)
  avg3 = func (2,q1) (1,q2)
  in binaryTest "close to q1"  (~=) (omegaDeg avg1 q1) (Deg 0) .&&.
     binaryTest "close to q2"  (~=) (omegaDeg avg2 q2) (Deg 0) .&&.
     binaryTest "closer to q1" (<)   (omegaDeg avg3 q1) (omegaDeg avg3 q2)

-- ================================================================================

testSymmAverageQuaternion :: TestTree
testSymmAverageQuaternion = testGroup "Simple average of two quaternions with symmetry"
  [ QC.testProperty "averageQuaternionWithSymm" (testSymmQAvg (\a b -> averageQuaternionWithSymm Cubic [a, b]))
  ]

testSymmQAvg :: (Quaternion -> Quaternion -> Quaternion) -> Quaternion -> Quaternion -> Property
testSymmQAvg func q1 q2 = let
  avg = func q1 q2
  in binaryTest "same angular dist" (~=) (omegaSymmDeg Cubic avg q1) (omegaSymmDeg Cubic avg q2)
