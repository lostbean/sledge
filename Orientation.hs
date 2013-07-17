{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

-- |
-- Module      : Hammer.Texture.Orientation
-- Copyright   : (c) 2013 Edgar Gomes
-- License     : Privete-style (see LICENSE)
--
-- Maintainer  : Edgar Gomes <talktoedgar@gmail.com>
-- Stability   : experimental
-- Portability : tested on GHC only
--
-- This module defines orientation representations and its operations.
-- The general convention is
--
--  * e1 (1,0,0) -> RD
--
--  * e2 (0,1,0) -> TD
--
--  * e3 (0,0,1) -> ND
--
-- References:
--
--  * Conversion of EBSD data by a quaternion based algorithm to be used for grain
-- structure simulations.
--
--  * Orientation Library Manual by Romain Quey
--
module Hammer.Texture.Orientation
       (  -- * Reference Frame
         RefFrame (..)
          -- * Angles
       , Angle      (..)
       , Deg        (..)
       , Rad        (..)
       , toDeg
         -- * Rotation class
       , Rot        (..)
       , (=>#)
         -- * Quaternions
       , Quaternion (quaterVec)
       , mkUnsafeQuaternion
       , mkQuaternion
       , splitQuaternion
       , unsafeMergeQuaternion
       , mergeQuaternion
       , antipodal
         -- * Euler
       , Euler      (phi1, phi, phi2)
       , mkEuler
         -- * Axis-Angle
       , AxisPair   (axisAngle)
       , mkAxisPair
         -- * Frank-Rodrigues
       , Rodrigues  (rodriVec)
       , mkRodrigues
         -- * Matrix
       , RotMatrix  (rotMat)
       , mkUnsafeRotMatrix
         -- * Vector rotation
       , activeVecRotation
       , passiveVecRotation
         -- * Other functions
       , quaterInterpolation
       , aproxToIdealAxis
       ) where

import Data.Ratio
import Numeric

import Hammer.Math.Algebra
import System.Random

-- import Debug.Trace
-- dbg s x = trace (s ++ show x) x

-- ==================================  Types ====================================

-- | External reference frame.
--
--  * RD -> x direction
--
--  * TD -> y direction
--
--  * ND -> z direction
--
data RefFrame = ND
              | TD
              | RD
              deriving (Show, Eq)

-- | Unit quaternion representation of orientation (rotation).
-- It is the basic and central representation used in this module.
newtype Quaternion =
  Quaternion
  { quaterVec :: Vec4
  } deriving (Eq)

instance Random Quaternion where
  random             = randomR (zerorot, zerorot)
  randomR (_, _) gen = let
    (Vec3 x y z, gen1) = randomR (zero, Vec3 1 1 1) gen
    -- Following : K. Shoemake., Uniform random rotations.
    -- In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992.
    q0 = sqrt(1-x) * sin (2 * pi * y)
    q1 = sqrt(1-x) * cos (2 * pi * y)
    q2 = sqrt(x)   * sin (2 * pi * z)
    q3 = sqrt(x)   * cos (2 * pi * z)
    in (mkQuaternion $ Vec4 q0 q1 q2 q3, gen1)

-- | Angles in degrees.
newtype Deg = Deg { unDeg :: Double } deriving (Eq, Num, Ord)

-- | Angles in radians.
newtype Rad = Rad { unRad :: Double } deriving (Eq, Num, Ord)

-- | Euler angles represents the composition of three independent
-- rotations: phi1 -> PHI -> phi2 where
--
--  * phi1 ~ ND   (Z)
--
--  * PHI  ~ RD\'  (X\')
--
--  * phi2 ~ ND\'\' (Z\'\')
--
data Euler =
  Euler
  { phi1 :: Double
  , phi  :: Double
  , phi2 :: Double
  } deriving (Eq)

-- | Axis-angle representation. The axis is a normalized direction
-- and the angle in radians.
newtype AxisPair =
  AxisPair
  { axisAngle :: (Vec3, Double)
  } deriving (Eq)

-- | Frank-Rodrigues representation.
newtype Rodrigues =
  Rodrigues
  { rodriVec :: Vec3
  } deriving (Eq)

-- | Frank-Rodrigues representation.
newtype RotMatrix =
  RotMatrix
  { rotMat :: Mat3
  } deriving (MultSemiGroup, Matrix, Eq)

-- ==================================  Angle class ===================================

class Angle a where
  -- | Converts a value in /radians/ to its 'Angle' format.
  toAngle   :: Double -> a
  -- | Converts a 'Angle' to its correspondent value in *radians*.
  fromAngle :: a -> Double

instance Angle Deg where
  toAngle   = Deg . (* (180 / pi))
  fromAngle = (* (pi / 180)) . unDeg

instance Angle Rad where
  toAngle   = Rad
  fromAngle = unRad

-- | Converts a value in radians to 'Deg'. It is the same as 'toAngle' for a 'Deg' type.
toDeg :: Double -> Deg
toDeg = toAngle

-- ===============================  Data constructors =================================

-- | Safe constructor for axis-angle.
mkAxisPair :: (Angle ang)=> Vec3 -> ang -> AxisPair
mkAxisPair r omega
  | rlen > 0 = AxisPair ((1 / rlen) *& r, fromAngle omega)
  | otherwise = zerorot
  where rlen = norm r

-- | Safe constructor for Euler angles in the Bunge convension.
-- Its sequence is: phi1 -> PHI -> phi2
mkEuler :: (Angle a1, Angle a2, Angle a3)=> a1 -> a2 -> a3 -> Euler
mkEuler p1 p p2 = Euler
  { phi1 = fromAngle p1
  , phi  = fromAngle p
  , phi2 = fromAngle p2 }

-- | Safe constructor for quaternions
mkQuaternion :: Vec4 -> Quaternion
mkQuaternion = Quaternion . normalize

-- | Constructor for Frank-Rodrigues.
mkRodrigues :: (Angle ang)=> Vec3 -> ang -> Rodrigues
mkRodrigues v omega
  | vlen > 0  = Rodrigues $ v &* (k / vlen)
  | otherwise = zerorot
  where
    vlen = norm v
    k    = tan ((fromAngle omega) / 2)

-- | Splits the quaternion in two components: angle related and direction related.
splitQuaternion :: Quaternion -> (Double, Vec3)
splitQuaternion (Quaternion (Vec4 q0 q1 q2 q3)) = (q0, Vec3 q1 q2 q3)

-- | Merges the two quaternion components: angle related and direction related.
mergeQuaternion :: (Double, Vec3) -> Quaternion
mergeQuaternion (q0, Vec3 q1 q2 q3) = mkQuaternion (Vec4 q0 q1 q2 q3)

-- | /Unsafe/ quaternion constructor. It assumes a normalized input.
mkUnsafeQuaternion :: Vec4 -> Quaternion
mkUnsafeQuaternion = Quaternion

-- | /Unsafe/ rotation matrix constructor. It assumes a orthonormal matrix with
-- unitary vector as input.
mkUnsafeRotMatrix :: Mat3 -> RotMatrix
mkUnsafeRotMatrix = RotMatrix

-- | /Unsafe/ quaternion merging constructor. It assumes a normalized input where
-- the single value and the vector in @R3@ form a unit vector in @R4@ (unit quaternion).
unsafeMergeQuaternion :: (Double, Vec3) -> Quaternion
unsafeMergeQuaternion (q0, Vec3 q1 q2 q3) = mkUnsafeQuaternion (Vec4 q0 q1 q2 q3)

-- =====================================  Rot class ======================================

-- | This class defines basic orientation (rotation) operations.
class Rot a where
  -- | Orientation composition. Given two orientations (passive rotations)
  -- @g1@ and @g2@, the operation
  --
  --  >  g3 = g1 #<= g2
  --
  -- applies a passive rotation @g2@ over the @g1@. In other words, it is the
  -- rotation @g1@ followed by the rotation @g2@. Note that the @g2@ rotation is
  -- applied on the reference frame of @g1@. This operation is not commutative.
  (#<=) :: a -> a -> a

  -- | Misorientation operator in the crystal frame. Given two orientations
  -- (passive rotations) g1 and g2, the operation
  --
  --  > g12 = g2 -#- g1
  --
  -- finds the rotation difference between g1 and g2 in the /reference frame
  -- of g1/ that brings g1 to coincide to g2 . Also note following rules:
  --
  --  >  g2 -#- g1 == invert (g1 -#- g2)
  --  >  g2 == g1 #<= (g2 -#- g1)
  --
  (-#-) :: a -> a -> a
  a -#- b = (invert b) #<= a

  -- | Misorientation operator in the sample reference frame. Given two orientations
  -- (passive rotations) @g1@ and @g2@, the operation
  --
  --  >  g12 = g2 -@- g1
  --
  -- finds the rotation difference between @g1@ and @g2@ in the reference global
  -- (sample) frame that brings @g1@ to coincide to @g2@. Also note following rule:
  --
  --  >  g2 -@- g1 == invert (g1 -@- g2)
  --
  (-@-) :: a -> a -> a
  a -@- b = a #<= (invert b)

  -- | Opposite orientation (inverse rotation).
  invert :: a -> a

  -- | Identity rotation.
  zerorot :: a

  -- | Converts a rotation to its quaternion form.
  toQuaternion :: a -> Quaternion

  -- | Converts a quaternion to another representation form.
  fromQuaternion :: Quaternion -> a

  -- | Get rotation angle in radians. Normally from [0 .. 2*pi].
  getOmega :: a -> Double

-- | The same as in '#<=' but with inverted parameters. Therefore, @g3 = g2 =># g1@ where
-- @g1@ is applied over @g2@ in the reference frame of @g2@ itself, resulting in @g3@.
(=>#) :: (Rot a)=> a -> a -> a
(=>#) = flip (#<=)

instance Rot Quaternion where
  p #<= q = let
    (p0, pv) = splitQuaternion p
    (q0, qv) = splitQuaternion q
    pq0 = p0 * q0 - pv &. qv
    pq  = p0 *& qv &+ q0 *& pv &+ pv &^ qv
    in unsafeMergeQuaternion (pq0, pq)
  invert q = let
    (q0, qv) = splitQuaternion q
    in unsafeMergeQuaternion (q0, neg qv)
  toQuaternion   = id
  fromQuaternion = id
  zerorot        = mkQuaternion $ Vec4 1 0 0 0
  getOmega       = (2 *) . acosSafe . fst . splitQuaternion

instance Rot Rodrigues where
  (Rodrigues rA) #<= (Rodrigues rB) = let
    r = (rB &+ rA) &- (rB &^ rA)
    k = 1.0 - (rB &. rA)
    in Rodrigues $ r &* (1 / k)
  invert (Rodrigues r) = Rodrigues $ neg r
  zerorot  = Rodrigues zero
  getOmega = (2 *) . atan . norm . rodriVec
  toQuaternion (Rodrigues r)
    | rlen > 0  = unsafeMergeQuaternion (q0, qv)
    | otherwise = zerorot
    where
      rlen  = norm r
      halfW = atan rlen  -- half omega
      q0    = cos halfW
      qv    = r &* ((sin halfW) / rlen)
  fromQuaternion q
    | qvlen > 0 = Rodrigues (k *& qv)
    | otherwise = zerorot
    where
      (q0, qv) = splitQuaternion q
      halfW    = acosSafe q0  -- half omega
      qvlen    = norm qv
      k        = (tan halfW) / qvlen

instance Rot AxisPair where
  a #<= b = fromQuaternion $ (toQuaternion a) #<= (toQuaternion b)
  invert (AxisPair (v, a)) = AxisPair (neg v, a)
  zerorot  = AxisPair (Vec3 1 0 0, 0)
  getOmega = snd . axisAngle
  toQuaternion (AxisPair (r, omega))= let
    q0 = cos (0.5 * omega)
    q  = r &* sin (0.5 * omega)
    in unsafeMergeQuaternion (q0, q)
  fromQuaternion q
    | qvlen > 0 = AxisPair ((1 / qvlen) *& qv, omega)
    | otherwise = zerorot
    where
      (q0, qv) = splitQuaternion q
      omega    = 2 * acosSafe q0
      qvlen    = norm qv

instance Rot Euler where
  a #<= b  = fromQuaternion $ (toQuaternion a) #<= (toQuaternion b)
  invert   = fromQuaternion . invert . toQuaternion
  zerorot  = Euler 0 0 0
  getOmega = getOmega . toQuaternion
  toQuaternion Euler{..} = let
    dphi = (phi1 - phi2) / 2  -- the reference was wrong, phi1 and phi2 were swapped.
    aphi = (phi1 + phi2) / 2
    s  = sin (phi / 2)
    c  = cos (phi / 2)
    q0 = c * (cos aphi)
    q1 = s * (cos dphi)
    q2 = s * (sin dphi)
    q3 = c * (sin aphi)
    in Quaternion (Vec4 q0 q1 q2 q3)
  fromQuaternion (Quaternion (Vec4 q0 q1 q2 q3))
    | phi == 0  = Euler (2 * acosSafe q0) phi 0 -- zero by convesion
    | phi == pi = Euler (2 * acosSafe q1) phi 0 -- zero by convesion
    | otherwise = Euler phi1          phi phi2
    where
      -- Uses atan2  to avoid problem with div by zero and zero/zero
      phi1 = atan2 q3 q0 + atan2 q2 q1
      phi  = 2 * atan2 (sqrt (q1 * q1 + q2 * q2)) (sqrt (q0 * q0 + q3 * q3))
      phi2 = atan2 q3 q0 - atan2 q2 q1

instance Rot RotMatrix where
  a #<= b  = b .*. a  -- In matrix composition the multiplication is commuted
  invert   = transpose
  zerorot  = idmtx
  getOmega = acosSafe . (\x -> 0.5 * (x - 1)). vecFoldr (+) . diagVec . rotMat
  toQuaternion = mat2quat
  fromQuaternion r = let
    Vec4 q0 q1 q2 q3 = quaterVec r

    g11 = q0*q0 + q1*q1 - q2*q2 - q3*q3
    g12 = 2 * (q1 * q2 + q0 * q3)
    g13 = 2 * (q1 * q3 - q0 * q2)
    g21 = 2 * (q1 * q2 - q0 * q3)
    g22 = q0*q0 - q1*q1 + q2*q2 - q3*q3
    g23 = 2 * (q2 * q3 + q0 * q1)
    g31 = 2 * (q1 * q3 + q0 * q2)
    g32 = 2 * (q2 * q3 - q0 * q1)
    g33 = q0*q0 - q1*q1 - q2*q2 + q3*q3

    in RotMatrix (Mat3 (Vec3 g11 g12 g13) (Vec3 g21 g22 g23) (Vec3 g31 g32 g33))

mat2quat :: RotMatrix -> Quaternion
mat2quat (RotMatrix m)
  | q0 > 0 = let
    q1 = (g23 - g32) / (4 * q0)
    q2 = (g31 - g13) / (4 * q0)
    q3 = (g12 - g21) / (4 * q0)
    in mkUnsafeQuaternion $ Vec4 q0 q1 q2 q3
  | q1_0 > q2_0 && q1_0 > q3_0 = let
    v = Vec4 q0 q1_0 (sgn g21 q2_0) (sgn g31 q3_0)
    in mkUnsafeQuaternion v
  | q2_0 > q1_0 && q2_0 > q3_0 = let
    v = Vec4 q0 (sgn g12 q1_0) q2_0 (sgn g32 q3_0)
    in mkUnsafeQuaternion v
  | otherwise = let
    v = Vec4 q0 (sgn g13 q1_0) (sgn g23 q2_0) q3_0
    in mkUnsafeQuaternion v
  where
    sgn ref x = if ref >= 0 then x else -x
    (Mat3 (Vec3 g11 g12 g13) (Vec3 g21 g22 g23) (Vec3 g31 g32 g33)) = m
    q0 = 0.5 * (sqrt $ g11 + g22 + g33 + 1)
    q1_0 = sqrt $ (g11 + 1) / 2
    q2_0 = sqrt $ (g22 + 1) / 2
    q3_0 = sqrt $ (g33 + 1) / 2

-- ================================== Other functions ==================================

-- | Find the antipodal (-q) quaternion that represents the same rotation.
antipodal :: Quaternion -> Quaternion
antipodal = mkUnsafeQuaternion . neg . quaterVec

-- | Apply a /active/ rotation on vector by quaternion. The vector is rotated around the
-- quaternion axis.
activeVecRotation :: Vec3 -> Quaternion -> Vec3
activeVecRotation v q = let
  (q0, qv) = splitQuaternion q
  qvl = normsqr qv
  k1 = (q0*q0 - qvl) *& v
  k2 = (2 * (v &. qv)) *& qv
  k3 = (2 * q0) *& (qv &^ v)
  in k1 &+ k2 &+ k3

-- | Apply a /passive/ rotation (change of basis) on vector by quaternion. The passive
-- rotation returns the equivalent input vector in the reference frame (basis) rotated by
-- the quaternion.
passiveVecRotation :: Vec3 -> Quaternion -> Vec3
passiveVecRotation v q = let
  (q0, qv) = splitQuaternion q
  k1 = (2*q0*q0 - 1) *& v
  k2 = (2 * (v &. qv)) *& qv
  k3 = (2 * q0) *& (v &^ qv)
  in k1 &+ k2 &+ k3

-- | This is shortest path interpolation in the space of rotations; however
-- this is achieved by possibly flipping the first endpoint in the space of
-- quaternions. Thus @slerpU 0.001 q1 q2@ may be very far from @q1@ (and very
-- close to @negU q1@) in the space of quaternions (but they are very close
-- in the space of rotations).
quaterInterpolation :: Double -> Quaternion -> Quaternion -> Quaternion
quaterInterpolation t (Quaternion pa) (Quaternion pb) = Quaternion v
  where
    v = (p &* y) &+ (pb &* yb)

    dab = pa &. pb
    (d, p) = if dab >= 0
             then ( dab,     pa)
             else (-dab, neg pa)

    omega = acosSafe d
    s = sin omega
    y  = sin (omega * (1-t)) / s
    yb = sin (omega *    t ) / s

-- | Calculates the approximated integer representation of poles (e.g. < 10 2 2 >).
aproxToIdealAxis :: Vec3 -> Double -> (Int, Int, Int)
aproxToIdealAxis (Vec3 x y z) err = let
  ls    = [x, y, z]
  ns    = map (numerator   . ratio) ls
  ds    = map (denominator . ratio) ls
  mmc   = foldr lcm 1 ds
  mdc   = foldr gcd 0 nzIntVaules
  ratio = (flip approxRational) err
  intValues    = zipWith (\n d -> fromIntegral $ n * (div mmc d)) ns ds
  nzIntVaules  = filter (/= 0) intValues
  [x', y', z'] = case nzIntVaules of
    [] -> intValues
    _  -> map (`div` mdc) intValues
  in (x', y', z')

-- ================================ Show instances ===============================

instance Show AxisPair where
  show (AxisPair (vec, theta)) = let
    --(x, y, z) = aproxToIdealAxis vec 0.05
    in showPretty vec ++ " @ " ++ (show $ (toAngle theta :: Deg))

instance Show Quaternion where
  show (Quaternion v) = "Q" ++ showPretty v

instance Show Rodrigues where
  show (Rodrigues v) = "FR" ++ showPretty v

instance Show Euler where
  show Euler{..} = let
    foo a = show $ (toAngle a :: Deg)
    in "Euler(" ++ foo phi1 ++ ", " ++ foo phi ++ ", " ++ foo phi2 ++ ")"

instance Show RotMatrix where
  show (RotMatrix m) = showPretty m

instance Show Deg where
  show (Deg theta) = showFFloat (Just 1) theta " deg"

instance Show Rad where
  show (Rad theta) = showFFloat (Just 1) theta " rad"


-- ================================== Test ======================================

-- | Verify the correctness of this module
testOrientationModule :: IO ()
testOrientationModule = let
  an = Deg 90
  test =
    [ toQuaternion $ mkAxisPair (Vec3 1 0 0) (Rad 0)
    , toQuaternion $ mkAxisPair (Vec3 1 0 0) (Rad $ 0.5*pi)
    , toQuaternion $ mkAxisPair (Vec3 0 1 0) (Rad $ 0.5*pi)
    , toQuaternion $ mkAxisPair (Vec3 0 0 1) (Rad $ 0.5*pi)
    , toQuaternion $ mkAxisPair (Vec3 1 1 1) (Rad $ 1.5*pi)
    , toQuaternion $ mkAxisPair (Vec3 1 0 0) (Rad $ pi)
    , toQuaternion $ mkEuler (Deg 0.0)  (Deg 54.7) (Deg 45.0)
    , toQuaternion $ mkEuler (Deg 30.0) (Deg 54.7) (Deg 45.0)
    ]
  in do
    putStrLn  "-- Checking RotationSystem module..."
    testAngleConv an
    mapM_ showAll test
    mapM_ testConv test
    testComposition
    testMisso
    putStrLn $ "-- Checking end --"

showAll :: Quaternion -> IO ()
showAll q = do
  putStrLn "------------------------------------------------"
  putStrLn $ show (fromQuaternion q :: Quaternion)
  putStrLn $ show (fromQuaternion q :: Rodrigues)
  putStrLn $ show (fromQuaternion q :: AxisPair)
  putStrLn $ show (fromQuaternion q :: Euler)
  putStrLn $ show (fromQuaternion q :: RotMatrix)
  putStrLn "\n"

testAngleConv :: Deg -> IO ()
testAngleConv an = do
  let
    an1 = (toAngle $ fromAngle an) :: Deg
  putStrLn "-- Angle convertion"
  print $ an == an1

testConv :: Quaternion -> IO ()
testConv q = do
  let
    ap  = (fromQuaternion q)::AxisPair
    qap = toQuaternion ap
    e   = (fromQuaternion q)::Euler
    qe  = toQuaternion e
    r   = (fromQuaternion q) :: Rodrigues
    qr  = toQuaternion r
    m   = (fromQuaternion q) :: RotMatrix
    qm  = toQuaternion m
    outm  = show $ getOmega (q -@- qm)
    outap = show $ getOmega (q -@- qap)
    oute  = show $ getOmega (q -@- qe)
    outr  = show $ getOmega (q -@- qr)
  putStrLn "------------------------------------------------"
  putStrLn $ "-- Quaternion: "           ++ show q
  putStrLn $ "-- Euler conversion: "     ++ oute
  putStrLn $ "-- AxisPair conversion: "  ++ outap
  putStrLn $ "-- Rodrigues conversion: " ++ outr
  putStrLn $ "-- RotMatrix conversion: " ++ outm

testComposition :: IO ()
testComposition = let
  e1 = mkEuler (Deg 30) (Deg 0)  (Deg 0)
  e2 = mkEuler (Deg 0)  (Deg 30) (Deg 0)
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

  in do
    putStrLn $ "\n-- Testing composition "
    putStrLn $ " Rotation e1 = " ++ show e1
    putStrLn $ " followed by e2 = " ++ show e2
    putStrLn $ " where e1 #<= e2 = Euler(30.0 deg, 30.0 deg, 0.0 deg)"
    putStrLn $ " Results:"
    putStrLn $ "Euler      => " ++ show ((fromQuaternion $ toQuaternion ce) :: Euler)
    putStrLn $ "Quaternion => " ++ show ((fromQuaternion $ toQuaternion cq) :: Euler)
    putStrLn $ "RotMatrix  => " ++ show ((fromQuaternion $ toQuaternion cm) :: Euler)
    putStrLn $ "AxisPair   => " ++ show ((fromQuaternion $ toQuaternion ca) :: Euler)
    putStrLn $ "Rodrigues  => " ++ show ((fromQuaternion $ toQuaternion cr) :: Euler)

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

  in do
    putStrLn $ "\n-- Testing missorientations."
    putStrLn $ " Rotation e1 = " ++ show e1
    putStrLn $ " Rotation e2 = " ++ show e2
    putStrLn $ " Results:"
    putStrLn $ " e1 -@- e2 => " ++ show ((fromQuaternion $ toQuaternion m)   :: AxisPair)
    putStrLn $ " e2 -@- e1 => " ++ show ((fromQuaternion $ toQuaternion mi)  :: AxisPair)
    putStrLn $ " e1 -#- e2 => " ++ show ((fromQuaternion $ toQuaternion cm)  :: AxisPair)
    putStrLn $ " e2 -#- e1 => " ++ show ((fromQuaternion $ toQuaternion cmi) :: AxisPair)
    putStrLn $ " e1 #<= (e2 -#- e1) == e2 == " ++ show ((fromQuaternion testm) :: Euler)
