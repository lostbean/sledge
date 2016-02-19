{-# LANGUAGE
    FlexibleContexts
  , GeneralizedNewtypeDeriving
  , MultiParamTypeClasses
  , NamedFieldPuns
  , RecordWildCards
  , TypeFamilies
  #-}

-- |
-- Module      : Texture.Orientation
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
module Texture.Orientation
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
       , getShortOmega
       , getAbsShortOmega
         -- * Quaternions
       , Quaternion (quaterVec)
       , mkUnsafeQuaternion
       , mkQuaternion
       , splitQuaternion
       , unsafeMergeQuaternion
       , mergeQuaternion
       , antipodal
       , composeQ0
         -- * Euler
       , Euler      (phi1, phi, phi2)
       , mkEuler
         -- * Axis-Angle
       , AxisPair   (axisAngle)
       , mkAxisPair
         -- * Frank-Rodrigues
       , Rodrigues  (rodriVec)
       , mkRodrigues
       , mkUnsafeRodrigues
         -- * Matrix
       , RotMatrix  (rotMat)
       , mkUnsafeRotMatrix
         -- * Vector rotation
       , activeVecRotation
       , passiveVecRotation
         -- * Other functions
       , getRTNdir
       , get100dir
       , quaterInterpolation
       , averageQuaternion
       , getScatterMatrix
       , aproxToIdealAxis
       , getAbsShortAngle
       , getShortAngle
       ) where

import qualified Data.List                   as L
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Generic.Base    as G
import qualified Data.Vector.Generic.Mutable as M

import           Control.Monad               (liftM)
import           Data.Vector                 (Vector)

import           Control.DeepSeq
import           Foreign
import           Data.Ratio
import           Numeric

import           Hammer.Math.Algebra
import           System.Random

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
-- Uses all directions with rotations from 0 to PI.
newtype Quaternion =
  Quaternion
  { quaterVec :: Vec4
  } deriving (Eq, NFData)

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
  } deriving (Eq, NFData)

-- | Frank-Rodrigues representation.
newtype Rodrigues =
  Rodrigues
  { rodriVec :: Vec3
  } deriving (Eq, NFData)

-- | Frank-Rodrigues representation.
newtype RotMatrix =
  RotMatrix
  { rotMat :: Mat3
  } deriving (MultSemiGroup, Transposable, Inversable, IdMatrix, Eq, NFData)

instance Matrix RotMatrix

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
mkQuaternion v
  | l == 0    = Quaternion (Vec4 0 0 0 1)
  | q0 > 0    = Quaternion (v &* (1/l))
  | otherwise = Quaternion (v &* (-1/l))
  where
    l  = len v
    q0 = _1 v

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
mkUnsafeQuaternion v
  | (_1 v) > 0 = Quaternion v
  | otherwise  = Quaternion $ neg v

-- | /Unsafe/ rotation matrix constructor. It assumes a orthonormal matrix with
-- unitary vector as input.
mkUnsafeRotMatrix :: Mat3 -> RotMatrix
mkUnsafeRotMatrix = RotMatrix

-- | /Unsafe/ Frank-Rodrigues constructor.
mkUnsafeRodrigues :: Vec3 -> Rodrigues
mkUnsafeRodrigues = Rodrigues

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

-- | Get rotation angle in radians in [ -pi, pi ] range.
getShortOmega :: (Rot a)=> a -> Double
getShortOmega = getShortAngle . getOmega

-- | Get the absolute rotation angle in radians in [ 0, pi ] range. This notation doesn't
-- specify if the rotation is clockwise or anti-clockwise
getAbsShortOmega :: (Rot a)=> a -> Double
getAbsShortOmega = getAbsShortAngle . getOmega


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
    in mkUnsafeQuaternion (Vec4 q0 q1 q2 q3)
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

-- | Get the RD, TD and ND directions in the crystal frame.
getRTNdir :: RotMatrix -> (Vec3, Vec3, Vec3)
getRTNdir (RotMatrix m) = let
  ( Mat3 rd td nd ) = transpose m
  in (rd, td, nd)

-- | Get the axis [100] (crystal base) in the sample frame.
get100dir :: RotMatrix -> (Vec3, Vec3, Vec3)
get100dir (RotMatrix (Mat3 e1 e2 e3)) = (e1, e2, e3)

-- | Converts an angle in Radians from the @[0, 2*pi]@ range to @[-pi, pi]@ range.
getAbsShortAngle :: Double -> Double
getAbsShortAngle a
  | a > pi    = 2 * pi - a
  | otherwise = a

-- | Converts an angle in Radians from the @[0, 2*pi]@ range to @[-pi, pi]@ range.
getShortAngle :: Double -> Double
getShortAngle a
  | a > pi    = a - 2 * pi
  | otherwise = a

-- | Fast calculation of @q0@ component of a quaternion composition. It is the same as
-- '#<=' or '=>#' for @q0@ which gives information about the rotation angle. The input
-- order doesn't matter in this case. This function is used for fast calculation of the
-- minimum misorientation angle or the misorientation angle itself.
composeQ0 :: Quaternion -> Quaternion -> Double
composeQ0 p q = let
  (p0, pv) = splitQuaternion p
  (q0, qv) = splitQuaternion q
  in p0 * q0 - pv &. qv
{-# INLINE composeQ0 #-}

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
quaterInterpolation t (Quaternion pa) (Quaternion pb) = mkUnsafeQuaternion v
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

-- | Calculates the scatter matrix that describes a given distribution.
getScatterMatrix :: U.Vector Quaternion -> Mat4
getScatterMatrix qs
  | n > 0     = total &* (1/n)
  | otherwise = zero
  where
    n     = fromIntegral (U.length qs)
    total = U.foldl' func zero qs
    func acc q = let
      v = quaterVec q
      in acc &+ (outer v v)

-- | Calculates the average orientation from a distribution.
-- See "Averaging Quaternions", FL Markley, 2007
averageQuaternion :: U.Vector Quaternion -> Quaternion
averageQuaternion vq = mkQuaternion $ case i of
  0 -> _1 v
  1 -> _2 v
  2 -> _3 v
  _ -> _4 v
  where
    scatter  = getScatterMatrix vq
    (ev, ex) = symmEigen scatter
    -- sort decomposition by eigenvalues
    -- (highest eigenvalue = highest concentration value = distribution mode)
    v = transpose ev
    i = U.maxIndex $ U.fromList [_1 ex, _2 ex, _3 ex, _4 ex]

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

-- -------------------------------------------- Unbox Quaternion ----------------------------------------------------

newtype instance U.MVector s Quaternion = MV_Quaternion (U.MVector s Vec4)
newtype instance U.Vector    Quaternion = V_Quaternion  (U.Vector    Vec4)

instance U.Unbox Quaternion

instance M.MVector U.MVector Quaternion where
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicOverlaps #-}
  {-# INLINE basicUnsafeNew #-}
  {-# INLINE basicUnsafeReplicate #-}
  {-# INLINE basicUnsafeRead #-}
  {-# INLINE basicUnsafeWrite #-}
  {-# INLINE basicClear #-}
  {-# INLINE basicSet #-}
  {-# INLINE basicUnsafeCopy #-}
  {-# INLINE basicUnsafeGrow #-}
  basicLength (MV_Quaternion v)                         = M.basicLength v
  basicUnsafeSlice i n (MV_Quaternion v)                = MV_Quaternion $ M.basicUnsafeSlice i n v
  basicOverlaps (MV_Quaternion v1) (MV_Quaternion v2)   = M.basicOverlaps v1 v2
  basicInitialize (MV_Quaternion v)                     = M.basicInitialize v
  basicUnsafeNew n                                      = MV_Quaternion `liftM` M.basicUnsafeNew n
  basicUnsafeReplicate n (Quaternion x)                 = MV_Quaternion `liftM` M.basicUnsafeReplicate n x
  basicUnsafeRead (MV_Quaternion v) i                   = M.basicUnsafeRead v i >>= (return . Quaternion)
  basicUnsafeWrite (MV_Quaternion v) i (Quaternion x)   = M.basicUnsafeWrite v i x
  basicClear (MV_Quaternion v)                          = M.basicClear v
  basicSet (MV_Quaternion v) (Quaternion x)             = M.basicSet v x
  basicUnsafeCopy (MV_Quaternion v1) (MV_Quaternion v2) = M.basicUnsafeCopy v1 v2
  basicUnsafeGrow (MV_Quaternion v) n                   = MV_Quaternion `liftM` M.basicUnsafeGrow v n

instance G.Vector U.Vector Quaternion where
  {-# INLINE basicUnsafeFreeze #-}
  {-# INLINE basicUnsafeThaw #-}
  {-# INLINE basicLength #-}
  {-# INLINE basicUnsafeSlice #-}
  {-# INLINE basicUnsafeIndexM #-}
  {-# INLINE elemseq #-}
  basicUnsafeFreeze (MV_Quaternion v)                 = V_Quaternion `liftM` G.basicUnsafeFreeze v
  basicUnsafeThaw (V_Quaternion v)                    = MV_Quaternion `liftM` G.basicUnsafeThaw v
  basicLength (V_Quaternion v)                        = G.basicLength v
  basicUnsafeSlice i n (V_Quaternion v)               = V_Quaternion $ G.basicUnsafeSlice i n v
  basicUnsafeIndexM (V_Quaternion v) i                = G.basicUnsafeIndexM v i >>= (return . Quaternion)
  basicUnsafeCopy (MV_Quaternion mv) (V_Quaternion v) = G.basicUnsafeCopy mv v
  elemseq _ (Quaternion x) t                          = G.elemseq (undefined :: Vector a) x t
