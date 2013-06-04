{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Hammer.Texture.Orientation
       ( Euler
       , AxisPair
       , Rodrigues  (..)
       , Quaternion (..)
       , Angle      (..)
       , RotMatrix  (..)
       , mkAxisPair
       , mkEuler
       , toRodrigues
       , fromRodrigues
       , (+@+)
       , (-@-)
       , invert
       , getTheta
       , toRad
       , toDeg
       , checkRotationSystemModule
       ) where

import Data.Ratio
import Numeric

import Hammer.Math.Algebra

--import Debug.Trace
--dbg s x = trace (s ++ show x) x

data Angle = Deg { unDeg :: Double }
           | Rad { unRad :: Double }
             deriving (Eq)

data Euler =
  Euler
  { phi1 :: Double
  , phi  :: Double
  , phi2 :: Double
  } deriving (Show, Eq)
             
newtype Rodrigues  = Rodrigues  Vec3           deriving (Eq)
newtype RotMatrix  = RotMatrix  Mat3           deriving (MultSemiGroup, Matrix, Eq)
newtype Quaternion = Quaternion (Double, Vec3) deriving (Eq, Show)
newtype AxisPair   = AxisPair   (Vec3, Double) deriving (Eq)

mkAxisPair :: Vec3 -> Angle -> AxisPair
mkAxisPair r omega = let
  Rad theta = toRad omega
  in AxisPair (normalize r, theta)

mkEuler :: Angle -> Angle -> Angle -> Euler
mkEuler p1 p p2 = Euler
  { phi1 = unRad $ toRad p1
  , phi  = unRad $ toRad p
  , phi2 = unRad $ toRad p2 }

-- | This class defines rotation composition
class Rot a where
  (+@+)          :: a -> a -> a
  (-@-)          :: a -> a -> a
  compose        :: a -> a -> a
  misorientation :: a -> a -> a
  invert         :: a -> a
  
  toRodrigues    :: a -> Rodrigues
  fromRodrigues  :: Rodrigues -> a
  
  compose a b        = a +@+ b
  misorientation a b = a -@- b
  a -@- b            = a +@+ (invert b)

-- =========================== Instances of Rot class ==========================

instance Rot Rodrigues where
  (Rodrigues rA) +@+ (Rodrigues rB) = let
    r = (rA &+ rB) &- (rA &^ rB)
    k = 1.0 - (rA &. rB)
    in Rodrigues $ r &* (1 / k)
  invert (Rodrigues r) = Rodrigues $ neg r
  toRodrigues   = id
  fromRodrigues = id  
    
instance Rot AxisPair where
  a +@+ b = fromRodrigues $ (toRodrigues a) +@+ (toRodrigues b)
  invert (AxisPair (v, a)) = AxisPair (neg v, a)
  toRodrigues (AxisPair (v, theta))
    | l == 0    = Rodrigues zero
    | otherwise = Rodrigues $ (k / l) *& v
    where
      l = norm v
      k = tan (theta / 2)
  fromRodrigues (Rodrigues r) = AxisPair (v, theta)
    where
      v | l == 0    = zero
        | otherwise = r &* (l / 1)
      l     = norm r
      theta = 2 * atan l

instance Rot Euler where
  toRodrigues Euler{..} = let
    -- P. Neumann (1991). “Representation of orientations of symmetrical objects by Rodrigues vectors.”
    -- Textures and Microstructures 14-18: 53-58
    diffphi = (phi2 - phi1)/2
    addphi  = (phi2 + phi1)/2
    tanPHI  = tan (phi/2)
    cosaddphi = cos addphi

    r1 = tanPHI * ((cos diffphi) / cosaddphi)
    r2 = (-1) * tanPHI * ((sin diffphi) / cosaddphi)
    r3 = tan addphi
    in Rodrigues $ Vec3 r1 r2 r3
  fromRodrigues = undefined

instance Rot Quaternion where
  (Quaternion (p0, p)) +@+ (Quaternion (q0, q)) = let
    pq0 = p0 * q0 - p &. q
    pq  = p0 *& q &+ q0 *& p &+ p &^ q
    in Quaternion (pq0, pq)
  invert (Quaternion (q0, q)) = Quaternion (q0, neg q)
  toRodrigues (Quaternion (q0, q)) = let
    in undefined 
  fromRodrigues x = let
    q0 = cos (theta / 2)
    q  = r &* (sin $ theta / 2)
    AxisPair (r, theta) = fromRodrigues x
    in Quaternion (q0, q)

instance Rot RotMatrix where
  a +@+ b = a .*. b
  invert  = transpose
  toRodrigues (RotMatrix m) = let
    (Mat3 (Vec3 a11 a12 a13) (Vec3 a21 a22 a23) (Vec3 a31 a32 a33)) = m
    -- from  Introduction to texture analysis : macrotexture, microtexture, and orientation
    -- mapping, Olaf Engler and Valerie Randle, 2nd ed.
    -- Obs
    -- the conv. is not always straight forward, e.g. "toRodrigues $ RotMatrix ( (-1.0,0.0,0.0),(0,0,-1),(0,-1,0))"
    -- result in (0,0,0) which is, by symmetry, is equivalent to the matrix!
    trace = a11 + a22 + a33
    theta = acos (0.5*(trace-1))
    q = 2 * (sin theta)
    c1 = (a23-a32)/q
    c2 = (a31-a13)/q
    c3 = (a12-a21)/q
    c1'= sqrt ((a11 + 1)/2)
    c2'= sqrt ((a22 + 1)/2)
    c3'= sqrt ((a33 + 1)/2)
    axispair
      | q == 0    = AxisPair (Vec3 c1' c2' c3', theta)
      | otherwise = AxisPair (Vec3 c1 c2 c3   , theta)
    in toRodrigues axispair
       
  fromRodrigues r = let
    AxisPair (Vec3 r1 r2 r3, theta) = fromRodrigues r
    -- from  Introduction to texture analysis : macrotexture, microtexture, and orientation
    -- mapping, Olaf Engler and Valerie Randle, 2nd ed.
    costheta = cos theta
    sintheta = sin theta
    g11 = (1 - r1*r1) * costheta + r1*r1
    g12 = r1 * r2 * (1 - costheta) + r3 * sintheta
    g13 = r1 * r3 * (1 - costheta) - r2 * sintheta
    g21 = r1 * r2 * (1 - costheta) - r3 * sintheta
    g22 = (1 - r2*r2) * costheta + r2*r2
    g23 = r2 * r3 * (1 - costheta) + r1 * sintheta
    g31 = r1 * r3 * (1 - costheta) + r2 * sintheta
    g32 = r2 * r3 * (1 - costheta) - r1 * sintheta
    g33 = (1 - r3*r3) * costheta + r3*r3
    in RotMatrix (Mat3 (Vec3 g11 g12 g13) (Vec3 g21 g22 g23) (Vec3 g31 g32 g33))

getTheta :: Rodrigues -> Double
getTheta (Rodrigues r) = 2 * atan (norm r)

toDeg :: Angle -> Angle
toDeg deg@(Deg _) = deg
toDeg (Rad theta) = Deg $ theta * (180 / pi)

toRad :: Angle -> Angle
toRad rad@(Rad _) = rad
toRad (Deg theta) = Rad $ theta * (pi / 180)

-- | Scale a list of num in such way that the min non zero absolute value is 1.0
aproxToIdealAxis :: Vec3 -> Double -> (Integer, Integer, Integer)
aproxToIdealAxis (Vec3 x y z) err = let
  ls = [x, y, z]
  ratio num = approxRational num err
  ns = map (numerator . ratio) ls
  ds = map (denominator . ratio) ls
  mmc = foldr lcm 1 ds
  intValues = zipWith (\n d -> n*(div mmc d)) ns ds
  [x', y', z'] = case nzIntVaules of
    [] -> intValues
    _ -> map (\k -> (div k mdc)) intValues
    where
      nzIntVaules = filter (/= 0) intValues
      mdc = foldr gcd 0 nzIntVaules
  in (x', y', z')

-- =========================== Quaternion functions ============================

eulerToQuar :: Euler -> Quaternion
eulerToQuar Euler{..} = Quaternion (q0, Vec3 q1 q2 q3)
  where
    dphi = (phi2 - phi1) / 2
    aphi = (phi1 + phi2) / 2
    s1 = sin (phi / 2)
    c1 = cos (phi / 2)
    q0 = (cos aphi) * c1
    q1 = (sin dphi) * s1
    q2 = (cos dphi) * s1
    q3 = (sin aphi) * c1

quarToEuler :: Quaternion -> Euler
quarToEuler (Quaternion (q0, Vec3 q1 q2 q3)) = Euler
  { phi1 = atan (q3 / q0) - atan (q1 / q2)
  , phi  = 2 * (atan $ sqrt((q1 * q1 + q2 * q2) / (q3 * q3 + q0 * q0)))
  , phi2 = atan (q3 / q0) + atan (q1 / q2) }

-- ================================ Show instances ===============================

instance Show AxisPair where
  show (AxisPair (vec, theta)) = let
    (x, y, z) = aproxToIdealAxis vec 0.05
    in "< " ++ show x ++ ", " ++ show y ++
       ", " ++ show z ++ " > " ++ (show $ toDeg $ Rad theta)

instance Show Rodrigues where
  show (Rodrigues (Vec3 x y z)) = let
    foo n = showEFloat (Just 1) n " "
    in "FR(" ++ (foo x) ++ ", " ++ (foo y) ++
       ", " ++ (foo z) ++ ")"

instance Show RotMatrix where
  show (RotMatrix m) = let
    (Mat3 (Vec3 a11 a12 a13) (Vec3 a21 a22 a23) (Vec3 a31 a32 a33)) = m
    foo x = showEFloat (Just 1) x " "
    in "| " ++ foo a11 ++ "\t" ++ foo a12 ++ "\t" ++ foo a13 ++ " |\n" ++
       "| " ++ foo a21 ++ "\t" ++ foo a22 ++ "\t" ++ foo a23 ++ " |\n" ++
       "| " ++ foo a31 ++ "\t" ++ foo a32 ++ "\t" ++ foo a33 ++ " |\n"

instance Show Angle where
  show (Deg theta) = showFFloat (Just 1) theta " deg"
  show (Rad theta) = showFFloat (Just 1) theta " rad"


-- | Verify the correctness of this module
checkRotationSystemModule :: IO ()
checkRotationSystemModule = let
  an = 90
  a = toRodrigues $ mkAxisPair (Vec3 1 1 1) (Rad $ pi/3)
  b = toRodrigues $ mkAxisPair (Vec3 0 0 1) (Rad $ pi/2)
  c = toRodrigues $ mkAxisPair (Vec3 1 0 0) (Rad $ pi)
  d = toRodrigues $ mkAxisPair (Vec3 0 (-1) (-1)) (Rad $ pi/2)
  e = toRodrigues $ mkAxisPair (Vec3 (-1) 0 0) (Rad $ pi/4)
  f = toRodrigues $ mkAxisPair (Vec3 1 0 0) (Rad 0)
  g = toRodrigues $ mkEuler (Deg 0.0)  (Deg 54.7) (Deg 45.0)
  h = toRodrigues $ mkEuler (Deg 30.0) (Deg 54.7) (Deg 45.0)
  test = [a,b,c,d,e,f,g,h]
  in do
    putStrLn  "-- Checking RotationSystem module..."
    testAngleConv an
    mapM testAxisPairConv test
    --mapM testEulerConv test
    mapM (\x -> mapM (foldRot x) test) test
    mapM (\x -> mapM (foldMis x) test) test
    putStrLn $ show $ (Rodrigues $ Vec3 (-1) 0 0) +@+ (Rodrigues $ Vec3 0 0.5 0)
    putStrLn $ show $ toRodrigues $ (fromRodrigues $ Rodrigues (Vec3 0 0.5 0) :: RotMatrix) +@+ (fromRodrigues $ Rodrigues (Vec3 (-1) 0 0))
    putStrLn $ "--- check end"

    where

        testAngleConv an = do
            let a'  = toRad (Deg an)
                a'' = toDeg $ toRad (Deg an)
            putStrLn "\n-- Angle convertion"
            putStrLn $ ">> " ++ show an ++ " -> " ++ show a' ++ " -> " ++ show a''

        testAxisPairConv a = do
            let a'  = (fromRodrigues a)::AxisPair
                a'' = toRodrigues a'
            putStrLn "\n-- AxisPair convertion"
            putStrLn $ ">> " ++ show a ++ " -> " ++ show a' ++ " -> " ++ show a''

        testEulerConv a = do
            let a'  = (fromRodrigues a)::Euler
                a'' = toRodrigues a'
            putStrLn "\n-- Euler convertion"
            putStrLn $ ">> " ++ show a ++ " -> " ++ show a' ++ " -> " ++ show a''

        foldRot::Rodrigues -> Rodrigues -> IO Rodrigues
        foldRot a b = do
            putStrLn "\n-- Compose rotation"
            putStrLn $ ">> a = " ++ show ((fromRodrigues a)::AxisPair)
            putStrLn $ ">> b = " ++ show ((fromRodrigues b)::AxisPair)
            let c = a +@+ b
            putStrLn $ ">> a.b (compose) = " ++ show ((fromRodrigues c)::AxisPair)
            return c

        foldMis::Rodrigues -> Rodrigues -> IO Rodrigues
        foldMis a b = do
            putStrLn "\n-- Compose misorientation"
            putStrLn $ ">> a = " ++ show ((fromRodrigues a)::AxisPair)
            putStrLn $ ">> b = " ++ show ((fromRodrigues b)::AxisPair)
            let c = a -@- b
            putStrLn $ ">> a[-1].b (compose) = " ++ show ((fromRodrigues c)::AxisPair)
            return c



