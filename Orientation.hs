-----------------------------------------------------------------------------
--
-- Module      :  RotationSystem
-- Copyright   :
-- License     :  Private
--
-- Maintainer  :  Edgar Gomes de Araujo
-- Stability   :  No
-- Portability :
--
-- |  Module to define differnt rotation representations,
-- coverntions and operations between them
--
-----------------------------------------------------------------------------

module Hammer.Texture.Orientation
(Euler(..),
Rodrigues(..),
Quartenion(..),
AxisPair(..),
Angle(..),
RotMatrix(..),
toRodrigues,
fromRodrigues,
(*@*),
(-@-),
misM,
compM,
getTheta,
toRad,
toDeg,
checkRotationSystemModule
)where

import Data.Ratio
import Data.List (intercalate)
import Control.Monad
import Numeric

import Hammer.Math.Vector

import Debug.Trace
debug :: Show a => String -> a -> a
debug s x = trace (s ++ show x) x

data Angle = Deg Double
           | Rad Double
             deriving (Eq)

data Euler      = Euler (Angle, Angle, Angle)
                  deriving (Show, Eq)

data Rodrigues  = Rodrigues Vec3
                  deriving (Eq)

data Quartenion = Quartenion Vec4
                  deriving (Show)

data RotMatrix  = RotMatrix Mat3
                  deriving (Eq)

data AxisPair   = AxisPair Vec3 Angle
                  deriving (Eq)


--------------------------
-- Rot class and instances
class Rot a where
  (+@+)          ::a -> a -> a
  (-@-)          ::a -> a -> a
  compose        ::a -> a -> a
  misorientation ::a -> a -> a
  inverse        ::a -> a -> a
  
  toRodrigues    ::a -> Rodrigues
  fromRodrigues  ::Rodrigues -> a
  
  compose a b        = a +@+ b
  misorientation a b = a -@- b
  a -@- b            = a +@+ (inverse b)




---------------------------
-- Instances of Rot class
instance Rot Rodrigues where
  (Rodrigues rA) +@+ (Rodrigues rB) = Rodrigues (r1 >/< r2)
    where
      r1 = ((rA &+ rB) &- (rA * rB))
      r2 = (1.0 - (rA &. rB))

  inverse (Rodrigues rA) = Rodrigues ((-1) *& rA)
  toRodrigues   = id
  fromRodrigues = id  
  
  
instance Rot AxisPair where
  a +@+ b = toRodrigues a +@+ toRodrigues b
  inverse (AxisPair vec a) = AxisPair (neg v) a 
  toRodrigues (AxisPair v alpha) = Rodrigues r
    where   len = sqrt (v >.< v)
            (Rad theta) = toRad alpha
            k = tan (theta/2)
            r   | len == 0 = (0,0,0)
                | otherwise = (k/len) >*< v

  fromRodrigues (Rodrigues r) = AxisPair v (Rad theta)
    where   len = sqrt (r >.< r)
            theta = 2*atan (len)
            v   | len == 0 = (0,0,0)
                | otherwise = r >/< len


instance Rot Euler where
    toRodrigues (Euler (phi1, phi, phi2)) = Rodrigues (r1, r2, r3)
        where
            (Rad phi1', Rad phi', Rad phi2') = (toRad phi1, toRad phi, toRad phi2)
            -- P. Neumann (1991). “Representation of orientations of symmetrical objects by Rodrigues vectors.”
            -- Textures and Microstructures 14-18: 53-58
            diffphi = (phi2' - phi1')/2
            addphi = (phi2' + phi1')/2
            tanPHI = tan (phi'/2)
            cosaddphi = cos addphi

            r1 = tanPHI*((cos diffphi) / cosaddphi)
            r2 = -1*tanPHI*((sin diffphi) / cosaddphi)
            r3 = tan addphi

    fromRodrigues r = undefined

instance Rot Quartenion where
    toRodrigues e = undefined
    fromRodrigues r = undefined

instance Rot RotMatrix where
  a +@+ b = a *&* b
  inverse = transpose
  toRodrigues (RotMatrix m) = toRodrigues axispair
    where   ((a11, a12, a13),(a21, a22, a23),(a31, a32, a33)) = m
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
            axispair    | q == 0 = AxisPair (c1',c2',c3') (Rad theta)
                        | otherwise = AxisPair (c1,c2,c3) (Rad theta)

  fromRodrigues r = RotMatrix ((g11, g12, g13),(g21, g22, g23),(g31, g32, g33))
    where   (AxisPair (r1,r2,r3) alpha) = fromRodrigues r
            Rad theta = toRad alpha
            -- from  Introduction to texture analysis : macrotexture, microtexture, and orientation
            -- mapping, Olaf Engler and Valerie Randle, 2nd ed.
            costheta = cos theta
            sintheta = sin theta
            g11 = (1 - r1^2) * costheta + r1^2
            g12 = r1 * r2 * (1 - costheta) + r3 * sintheta
            g13 = r1 * r3 * (1 - costheta) - r2 * sintheta
            g21 = r1 * r2 * (1 - costheta) - r3 * sintheta
            g22 = (1 - r2^2 ) * costheta + r2^2
            g23 = r2 * r3 * (1 - costheta) + r1 * sintheta
            g31 = r1 * r3 * (1 - costheta) + r2 * sintheta
            g32 = r2 * r3 * (1 - costheta) - r1 * sintheta
            g33 = (1 - r3^2 ) * costheta + r3^2


getTheta::Rodrigues -> Double
getTheta (Rodrigues rA) = 2 * atan nor
    where
    nor = sqrt( rA >.< rA)


toDeg deg@(Deg theta) = deg
toDeg (Rad theta) = Deg $ theta*(180/pi)


toRad rad@(Rad theta) = rad
toRad (Deg theta) = Rad $ theta*(pi/180)



-- Quartenion functions
----------------------

eulerToQuar::Euler -> Quartenion
eulerToQuar (Euler (phi1, phi, phi2)) = Quartenion (q1, q2, q3, q4)
    where
        (Rad phi1', Rad phi', Rad phi2') = (toRad phi1, toRad phi, toRad phi2)
        dphi = (phi2' - phi1')/2
        aphi = (phi1' + phi2')/2
        s1 = sin (phi'/2)
        c1 = cos (phi'/2)
        q1 = (cos aphi) * c1
        q2 = (sin dphi) * s1
        q3 = (cos dphi) * s1
        q4 = (sin aphi) * c1


quarToEuler::Quartenion -> Euler
quarToEuler (Quartenion (q1, q2, q3, q4)) = Euler (Rad phi1, Rad phi, Rad phi2)
    where
        phi1 = atan(q4/q1)-atan(q2/q3)
        phi = 2 * (atan $ sqrt((q2*q2+q3*q3)/(q4*q4+q1*q1)))
        phi2 = atan(q4/q1)+atan(q2/q3)











-----------------
-- Show instances
instance Show AxisPair where
    show (AxisPair vec theta) = "< " ++ show x ++ ", " ++ show y ++ ", " ++ show z ++ " > " ++ (show $ toDeg theta)
        where (x, y, z) = aproxToIdealAxis vec 0.05

instance Show Rodrigues where
    show (Rodrigues (x,y,z)) = "FR(" ++ (showNumber x)++ ", " ++ (showNumber y) ++ ", " ++ (showNumber z) ++ ")"
        where showNumber x = showEFloat (Just 1) x " "

instance Show RotMatrix where
    show (RotMatrix m) =
        let ((a11, a12, a13),(a21, a22, a23),(a31, a32, a33)) = m
            showNumber x = showEFloat (Just 1) x " "
        in     "| " ++ showNumber a11 ++ "\t" ++ showNumber a12 ++ "\t" ++ showNumber a13 ++ " |\n"
            ++ "| " ++ showNumber a21 ++ "\t" ++ showNumber a22 ++ "\t" ++ showNumber a23 ++ " |\n"
            ++ "| " ++ showNumber a31 ++ "\t" ++ showNumber a32 ++ "\t" ++ showNumber a33 ++ " |\n"

instance Show Angle where
    show (Deg theta) = showFFloat (Just 1) theta " deg"
    show (Rad theta) = showFFloat (Just 1) theta " rad"

-- | Scale a list of num in such way that the min non zero absolute value is 1.0
aproxToIdealAxis::(Double, Double, Double) -> Double -> (Integer, Integer, Integer)
aproxToIdealAxis (x,y,z) err =  (x', y', z')
    where
        ls = [x,y,z]
        ratio num = approxRational num err
        ns = map (numerator.ratio) ls
        ds = map (denominator.ratio) ls
        mmc = foldr lcm 1 ds
        intValues = zipWith (\n d -> n*(div mmc d)) ns ds
        [x', y', z'] = case nzIntVaules of
            [] -> intValues
            _ -> map (\x -> (div x mdc)) intValues
            where
                nzIntVaules = filter (/= 0) intValues
                mdc = foldr gcd 0 nzIntVaules















-- | Verify the correctness of this module
checkRotationSystemModule::IO ()
checkRotationSystemModule = do
    let an = 90
        a = toRodrigues $ AxisPair (1,1,1) (Rad $ pi/3)
        b = toRodrigues $ AxisPair (0,0,1) (Rad $ pi/2)
        c = toRodrigues $ AxisPair (1,0,0) (Rad $ pi)
        d = toRodrigues $ AxisPair (0,-1,-1) (Rad $ pi/2)
        e = toRodrigues $ AxisPair (-1,0,0) (Rad $ pi/4)
        f = toRodrigues $ AxisPair (1,0,0) (Rad 0)
        g = toRodrigues $ Euler (Deg 0.0, Deg 54.7, Deg 45.0)
        h = toRodrigues $ Euler (Deg 30.0, Deg 54.7, Deg 45.0)
        test = [a,b,c,d,e,f,g,h]

    putStrLn  "-- Checking RotationSystem module..."
    testAngleConv an
    mapM testAxisPairConv test
    --mapM testEulerConv test
    mapM (\x -> mapM (foldRot x) test) test
    mapM (\x -> mapM (foldMis x) test) test
    putStrLn $ show $ (Rodrigues (-1,0,0)) *@* (Rodrigues (0,0.5,0))
    putStrLn $ show $ toRodrigues $ compM (fromRodrigues $ Rodrigues (0,0.5,0)) (fromRodrigues $ Rodrigues (-1,0,0))
    putStrLn $ "--- check end"

    where

        testAngleConv an = do
            let a' = toRad (Deg an)
                a'' = toDeg $ toRad (Deg an)
            putStrLn "\n-- Angle convertion"
            putStrLn $ ">> " ++ show an ++ " -> " ++ show a' ++ " -> " ++ show a''

        testAxisPairConv a = do
            let a' = (fromRodrigues a)::AxisPair
                a'' = toRodrigues a'
            putStrLn "\n-- AxisPair convertion"
            putStrLn $ ">> " ++ show a ++ " -> " ++ show a' ++ " -> " ++ show a''

        testEulerConv a = do
            let a' = (fromRodrigues a)::Euler
                a'' = toRodrigues a'
            putStrLn "\n-- Euler convertion"
            putStrLn $ ">> " ++ show a ++ " -> " ++ show a' ++ " -> " ++ show a''

        foldRot::Rodrigues -> Rodrigues -> IO Rodrigues
        foldRot a b = do
            putStrLn "\n-- Compose rotation"
            putStrLn $ ">> a = " ++ show ((fromRodrigues a)::AxisPair)
            putStrLn $ ">> b = " ++ show ((fromRodrigues b)::AxisPair)
            let c = a *@* b
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




