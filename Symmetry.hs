-----------------------------------------------------------------------------
--
-- Module      :  SymmetricSystem
-- Copyright   :
-- License     :  Private
--
-- Maintainer  :  Edgar Gomes de Araujo
-- Stability   :  No
-- Portability :
--
-- |  Module to define symmetry for crystal orientations
--
-----------------------------------------------------------------------------

module Hammer.Texture.Symmetry
(Symm(..), toFZ, checkSymmetricSystemModule,mSymm) where

import Data.List
import Hammer.Math.Vector
import Hammer.Texture.Orientation


import Debug.Trace
debug :: Show a => String -> a -> a
debug s x = trace (s ++ show x) x


data Symm = Cubic | Hexagonal | Custom [SymmAxis]

data SymmAxis = Axis (Double, Double, Double) Int
data FZPlane = Plane {planeNormal::(Double, Double, Double), planeDist::Double}
data SymmOp = SymmOp (Double, Double, Double)

getVectorList::Symm -> [SymmAxis]
getVectorList Cubic =   [--Axis (0,0,1) 2, Axis (0,1,0) 2, Axis (1,0,0) 2,
                        --Axis (1,1,0) 2, Axis (1,0,1) 2, Axis (0,1,1) 2,
                        Axis (0,0,1) 4, Axis (0,1,0) 4, Axis (1,0,0) 4,
                        Axis (1,1,1) 3, Axis (1,1,-1) 3, Axis (1,-1,1) 3, Axis (1,-1,-1) 3]
getVectorList Hexagonal =   [Axis (0,0,1) 6, Axis (1,1,1) 3, Axis (1,-1,1) 3, Axis (1,1,-1) 3, Axis (1,-1,-1) 3]
getVectorList (Custom xs) = xs
-- getVectorList _ =   [Axis (1,1,1) 1]


getFZPlane::SymmAxis -> FZPlane
getFZPlane (Axis v n) = Plane vecnor dist
    where
        vecnor = normalVector v
        dist = abs $  tan (pi/(2*(fromIntegral n)))

getSymmOp::SymmAxis -> SymmOp
getSymmOp (Axis v n) = SymmOp c
    where
        lengthV = lengthVector v
        dist = tan (pi/(fromIntegral n))
        c = (dist/lengthV) >*< v

rotate::Double -> SymmOp -> Rodrigues -> Rodrigues
rotate dir (SymmOp cVec) r = s *@* r  -- add (debug "out" $) to check
    where s = Rodrigues $ (-1*dir) >*< cVec

isInFZ::Rodrigues -> FZPlane -> Double
isInFZ (Rodrigues r) (Plane nor dist)
    | isFZ = 0
    | rproj > 0 = 1
    | rproj < 0 = (-1)
    where
        rproj = (nor >.< r)
        isFZ = ((-dist) <= rproj) && (rproj <= dist)

getInFZ::[(FZPlane,SymmOp)] -> Int -> Rodrigues -> Rodrigues
getInFZ _ 0 r = debug "Warning!!! Numeric Instability or wrong symmetry definition. " r
getInFZ symmInfo limit r = case find (\x -> (snd x) /= 0) result of
    (Just ((_, symmOp), dir)) -> getInFZ symmInfo (limit-1) $ rotate dir symmOp r
    Nothing -> r
    where
        testFZ = map (\x-> isInFZ r $ fst x) symmInfo
        result = zip symmInfo testFZ

toFZ::Symm -> [Rodrigues] -> [Rodrigues]
toFZ symm rs = map toFZ rs
    where
        axisSymmList = getVectorList symm
        opsSymm = map getSymmOp axisSymmList
        fzPlanes = map getFZPlane axisSymmList
        maxNFold = foldr (\(Axis _ n) acc -> max n acc) 1 axisSymmList
        minSphereRadius = foldr (\(Plane _ dist) acc -> max dist acc) 0 fzPlanes
        symmInfo = zip fzPlanes opsSymm
        isInsideSphere (Rodrigues r) = (lengthVector r) <= minSphereRadius
        toFZ r  | isInsideSphere r = r
                | otherwise = getInFZ symmInfo (2*maxNFold) r


-- | Check the correctness of the module
checkSymmetricSystemModule::IO ()
checkSymmetricSystemModule = do
    let
        a = toRodrigues $ AxisPair (1,1,1) (Rad $ pi/3)
        b = toRodrigues $ AxisPair (0,0,1) (Rad $ pi/2)
        c = toRodrigues $ AxisPair (1,0,0) (Rad $ pi)
        d = toRodrigues $ AxisPair (0,-1,-1) (Rad $ pi/2)
        e = toRodrigues $ AxisPair (-1,0,0) (Rad $ pi/4)
        f = toRodrigues $ AxisPair (1,0,0) (Rad 0)
        g = toRodrigues $ Euler (Deg 0.0, Deg 54.7, Deg 45.0)
        h = toRodrigues $ Euler (Deg 30.0, Deg 54.7, Deg 45.0)
        i = toRodrigues $ Euler (Deg 120.0, Deg 54.7, Deg 45.0)
        j = Rodrigues (1.0, 0.64, -0.19)
        test = [a,b,c,d,e,f,g,h,i,j]

    putStrLn  "-- Checking SymmetricSystem module..."
    mapM testFZ test
    -- mapM (\x -> mapM (foldMis x) test) test
    putStrLn $ "--- check end"

    where

        testFZ a = do
            let syList = map getFZPlane (getVectorList Cubic)
                test = map (isInFZ a) syList
                [a'] = toFZ Cubic [a]
                test' = map (isInFZ a') syList
            putStrLn "\n-- Test Correction"
            putStrLn $ ">> " ++ show a ++ " -> " ++ show a'
            putStrLn $ ">> " ++ show ((fromRodrigues a)::AxisPair) ++ " -> " ++ show ((fromRodrigues a')::AxisPair)
            putStrLn $ show test ++ " -> " ++ show test'

        foldMis::Rodrigues -> Rodrigues -> IO Rodrigues
        foldMis a b = do
            putStrLn "\n-- Compose misorientation"
            putStrLn $ ">> a = " ++ show ((fromRodrigues a)::AxisPair)
            putStrLn $ ">> b = " ++ show ((fromRodrigues b)::AxisPair)
            let c = a -@- b
            putStrLn $ ">> a[-1].b (compose) = " ++ show ((fromRodrigues c)::AxisPair)
            return c

mSymm::[RotMatrix]
mSymm = map RotMatrix ls
    where   ls = [
                ((1,0,0),
                 (0,1,0),
                 (0,0,1)),

                ((0,0,1),
                 (0,1,0),
                 (-1,0,0)),

                ((0,-1,0),
                 (1,0,0),
                 (0,0,1)),

                ((0,1,0),
                 (0,0,1),
                 (1,0,0)),

                ((0,0,-1),
                 (1,0,0),
                 (0,-1,0)),

                ((-1,0,0),
                 (0,0,1),
                 (0,1,0)),



                ((0,0,-1),
                 (0,-1,0),
                 (-1,0,0)),

                ((1,0,0),
                 (0,0,-1),
                 (0,1,0)),

                ((-1,0,0),
                 (0,-1,0),
                 (0,0,-1)),

                ((0,0,-1),
                 (-1,0,0),
                 (0,1,0)),

                ((0,0,1),
                 (-1,0,0),
                 (0,-1,0)),

                ((0,0,1),
                 (0,-1,0),
                 (1,0,0)),



                ((0,0,-1),
                 (0,1,0),
                 (1,0,0)),

                ((1,0,0),
                 (0,-1,0),
                 (0,0,-1)),

                ((0,1,0),
                 (-1,0,0),
                 (0,0,1)),

                ((0,-1,0),
                 (0,0,1),
                 (-1,0,0)),

                ((0,-1,0),
                 (0,0,-1),
                 (1,0,0)),

                ((0,-1,0),
                 (-1,0,0),
                 (0,0,-1)),



                 ((-1,0,0),
                 (0,1,0),
                 (0,0,-1)),

                ((1,0,0),
                 (0,0,1),
                 (0,-1,0)),

                ((0,0,1),
                 (1,0,0),
                 (0,1,0)),

                ((0,1,0),
                 (0,0,-1),
                 (-1,0,0)),

                ((0,1,0),
                 (1,0,0),
                 (0,0,-1)),

                ((-1,0,0),
                 (0,0,-1),
                 (0,-1,0))
                ]

test =   mapM (\x -> testRot x (Rodrigues (0.254,-0.235,-0.291))) mSymm
    where
        testRot::RotMatrix -> Rodrigues -> IO ()
        testRot a b = do
            putStrLn "\n-- Test Matrix Symm"
            let c = compM a (fromRodrigues b)
            putStrLn $ " >> By Matrix\n" ++ show c
            putStrLn $ " >> rot =\n" ++ show (a)
            putStrLn $ " >> = " ++ show (toRodrigues c)
            putStrLn $ " >> = " ++ show ((fromRodrigues $ toRodrigues c)::AxisPair)

            let c' = (toRodrigues a) *@* b
            putStrLn $ "\n >> By Frank\n" ++ show c'
            putStrLn $ " >> rot = " ++ show (toRodrigues a)
            putStrLn $ " >> = " ++ show ((fromRodrigues c')::AxisPair)




testRot::(Rodrigues, Rodrigues) -> IO ()
testRot (b, c) = do
    putStrLn "\n-- Test Matrix Symm"
    let d = map (\s -> compM s (fromRodrigues b)) mSymm
        d' = (map (fromRodrigues.toRodrigues) d)::[AxisPair]
    --putStrLn $ " >> By Matrix\n" ++ show c
    putStrLn $ " >> ini =\n" ++ show (b)
    putStrLn $ "->> " ++ show ((fromRodrigues b)::AxisPair)
    putStrLn $ " >> symm = " ++ show d'
    putStrLn $ "->> symm by module" ++ show ((fromRodrigues c)::AxisPair)

















