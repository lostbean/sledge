{-# LANGUAGE RecordWildCards #-}

module TestTexture
       ( runTextureChecker
       ) where

import qualified Data.Vector as V

import           Data.Vector (Vector)

import           Control.Applicative
import           Hammer.Math.Algebra
import           Texture.Orientation
import           Texture.Symmetry
import           Test.QuickCheck

import           Debug.Trace

-- ================================================================================

dbg :: (Show a1, Show a)=> a1 -> a -> a
dbg x a = trace (show x ++ " " ++ show a) a

runTextureChecker :: IO ()
runTextureChecker =  do
  let myArgs = Args { replay          = Nothing
                    , maxSuccess      = 1000
                    , maxDiscardRatio = 50
                    , maxSize         = 1000
                    , chatty          = True }
  putStrLn "Testing Orientation stuffs .."
  quickCheckWith myArgs testFundamentalZone
  quickCheckWith myArgs testFundamentalZoneMatrix

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
  fz_r  = V.map (fromQuaternion . symmOp) fz :: Vector Rodrigues
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
