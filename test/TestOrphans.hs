{-# LANGUAGE FlexibleInstances #-}
{-# OPTIONS_GHC -Wno-orphans #-}

module TestOrphans where

import Control.Applicative
import Test.QuickCheck

import Linear.Vect
import Texture.Orientation

instance Arbitrary Euler where
    arbitrary = liftA3 mkEuler x2 x1 x2
      where
        x1 = Deg <$> choose (0, 180)
        x2 = Deg <$> choose (0, 360)

instance Arbitrary Vec3D where
    arbitrary = normalize <$> liftA3 Vec3 p p p
      where
        p = choose (0, 1)

instance Arbitrary Quaternion where
    arbitrary = toQuaternion <$> (arbitrary :: Gen Euler)

instance Arbitrary Rodrigues where
    arbitrary = fromQuaternion <$> arbitrary

instance Arbitrary Deg where
    arbitrary = Deg <$> arbitrary
