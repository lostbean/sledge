module BinghamConstantTest (
    checkDerivatives,
) where

import qualified Data.Vector.Unboxed as U
import Texture.Bingham.Constant (dC, saddleC3)

checkDerivatives :: U.Vector Double -> Double
checkDerivatives xs =
    let
        gends n = U.generate size (\i -> if i == n then 1 else 0)
        size = U.length xs
        c1 = sum [dC (gends i) xs | i <- [0 .. size - 1]]
        c2 = saddleC3 xs
     in
        (c2 + c1) / c1
