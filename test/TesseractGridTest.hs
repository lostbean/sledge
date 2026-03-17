module TesseractGridTest (
    test,
) where

import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import Hammer.VTK
import System.Random (randomRIO)

import Texture.Orientation (Quaternion)
import Texture.TesseractGrid

test :: IO Bool
test = do
    let m = 10
    c <- randomRIO (1, 4)
    n <- randomRIO (0, m * m * m - 1)
    let
        t0 = emptyTesseract m 0
        p = getTesseractPoint m (c, ixToTessPos m n)
        q = tesseractToQuaternion p
        t1 :: TesseractGrid Double
        t1 = binningTesseract (V.singleton q) t0
        cs = V.fromList [cell1 t1, cell2 t1, cell3 t1, cell4 t1]
        cm = V.maxIndex (V.map V.maximum cs)
        nm = V.maxIndex (cs V.! cm)
        qm = tesseractToQuaternion $ maxTesseractPoint t1
        vtk = plotTesseractPoints (U.fromList [q, qm])
    printTesseract t1 "tessTest"
    writeUniVTKfile "/home/edgar/Desktop/tessTest.vtu" True vtk
    print (q, qm)
    putStr "(cell, linear index, pos) = " >> print (c, n, ixToTessPos m n)
    putStr "TessPoint = " >> print p
    putStr "TessPoint(recalc.) = " >> print (quaternionToTesseract q)
    putStr "linear index(recalc) = " >> print (getTesseractPos m p)
    putStr "(cell, linear index)(recalc) = " >> print (cm + 1, nm)
    return $ c == (cm + 1) && n == nm && q == qm
