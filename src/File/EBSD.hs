{-# LANGUAGE RecordWildCards   #-}
module File.EBSD
       ( EBSD
       , toANG
       , toCTF
       , writeANG
       , writeCTF
       , readEBSD
       , readEBSDToVoxBox
       ) where

import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import Control.Exception
import Control.Monad
import GHC.IO.Handle

import Hammer.VoxBox

import File.ANGReader as A
import File.ANGWriter
import File.CTFReader as C
import File.CTFWriter

class EBSD a where
  toANG :: a -> ANGdata
  toCTF :: a -> CTFdata

instance EBSD ANGdata where
  toANG = id
  toCTF = angdata2ctfdata

instance EBSD CTFdata where
  toANG = ctfdata2angdata
  toCTF = id

-- =================================== ANG -> CTF ========================================

angpoint2ctfpoint :: ANGpoint -> CTFpoint
angpoint2ctfpoint p = let
  x = round $ 10 * (A.ci p)
  in CTFpoint
     { C.phase        = A.phaseNum p
     , C.rotation     = A.rotation p
     , C.xpos         = A.xpos     p
     , C.ypos         = A.ypos     p
     , C.bands        = x
     , C.errorID      = x
     , C.angulardev   = A.fit p
     , C.bandcontrast = round $ A.iq p
     , C.bandslope    = round $ A.iq p
     }

angphase2ctfphase :: ANGphase -> CTFphase
angphase2ctfphase p = CTFphase
  { C.phaseID      = A.phase p
  , C.latticeCons  = A.latticeCons p
  , C.materialName = A.materialName p
  , C.phaseInfo    = [A.formula p, A.info p]
  }

anggrid2ctfgrid :: ANGgrid -> CTFgrid
anggrid2ctfgrid g = CTFgrid
  { C.rowCols = let (row, col, _) = A.rowCols g in (row, col)
  , C.xystep  = A.xystep g
  , C.origin  = A.origin g
  }

anginfo2ctfinfo :: ANGinfo -> CTFinfo
anginfo2ctfinfo i = CTFinfo
  { C.project = (A.sampleID i) ++ "(from ANG)"
  , C.jobMode = "Grid"
  , C.author  = A.operator i
  , C.info    = ["workdistace=" ++ show (A.workDist i)]
  }

angdata2ctfdata :: ANGdata -> CTFdata
angdata2ctfdata d = CTFdata
  { C.nodes    = V.map angpoint2ctfpoint (A.nodes d)
  , C.grid     = anggrid2ctfgrid         (A.grid d)
  , C.ebsdInfo = anginfo2ctfinfo         (A.ebsdInfo d)
  , C.phases   = map angphase2ctfphase   (A.phases d)
  }

-- =================================== CTF -> ANG ========================================

ctfpoint2angpoint :: CTFpoint -> ANGpoint
ctfpoint2angpoint p = ANGpoint
  { A.rotation = C.rotation p
  , A.xpos     = C.xpos p
  , A.ypos     = C.ypos p
  , A.iq       = fromIntegral $ C.bandcontrast p
  , A.ci       = (fromIntegral $ C.bands p) / 10
  , A.phaseNum = C.phase p
  , A.detecInt = 1
  , A.fit      = C.angulardev p
  }

ctfphase2angphase :: CTFphase -> ANGphase
ctfphase2angphase p = ANGphase
  { A.phase        = C.phaseID p
  , A.materialName = C.materialName p
  , A.formula      = ""
  , A.info         = unwords (C.phaseInfo p)
  , A.symmetry     = "0"
  , A.latticeCons  = C.latticeCons p
  , A.numFamilies  = 1
  , A.hklFamilies  = [(0, 0, 0, 0, 0, 0)]
  , A.elasticConst = [(0, 0, 0, 0, 0, 0)]
  , A.categories   = (0, 0, 0, 0, 0)
  }

ctfgrid2anggrid :: CTFgrid -> ANGgrid
ctfgrid2anggrid g = ANGgrid
  { A.rowCols = let (row, col) = C.rowCols g in (row, col, col)
  , A.xystep  = C.xystep g
  , A.origin  = C.origin g
  , A.hexGrid = False
  }

ctfinfo2anginfo :: CTFinfo -> ANGinfo
ctfinfo2anginfo i = ANGinfo
  { workDist = 16
  , pixperum = 1
  , operator = C.author i
  , sampleID = C.project i
  , scanID   = ""
  }

ctfdata2angdata :: CTFdata -> ANGdata
ctfdata2angdata d = ANGdata
  { A.nodes    = V.map ctfpoint2angpoint (C.nodes d)
  , A.grid     = ctfgrid2anggrid         (C.grid d)
  , A.ebsdInfo = ctfinfo2anginfo         (C.ebsdInfo d)
  , A.phases   = map ctfphase2angphase   (C.phases d)
  }

-- =================================== Reader ========================================

readEBSD :: FilePath -> IO (Either ANGdata CTFdata)
readEBSD f = let
  msg e1 e2 = "This file is neither a valid ANG file nor a valid CTF file.\n" ++ e1 ++ "\n" ++ e2
  in catch (liftM Left (parseANG f)) (\eANG -> do
     let errANG = show (eANG :: SomeException)
     catch (liftM Right (parseCTF f)) (\eCTF -> do
        let errCTF = show (eCTF :: SomeException)
        error (msg errANG errCTF)
        )
  )

readEBSDToVoxBox :: (U.Unbox a)=> (CTFpoint -> a) -> (ANGpoint -> a) -> Either ANGdata CTFdata -> VoxBox a
readEBSDToVoxBox fctf fang = either (either error id . (flip A.angToVoxBox) fang) ((flip C.ctfToVoxBox) fctf)

writeANG :: (EBSD a)=> FilePath -> a -> IO ()
writeANG f a = renderANGFile f (toANG a)

writeCTF :: (EBSD a)=> FilePath -> a -> IO ()
writeCTF f a = renderCTFFile f (toCTF a)
