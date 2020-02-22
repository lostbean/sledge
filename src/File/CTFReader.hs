{-# LANGUAGE
    RecordWildCards
  , OverloadedStrings
  #-}
-- | Read and load *.ctf files from CTF measure systems.
module File.CTFReader
  ( loadCTF
  , readFileCTF
  , ctfToVoxBox
  , CTFpoint (..)
  , CTFinfo  (..)
  , CTFphase (..)
  , CTFgrid  (..)
  , CTFdata  (..)
  ) where

import Control.Applicative ((<|>))
import Data.Attoparsec.Lazy
import Data.Vector (Vector)
import Hammer.VoxBox
import qualified Data.Attoparsec.Lazy       as A
import qualified Data.ByteString.Lazy       as BSL
import qualified Data.ByteString.Internal   as BS (w2c)
import qualified Data.Text                  as T
import qualified Data.Text.Encoding         as TE
import qualified Data.Vector                as V
import qualified Data.Vector.Unboxed        as U

import Texture.Orientation (Quaternion, toQuaternion, mkEuler, Deg(..))
import File.InternalParsers

-- ============================== CTF manipulation =======================================

ctfToVoxBox :: (U.Unbox a)=> CTFdata -> (CTFpoint -> a) -> VoxBox a
ctfToVoxBox CTFdata{..} func = box
  where
    CTFinfo{..}    = ebsdInfo
    CTFgrid{..}    = grid
    (xstep, ystep) = xystep
    (row, col)     = rowCols
    boxorg = VoxelPos 0 0 0
    boxdim = VoxBoxDim col row 1
    dime   = VoxBoxRange boxorg boxdim
    org    = VoxBoxOrigin 0 0 0
    spc    = VoxelDim xstep ystep ((xstep + ystep) / 2)
    box    = VoxBox dime org spc (V.convert $ V.map func nodes)

-- ============================== Data definition ========================================

-- | Information associated to each point.
-- Relation with OIM parameters:
-- CI  <-> (number of bands / max. number of bands)* or (error / max. error)
-- IQ  <-> band contrast* or band slope
-- fit <-> mad (mean angular deviation)
--
-- OBS. For hexagonal grid a rotation of 30 deg (phi 2) may be necessary. The standard
-- options are marked with (*).
data CTFpoint
  = CTFpoint
  { phase        :: {-# UNPACK #-} !Int     -- ^ phase id (0 means non-indexed)
  , rotation     :: {-# UNPACK #-} !Quaternion
  , xpos         :: {-# UNPACK #-} !Int     -- ^ x position in um
  , ypos         :: {-# UNPACK #-} !Int     -- ^ y position in um
  , bands        :: {-# UNPACK #-} !Int     -- ^ number of bands found
  , errorID      :: {-# UNPACK #-} !Int     -- ^ error
  , angulardev   :: {-# UNPACK #-} !Double  -- ^ MAD (mean angular deviation)
  , bandcontrast :: {-# UNPACK #-} !Int     -- ^ BC  (band contrast)
  , bandslope    :: {-# UNPACK #-} !Int     -- ^ BS  (band slope)
  } deriving Show

-- | Information describing the measuriment.
data CTFinfo =
  CTFinfo
  { project :: String
  , jobMode :: String
  , author  :: String
  , info    :: [String]
  } deriving (Show)

data CTFphase =
  CTFphase
  { phaseID      :: Int
  , latticeCons  :: (Double, Double, Double, Double, Double, Double)
  , materialName :: String
  , phaseInfo    :: [String]
  } deriving Show

-- | Information about the grid of point. Hexagonal or Square
data CTFgrid
  = CTFgrid
  { rowCols :: (Int, Int)
  , xystep  :: (Double, Double)
  , origin  :: (Double, Double, Double)
  } deriving Show

-- | Hold the whole CTF data strcuture
data CTFdata
  = CTFdata
  { nodes    :: Vector CTFpoint
  , grid     :: CTFgrid
  , ebsdInfo :: CTFinfo
  , phases   :: [CTFphase]
  } deriving Show

-- ============================== Parse functions ========================================

-- | Read the input CTF file. Rise an error mesage in case of bad format or acess.
loadCTF :: BSL.ByteString -> Either String CTFdata
loadCTF bs = do
  ebsd <- eitherResult (parse parseCTFdata bs)
  _ <- checkDataSize ebsd
  return ebsd

-- | Read the input CTF file. Rise an error mesage in case of bad format or acess.
readFileCTF :: String -> IO (Either String CTFdata)
readFileCTF fileName = loadCTF <$> BSL.readFile fileName

parseCTFdata :: Parser CTFdata
parseCTFdata = do
  _ <- stringInfo "Channel Text File"
  _proj   <- getInfo "Prj"      parseText
  _auth   <- getInfo "Author"   parseText
  _job    <- getInfo "JobMode"  parseText
  _grid   <- gridParse
  _info   <- parseFields
  _phases <- many1' phaseParse
  skipWhile (not . isEOL . BS.w2c) >> eol
  _points <- many1' (pointParse _grid)
  let _head = CTFinfo _proj _job _auth _info
  return $ CTFdata (V.fromList _points) _grid _head _phases

checkDataSize :: CTFdata -> Either String Int
checkDataSize CTFdata{..}
  | V.length nodes == n = Right n
  | otherwise = Left msg
  where
    CTFgrid{..} = grid
    (row, col)  = rowCols
    n = row * col
    msg = "[CTF] Invalid numbers of points. Expected " ++
          show n ++ " but found " ++
          show (V.length nodes) ++ ". The grid shape is given by (row, cEven, cOdd) = " ++
          show rowCols

-- ------------------------------------ SubParsers ---------------------------------------

phaseParse :: Parser CTFphase
phaseParse = parser <?> "Failed to parse phase info"
  where
    parser = do
      phn <- getInfo "Phases" parseInt
      eol
      lat <- latticeParse
      mat <- parseField
      rs  <- parseFields
      return $ CTFphase phn lat mat rs

gridParse :: Parser CTFgrid
gridParse = parser <?> "Failed to parse grid info"
  where
    parser = do
      row   <- getInfo "XCells" parseInt
      col   <- getInfo "YCells" parseInt
      xstep <- getInfo "XStep"  parseFloat
      ystep <- getInfo "YStep"  parseFloat
      x     <- getInfo "AcqE1"  parseFloat
      y     <- getInfo "AcqE2"  parseFloat
      z     <- getInfo "AcqE3"  parseFloat
      return $ CTFgrid (row, col) (xstep, ystep) (x, y, z)

latticeParse :: Parser (Double, Double, Double, Double, Double, Double)
latticeParse = parser <?> "Failed to parse lattice parameters"
  where
    parser = do
      a  <- parseFloat
      _  <- stringInfo ";"
      b  <- parseFloat
      _  <- stringInfo ";"
      c  <- parseFloat
      blanks
      w1 <- parseFloat
      _  <- stringInfo ";"
      w2 <- parseFloat
      _  <- stringInfo ";"
      w3 <- parseFloat
      eol
      return (a, b, c, w1, w2, w3)

pointParse :: CTFgrid -> Parser CTFpoint
pointParse g = parser <?> "Failed to parse EBSD point"
  where
    parser = do
      ph   <- parseInt
      x    <- parseFloat
      y    <- parseFloat
      bd   <- parseInt
      er   <- parseInt
      phi1 <- parseFloat
      phi  <- parseFloat
      phi2 <- parseFloat
      ma   <- parseFloat
      bc   <- parseInt
      bs   <- parseInt
      eol
      let rot     = toQuaternion $ mkEuler (Deg phi1) (Deg phi) (Deg phi2)
          (xi,yi) = getGridPoint g (x, y)
      return $ CTFpoint ph rot xi yi bd er ma bc bs

-- | Calculate colunm position (row,col) from ID sequence
-- Origin at (1,1)
getGridPoint :: CTFgrid -> (Double, Double) -> (Int, Int)
getGridPoint g (x, y) = let
  (xstep, ystep) = xystep g
  yi = round ((2*y)/ystep) `div` 2
  -- divide for 0.5*xstep to avoid error when round value
  -- ex. round 7.4/5.0 \= round 7.6/5.0;
  xi = round ((2*x)/xstep) `div` 2
  in (xi, yi)

-- -------------------------------------- Basic parsers ----------------------------------

parseFields :: Parser [String]
parseFields = do
  x <- manyTill' parseField (string "\r" <|> string "\n")
  eol
  return x

parseField :: Parser String
parseField = let
  isTheEnd = (\c -> not (isEOL c) && not (isEOField c)) . BS.w2c
  in blanks >> (T.unpack . TE.decodeUtf8 <$> A.takeWhile isTheEnd)

isEOField :: Char -> Bool
isEOField = (== '\t')
