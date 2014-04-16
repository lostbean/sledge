{-# OPTIONS_GHC -fno-warn-unused-do-bind #-}
{-# LANGUAGE RecordWildCards   #-}
{-# LANGUAGE OverloadedStrings #-}

-- | Read and load *.ang files from EBSD measure systems.
module File.ANGReader
       ( parseANG
       , ebsdToVoxBox
       , EBSDpoint (..)
       , EBSDinfo  (..)
       , EBSDphase (..)
       , Gridinfo  (..)
       , EBSDdata  (..)
       ) where

import qualified Data.Vector                      as V
import qualified Data.Vector.Unboxed              as U
import qualified Data.ByteString                  as B
import qualified Data.ByteString.Char8            as BC
import qualified Data.Attoparsec.ByteString.Char8 as AC

import           Hammer.VoxBox

import           Texture.Orientation (Quaternion, toQuaternion, mkEuler, Rad(..))
import           Data.Vector         (Vector)
import           Control.Applicative ((<|>), (<$>), pure)

import           Data.Attoparsec.ByteString.Char8

-- ============================== ANG manipulation =======================================

ebsdToVoxBox :: (U.Unbox a)=> EBSDdata -> (EBSDpoint -> a) -> VoxBox a
ebsdToVoxBox EBSDdata{..} func = let
  EBSDinfo{..}       = ebsdInfo
  Gridinfo{..}       = grid
  (xstep, ystep)     = xystep
  (row, col_even, _) = rowCols
  boxorg = VoxelPos 0 0 0
  boxdim = VoxBoxDim col_even row 1
  dime   = VoxBoxRange boxorg boxdim
  org    = VoxBoxOrigin xstart ystart zstart
  spc    = VoxelDim xstep ystep ((xstep + ystep) / 2)
  in VoxBox dime org spc (V.convert $ V.map func nodes)

-- ============================== Data definition ========================================

-- | Information associated to each point. For futher reference consult OIM manual.
data EBSDpoint
  = EBSDpoint
  { rotation :: {-# UNPACK #-} !Quaternion
  , xpos     :: {-# UNPACK #-} !Int
  , ypos     :: {-# UNPACK #-} !Int
  , qi       :: {-# UNPACK #-} !Double
  , ci       :: {-# UNPACK #-} !Double
  , phaseNum :: {-# UNPACK #-} !Int
  , detecInt :: {-# UNPACK #-} !Int
  } deriving Show

-- | Information describing the measuriment.
data EBSDinfo =
  EBSDinfo
  { xstart       :: Double
  , ystart       :: Double
  , zstart       :: Double
  , workDist     :: Double
  } deriving (Show)

data EBSDphase =
  EBSDphase
  { phase        :: Int
  , materialName :: String
  , formula      :: String
  , info         :: String
  , symmetry     :: String
  , latticeCons  :: (Double, Double, Double, Double, Double, Double)
  , numFamilies  :: Int
  , hklFamilies  :: [(Int, Int, Int, Int, Double, Int)]
  , elasticConst :: [(Double, Double, Double, Double, Double, Double)]
  , categories   :: (Int, Int, Int, Int, Int)
  } deriving Show

-- | Information about the grid of point. Hexagonal or Square
data Gridinfo
  = Gridinfo
  { rowCols :: (Int, Int, Int)
  , xystep  :: (Double, Double)
  , hexGrid :: Bool
  } deriving Show

-- | Hold the whole ANG data strcuture
data EBSDdata
  = EBSDdata
  { nodes    :: Vector EBSDpoint
  , grid     :: Gridinfo
  , ebsdInfo :: EBSDinfo
  , phases   :: [EBSDphase]
  } deriving Show

-- ============================== Parse functions ========================================

-- | Read the input ANG file. Rise an error mesage in case of bad format or acess.
parseANG :: String -> IO EBSDdata
parseANG fileName = do
  stream <- B.readFile fileName
  case parseOnly parseEBSDdata stream of
    Left err -> error (">> Error reading file " ++ fileName ++ "\n" ++ show err)
    Right xs -> return xs

parseEBSDdata :: Parser EBSDdata
parseEBSDdata = do
  head_info   <- headParse
  phases_info <- many1 (skipEmptyComment >> phaseParse)
  grid_info   <-        skipEmptyComment >> gridParse
  skipMany (skipComment)
  point_list  <- many1 (pointParse grid_info)
  return $ EBSDdata (V.fromList point_list) grid_info head_info phases_info

-- ------------------------------------ SubParsers ---------------------------------------

headParse :: Parser EBSDinfo
headParse = do
  _ <- getInfo "# TEM_PIXperUM"    parseFloat
  x <- getInfo "# x-star"          parseFloat
  y <- getInfo "# y-star"          parseFloat
  z <- getInfo "# z-star"          parseFloat
  w <- getInfo "# WorkingDistance" parseFloat
  return $ EBSDinfo x y z w

phaseParse :: Parser EBSDphase
phaseParse = do
  phn <- getInfo "# Phase"        parseInt
  mat <- getInfo "# MaterialName" parseText
  for <- getInfo "# Formula"      parseText
  inf <- getInfo "# Info"         parseText
  sym <- getInfo "# Symmetry"     parseText
  lat <- latticeParse
  fam <- option 0  $ getInfo "# NumberFamilies" parseInt
  hkl <- option [] $ many1 hklParse
  ela <- option [] $ many1 elasticParse
  cat <- categoryParse
  return $ EBSDphase phn mat for inf sym lat fam hkl ela cat

gridParse :: Parser Gridinfo
gridParse = do
  isHex <- getInfo "# GRID:"       gridType
  xstep <- getInfo "# XSTEP:"      parseFloat
  ystep <- getInfo "# YSTEP:"      parseFloat
  codd  <- getInfo "# NCOLS_ODD:"  parseInt
  ceven <- getInfo "# NCOLS_EVEN:" parseInt
  row   <- getInfo "# NROWS:"      parseInt
  return $ Gridinfo (row, ceven, codd) (xstep, ystep) isHex

elasticParse :: Parser (Double, Double, Double, Double, Double, Double)
elasticParse = "# ElasticConstants" .*> do
  i1 <- parseFloat
  i2 <- parseFloat
  i3 <- parseFloat
  i4 <- parseFloat
  i5 <- parseFloat
  i6 <- parseFloat
  eol
  return (i1, i2, i3, i4, i5, i6)

latticeParse :: Parser (Double, Double, Double, Double, Double, Double)
latticeParse = "# LatticeConstants" .*> do
  a  <- parseFloat
  b  <- parseFloat
  c  <- parseFloat
  w1 <- parseFloat
  w2 <- parseFloat
  w3 <- parseFloat
  eol
  return (a, b, c, w1, w2, w3)

hklParse :: Parser (Int, Int, Int, Int, Double, Int)
hklParse = "# hklFamilies" .*> do
  h <- parseInt
  k <- parseInt
  l <- parseInt
  a <- parseInt
  -- TODO Fix
  b <- (blanks >> manyTill anyChar space >> return 1)
  c <- parseInt <|> return 0
  eol
  return (h, k, l, a, b, c)

categoryParse :: Parser (Int,Int,Int,Int,Int)
categoryParse = "# Categories" .*> do
  a <- parseInt
  b <- parseInt
  c <- parseInt
  d <- parseInt
  f <- parseInt
  eol
  return (a,b,c,d,f)

pointParse :: Gridinfo -> Parser EBSDpoint
pointParse g = do
  p1 <- parseFloat
  p  <- parseFloat
  p2 <- parseFloat
  x  <- parseFloat
  y  <- parseFloat
  q  <- parseFloat
  c  <- parseFloat
  ph <- parseInt
  dc <- parseInt
  _  <- parseFloat
  eol
  let rot     = toQuaternion $ mkEuler (Rad p1) (Rad p) (Rad p2)
      (xi,yi) = getGridPoint g (x, y)
  return $ EBSDpoint rot xi yi q c ph dc

-- | Calculate colunm position (row,col) from ID sequence
-- Origin at (1,1)
getGridPoint :: Gridinfo -> (Double, Double) -> (Int, Int)
getGridPoint g (x, y) = let
  (xstep, ystep) = xystep g
  yi = round ((2*y)/ystep) `div` 2
  -- divide for 0.5*xstep to avoid error when round value
  -- ex. round 7.4/5.0 \= round 7.6/5.0;
  xi = round ((2*x)/xstep) `div` 2
  in (xi, yi)

getInfo :: B.ByteString -> Parser a -> Parser a
getInfo ident func = let
  parser = do
    string ident
    x <- func
    eol
    return x
  in parser <|> do
    found <- BC.unpack <$> AC.take (B.length ident)
    fail ("Can't find the following field -> " ++ BC.unpack ident ++ ", found -> " ++ found )

-- -------------------------------------- Basic parsers ----------------------------------

skipEmptyComment :: Parser ()
skipEmptyComment = getInfo "#" (pure ())

skipComment :: Parser ()
skipComment = getInfo "#" (skipWhile (not . isEOL))

gridType :: Parser Bool
gridType = blanks >> (getInfo "HexGrid" (pure True)) <|> (getInfo "SqrGrid" (pure False))

parseText :: Parser String
parseText = blanks >> ((unwords . words . BC.unpack) <$> AC.takeWhile (not . isEOL))

parseFloat :: Parser Double
parseFloat = blanks >> signed double

parseInt :: Parser Int
parseInt = blanks >> signed decimal

-- | Skips blank chars (space and tab)
blanks :: Parser ()
blanks = skipWhile (\c -> c == ' ' || c == '\t')

isEOL :: Char -> Bool
isEOL c = c == '\r' || c == '\n'

-- | Skips blanks chars till and including the eol (End Of Line - CR-LF or LF)
eol :: Parser ()
eol = blanks >> skipWhile isEOL
