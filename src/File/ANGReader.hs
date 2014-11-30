{-# OPTIONS_GHC -fno-warn-unused-do-bind #-}
{-# LANGUAGE RecordWildCards   #-}
{-# LANGUAGE OverloadedStrings #-}

-- | Read and load *.ang files from ANG measure systems.
module File.ANGReader
       ( parseANG
       , angToVoxBox
       , ANGpoint (..)
       , ANGinfo  (..)
       , ANGphase (..)
       , ANGgrid  (..)
       , ANGdata  (..)
       ) where

import qualified Data.Vector                      as V
import qualified Data.Vector.Unboxed              as U
import qualified Data.ByteString                  as B
import qualified Data.ByteString.Char8            as BC
import qualified Data.Attoparsec.ByteString.Char8 as AC

import           Hammer.VoxBox

import           Texture.Orientation (Quaternion, toQuaternion, mkEuler, Rad(..))
import           Data.Vector         (Vector)
import           Control.Applicative ((<|>), (<$>), (*>), pure)

import           Data.Attoparsec.ByteString.Char8

-- ============================== ANG manipulation =======================================

angToVoxBox :: (U.Unbox a)=> ANGdata -> (ANGpoint -> a) -> Either String (VoxBox a)
angToVoxBox ANGdata{..} func
  | not hexGrid = Right box
  | otherwise   = Left "[ANG] Can't produce VoxBox from ANG file with hexagonal grid."
  where
    ANGinfo{..} = ebsdInfo
    ANGgrid{..} = grid
    (xstep, ystep) = xystep
    (row, col_even, _) = rowCols
    (xstart, ystart, zstart) = origin
    boxorg = VoxelPos 0 0 0
    boxdim = VoxBoxDim col_even row 1
    dime   = VoxBoxRange boxorg boxdim
    org    = VoxBoxOrigin xstart ystart zstart
    spc    = VoxelDim xstep ystep ((xstep + ystep) / 2)
    box    = VoxBox dime org spc (V.convert $ V.map func nodes)

-- ============================== Data definition ========================================

-- | Information associated to each point. For futher reference consult OIM manual.
data ANGpoint
  = ANGpoint
  { rotation :: {-# UNPACK #-} !Quaternion
  , xpos     :: {-# UNPACK #-} !Int
  , ypos     :: {-# UNPACK #-} !Int
  , iq       :: {-# UNPACK #-} !Double
  , ci       :: {-# UNPACK #-} !Double
  , phaseNum :: {-# UNPACK #-} !Int
  , detecInt :: {-# UNPACK #-} !Int
  , fit      :: {-# UNPACK #-} !Double
  } deriving Show

-- | Information describing the measuriment.
data ANGinfo =
  ANGinfo
  { workDist :: Double
  , pixperum :: Double
  , operator :: String
  , sampleID :: String
  , scanID   :: String
  } deriving Show

data ANGphase =
  ANGphase
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
data ANGgrid
  = ANGgrid
  { rowCols :: (Int, Int, Int)
  , xystep  :: (Double, Double)
  , origin  :: (Double, Double, Double)
  , hexGrid :: Bool
  } deriving Show

-- | Hold the whole ANG data strcuture
data ANGdata
  = ANGdata
  { nodes    :: Vector ANGpoint
  , grid     :: ANGgrid
  , ebsdInfo :: ANGinfo
  , phases   :: [ANGphase]
  } deriving Show

-- ============================== Parse functions ========================================

-- | Read the input ANG file. Rise an error mesage in case of bad format or acess.
parseANG :: String -> IO ANGdata
parseANG fileName = do
  stream <- B.readFile fileName
  let
    ebsd = case parseOnly parseANGdata stream of
      Left err -> error ("[ANG] Error reading file " ++ fileName ++ "\n" ++ show err)
      Right x  -> x
  checkDataShape ebsd
  checkDataSize  ebsd
  return ebsd

parseANGdata :: Parser ANGdata
parseANGdata = do
  p <- getInfo "# TEM_PIXperUM"    parseFloat
  r <- origParse
  w <- getInfo "# WorkingDistance" parseFloat
  phases_info <- many1 (skipCommentLine >> phaseParse)
  grid_info   <- skipCommentLine >> gridParse r
  o <- skipCommentLine >> getInfo "# OPERATOR:" parseText
  s <- skipCommentLine >> getInfo "# SAMPLEID:" parseText
  c <- skipCommentLine >> getInfo "# SCANID:"   parseText
  skipCommentLine
  let head_info = ANGinfo w p o s c
  point_list  <- many1 (pointParse grid_info)
  return $ ANGdata (V.fromList point_list) grid_info head_info phases_info

checkDataShape :: ANGdata -> IO ()
checkDataShape ANGdata{..}
  | not hexGrid && cEven == cOdd = putStrLn "[ANG] Using square grid."
  | not hexGrid = msg
  | otherwise   = putStrLn "[ANG] Using hexagonal grid."
  where
    msg = error $ "[ANG] Improper square grid shape. (row, cEven, cOdd) = " ++ show rowCols
    ANGgrid{..}     = grid
    (_, cEven, cOdd) = rowCols

checkDataSize :: ANGdata -> IO ()
checkDataSize ANGdata{..}
  | V.length nodes == n = putStrLn $ "[ANG] Loaded numbers of points: " ++ show n
  | otherwise = msg
  where
    ANGgrid{..}       = grid
    (row, cEven, cOdd) = rowCols
    (halfrow, leftrow) = row `quotRem` 2
    n = halfrow * cEven + (halfrow + leftrow) * cOdd
    msg = error $ "[ANG] Invalid numbers of points. Expected " ++
          show n ++ " but found " ++
          show (V.length nodes) ++ ". The grid shape is given by (row, cEven, cOdd) = " ++
          show rowCols

-- ------------------------------------ SubParsers ---------------------------------------

phaseParse :: Parser ANGphase
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
  return $ ANGphase phn mat for inf sym lat fam hkl ela cat

origParse :: Parser (Double, Double, Double)
origParse = do
  x <- getInfo "# x-star" parseFloat
  y <- getInfo "# y-star" parseFloat
  z <- getInfo "# z-star" parseFloat
  return (x, y, z)

gridParse :: (Double, Double, Double) -> Parser ANGgrid
gridParse orig = do
  isHex <- getInfo "# GRID:"       gridType
  xstep <- getInfo "# XSTEP:"      parseFloat
  ystep <- getInfo "# YSTEP:"      parseFloat
  codd  <- getInfo "# NCOLS_ODD:"  parseInt
  ceven <- getInfo "# NCOLS_EVEN:" parseInt
  row   <- getInfo "# NROWS:"      parseInt
  return $ ANGgrid (row, ceven, codd) (xstep, ystep) orig isHex

elasticParse :: Parser (Double, Double, Double, Double, Double, Double)
elasticParse = "# ElasticConstants" *> do
  i1 <- parseFloat
  i2 <- parseFloat
  i3 <- parseFloat
  i4 <- parseFloat
  i5 <- parseFloat
  i6 <- parseFloat
  eol
  return (i1, i2, i3, i4, i5, i6)

latticeParse :: Parser (Double, Double, Double, Double, Double, Double)
latticeParse = "# LatticeConstants" *> do
  a  <- parseFloat
  b  <- parseFloat
  c  <- parseFloat
  w1 <- parseFloat
  w2 <- parseFloat
  w3 <- parseFloat
  eol
  return (a, b, c, w1, w2, w3)

hklParse :: Parser (Int, Int, Int, Int, Double, Int)
hklParse = "# hklFamilies" *> do
  h <- parseInt
  k <- parseInt
  l <- parseInt
  a <- parseInt
  b <- parseFloat
  c <- parseInt
  eol
  return (h, k, l, a, b, c)

categoryParse :: Parser (Int,Int,Int,Int,Int)
categoryParse = "# Categories" *> do
  a <- parseInt
  b <- parseInt
  c <- parseInt
  d <- parseInt
  f <- parseInt
  eol
  return (a,b,c,d,f)

pointParse :: ANGgrid -> Parser ANGpoint
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
  f  <- parseFloat <|> return 0
  eol
  let rot     = toQuaternion $ mkEuler (Rad p1) (Rad p) (Rad p2)
      (xi,yi) = getGridPoint g (x, y)
  return $ ANGpoint rot xi yi q c ph dc f

-- | Calculate colunm position (row,col) from ID sequence
-- Origin at (1,1)
getGridPoint :: ANGgrid -> (Double, Double) -> (Int, Int)
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

skipCommentLine :: Parser ()
skipCommentLine = getInfo "#" (skipWhile (not . isEOL))

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
