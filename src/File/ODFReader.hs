{-# LANGUAGE
    NamedFieldPuns
  , OverloadedStrings
  #-}

module File.ODFReader
  ( parseDiscODF
  , getPos
  , getLinPos
  ) where

import qualified Data.Vector as V

import Data.Vector         (Vector)
import Control.Applicative ((<|>))
import Data.ByteString     (ByteString)

import Data.Attoparsec.ByteString.Char8

data DiscODF =
  DiscODF
  { step  :: Double
  , nPHI1 :: Int
  , nPHI  :: Int
  , nPHI2 :: Int
  , sPHI1 :: Double
  , sPHI  :: Double
  , sPHI2 :: Double
  , odf   :: Vector Double
  } deriving (Show)

getLinPos :: DiscODF -> (Int, Int, Int) -> Int
getLinPos df = getLinPos' (nPHI1 df) (nPHI df)

getLinPos' :: Int -> Int -> (Int, Int, Int) -> Int
getLinPos' nPHI1 nPHI (i,j,k) = k * nPHI * nPHI1 + j * nPHI1 + i

getPos :: DiscODF -> Int -> (Int, Int, Int)
getPos df n = (i,j,k)
  where
    (k, rk) = quotRem n $ (nPHI1 df) * (nPHI df)
    (j,  i) = quotRem rk (nPHI1 df)

-- | Parse just one input file. Rise an error mesage in case of error.
parseDiscODF :: ByteString -> DiscODF
parseDiscODF odf =
  case parseOnly parseAll odf of
   Left err -> error ("[ODFReader] Error reading discrete ODF: " ++ show err)
   Right xs -> xs

parseAll :: Parser DiscODF
parseAll = do
  _     <- manyTill anyChar endOfLine
  step  <- func "Distance between grid points"
  nPHI1 <- func "Number of PHI1-values"
  sPHI1 <- func "First value of PHI1"
  nPHI  <- func "Number of PHI2-values"
  sPHI  <- func "First value of PHI2"
  nPHI2 <- func "Number of PHI-values"
  sPHI2 <- func "First value of PHI"
  grid  <- gridParse nPHI1 nPHI nPHI2
  return DiscODF { step, nPHI1, nPHI, nPHI2, sPHI1, sPHI, sPHI2, odf = grid }
  where func s = do { x <- parseNum; skipSpace >> string s >> endOfLine; return x }

gridParse :: Int -> Int -> Int -> Parser (Vector Double)
gridParse nPHI1 nPHI nPHI2 = do
  grid <- V.fromList <$> many1 parseNum
  if V.length ls == V.length grid
    then return $ newArr grid
    else error $ "[ODFReader] Unexpected number of points: " ++ show ( V.length ls, V.length grid)
  where
    newArr x = V.update (V.replicate (V.length ns) 0) (V.zip ns x)
    ns = V.map (getLinPos' nPHI1 nPHI) ls
    ls = V.fromList [(i,j,k) | j <- [0..nPHI-1], k <- [0..nPHI2-1], i <- [0..nPHI1-1] ]

parseNum :: (Read a) => Parser a
parseNum = try $ do
  skipSpace
  n <- many1 (digit <|> char '-' <|> char '.')
  -- correct double represantion when -1 < n < 1 (.xxx or -.xxx)
  let n' = case n of
        ('.':_) -> '0':n
        ('-':'.':xs) -> '-':'0':'.':xs
        _ -> n
  return (read n')
