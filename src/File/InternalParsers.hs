{-# OPTIONS_GHC -fno-warn-unused-do-bind #-}
{-# LANGUAGE OverloadedStrings #-}

module File.InternalParsers where

import qualified Data.ByteString                  as B
import qualified Data.ByteString.Char8            as BC
import qualified Data.Attoparsec.ByteString.Char8 as AC

import           Control.Applicative ((<|>), (<$>))

import           Data.Attoparsec.ByteString.Char8

getInfo :: B.ByteString -> Parser a -> Parser a
getInfo ident func = do
  stringInfo ident
  x <- func
  eol
  return x

stringInfo :: BC.ByteString -> Parser BC.ByteString
stringInfo s = AC.string s <|> do
    found <- BC.unpack <$> AC.take (BC.length s)
    fail ("expecting -> " ++ BC.unpack s ++ ", but found -> " ++ found)

parseText :: Parser String
parseText = blanks >> ((unwords . words . BC.unpack) <$> AC.takeWhile (not . isEOL))

parseFloat :: Parser Double
parseFloat = parseNumber (signed double <?> "Invalid float number")

parseInt :: Parser Int
parseInt = parseNumber (signed decimal <?> "Invalid integer number")

parseNumber :: Show a => Parser a -> Parser a
parseNumber parser = do
  blanks
  num  <- parser
  next <- AC.peekChar
  let func x = if (isBlank x || isEOL x) then return num else fail ("Invalid number -> " ++ show num ++ show x)
  maybe (return num) func next

-- | Skips blank chars (space and tab)
blanks :: Parser ()
blanks = skipWhile isBlank

skipRestOfTheLine :: Parser ()
skipRestOfTheLine = skipWhile (not . isEOL) >> skipWhile isEOL

isBlank :: Char -> Bool
isBlank c = c == ' ' || c == '\t'

isEOL :: Char -> Bool
isEOL c = c == '\r' || c == '\n'

-- | Skips blanks chars till and including the eol (End Of Line - CR-LF or LF)
eol :: Parser ()
eol = blanks >> skipWhile isEOL