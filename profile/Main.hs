{-# LANGUAGE RecordWildCards #-}
module Main where

import           Options.Applicative
import           Prelude

import           Texture.SphericalHarmonics
import           TestTexture

data Tester =
  Tester
  { run_test_HSH :: Bool
  } deriving (Show)

tester :: Parser Tester
tester = Tester
  <$> switch
      (  long "test-HSH"
      <> short 's'
      <> help "Run test on HSH." )

main :: IO ()
main = execParser opts >>= run
  where
    opts = info (helper <*> tester)
      ( fullDesc
      <> progDesc "Test and profile sledge library."
      <> header "Hammer library" )

run :: Tester -> IO ()
run Tester{..} = do
  if run_test_HSH
    then do
    testRotSymm
    --runTextureChecker
    else putStrLn "Skiping HSH test."
