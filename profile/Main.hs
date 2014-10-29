{-# LANGUAGE RecordWildCards #-}
module Main where

import Options.Applicative
import Control.Monad

--import TestTexture
import TestKernel
import TestSampler
import TestKernelSampling


data Tester =
  Tester
  { run_ker_est :: Bool
  , run_sap_fit :: Bool
  , run_sap_mul :: Bool
  , run_ker_sap :: Bool
  } deriving (Show)

tester :: Parser Tester
tester = Tester
  <$> switch
      (  long "ker-est"
      <> help "test kernel estimation" )
  <*> switch
      (  long "samp-fit"
      <> help "sample and fit" )
  <*> switch
      (  long "samp-mult"
      <> help "sample multi-modal distribution" )
  <*> switch
      (  long "ker-samp"
      <> help "sampling from a kernel distribution" )

main :: IO ()
main = execParser opts >>= run
  where
    opts = info (helper <*> tester)
      ( fullDesc
      <> progDesc "Test and profile sledge library."
      <> header "Hammer library" )

run :: Tester -> IO ()
run Tester{..} = do
  when run_ker_est (testKernel)
  when run_sap_fit (testSampFit 10000)
  when run_sap_mul (testSampMulti 10000)
  when run_ker_sap (testKernelSampling 1000)
