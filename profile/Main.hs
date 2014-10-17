{-# LANGUAGE RecordWildCards #-}
module Main where

import qualified Data.Vector.Unboxed as U

import Options.Applicative
import Control.Monad

import Hammer.VTK
import Hammer.Math.Algebra

import Texture.Bingham
import Texture.Orientation
import Texture.HyperSphere
import Texture.SphericalHarmonics
import Texture.Sampler
--import TestTexture
import TestKernel


data Tester =
  Tester
  { run_rot_SH  :: Bool
  , run_rot_HSH :: Bool
  , run_sym_SH  :: Bool
  , run_sym_HSH :: Bool
  , run_fam_SH  :: Bool
  , run_fam_HSH :: Bool
  , run_fit_HSH :: Bool
  , run_ker_est :: Bool
  } deriving (Show)

tester :: Parser Tester
tester = Tester
  <$> switch
      (  long "rot-SH"
      <> help "Run test on SH rotation." )
  <*> switch
      (  long "rot-HSH"
      <> help "Run test on HSH rotation." )
  <*> switch
      (  long "symm-SH"
      <> help "Run test on SH symmetry." )
  <*> switch
      (  long "symm-HSH"
      <> help "Run test on HSH symmetry." )
  <*> switch
      (  long "base-SH"
      <> help "Plot SH base functions." )
  <*> switch
      (  long "base-HSH"
      <> help "Plot HSH base functions." )
  <*> switch
      (  long "fit-HSH"
      <> help "fit HSH on 1000 points" )
  <*> switch
      (  long "ker-est"
      <> help "test kernel estimation" )

main :: IO ()
main = execParser opts >>= run
  where
    opts = info (helper <*> tester)
      ( fullDesc
      <> progDesc "Test and profile sledge library."
      <> header "Hammer library" )

run :: Tester -> IO ()
run Tester{..} = do
  when run_rot_SH  testRotSH
  when run_rot_HSH testRotHSH
  --when run_sym_SH  testSymmSH  -- no implemented yet
  when run_sym_HSH testSymmHSH
  when run_fam_SH  (plotSHFuncFamily 10)
  when run_fam_HSH (plotHSHFuncFamily 10)
  when run_fit_HSH (testSamplerAndFit 1000)
  when run_ker_est (testKernel)

testSamplerAndFit :: Int -> IO ()
testSamplerAndFit n = let
  da1 = (1, mkQuaternion  (Vec4 0 0 1 0))
  da2 = (1, mkQuaternion  (Vec4 0 1 0 0))
  da3 = (20, mkQuaternion (Vec4 1 0 0 1))
  da  = mkBingham da1 da2 da3
  db1 = (1, mkQuaternion  (Vec4 1 0 0 (-1)))
  db2 = (20, mkQuaternion (Vec4 0 1 0 0))
  db3 = (1, mkQuaternion  (Vec4 1 0 0 1))
  db  = mkBingham db1 db2 db3
  in do
    xs <- hitAndRunSlice (\q -> binghamPDF da q + binghamPDF db q) (zerorot) (2*pi) n
    writeQuater "Bing-PDF-A-testSamplerMultiModal" $ renderBingham da
    writeQuater "Bing-PDF-B-testSamplerMultiModal" $ renderBingham db
    let ss = map quaternionToSO3 xs
    writeQuater "Bing-Samples-testSamplerMultiModal" $ renderSO3PointsVTK $ U.fromList ss
    writeQuater "HSH-Samples-testSamplerMultiModal" $ plotHSH 20 ss id

writeQuater :: (RenderElemVTK a)=> String -> VTK a -> IO ()
writeQuater name = writeUniVTKfile (name ++ ".vtu") True
