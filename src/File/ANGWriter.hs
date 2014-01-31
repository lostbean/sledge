{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}

module File.ANGWriter
       ( renderANGFile
       , renderEBSDdata
       ) where

import qualified Data.Vector          as V
import qualified Data.ByteString      as BS
import qualified Data.ByteString.Lazy as BSL

import           Data.Monoid            ((<>), mconcat)
import           Data.Text              (Text)
import           Data.Text.Lazy.Builder (toLazyText)
import           Data.Vector            (Vector)

import           Blaze.ByteString.Builder
import           Blaze.ByteString.Builder.Char8
import           Data.Text.Lazy.Builder.RealFloat
import           Data.Text.Lazy.Builder.Int

import           File.ANGReader
import           Texture.Orientation

renderANGFile :: String -> EBSDdata -> IO ()
renderANGFile fileName ang = let
  txt = toLazyByteString $ renderEBSDdata ang
  in BSL.writeFile fileName txt

renderEBSDdata :: EBSDdata -> Builder
renderEBSDdata EBSDdata{..} = let
  comm = render ("\n#" :: Text)
  in headRender ebsdInfo
  <> comm
  <> (mconcat $ map ((comm <>) . phaseRender) phases)
  <> comm
  <> gridRender grid
  <> (mconcat $ replicate 6 comm)
  <> dataRender grid nodes

headRender :: EBSDinfo -> Builder
headRender EBSDinfo{..} =
  render ("# TEM_PIXperUM" :: Text)
  <> setInfo "# x-star"          xstart
  <> setInfo "# y-star"          ystart
  <> setInfo "# z-star"          zstart
  <> setInfo "# WorkingDistance" workDist

phaseRender :: EBSDphase -> Builder
phaseRender EBSDphase{..} =
     setInfo "# Phase"          phase
  <> setInfo "# MaterialName"   materialName
  <> setInfo "# Formula"        formula
  <> setInfo "# Info"           info
  <> setInfo "# Symmetry"       symmetry
  <> latticeRender              latticeCons
  <> setInfo "# NumberFamilies" numFamilies
  <> hklRender hklFamilies
  <> (mconcat $ map elasticRender elasticConst)
  <> categoryRender categories

latticeRender :: (Double, Double, Double, Double, Double, Double) -> Builder
latticeRender (a, b, c, alpha1, alpha2, alpha3) = let
  txt = a <-> b <-> c <-> alpha1 <-> alpha2 <-> alpha3
  in setInfo "# LatticeConstants" txt

elasticRender :: (Double, Double, Double, Double, Double, Double) -> Builder
elasticRender (i1,i2,i3,i4,i5,i6) = let
  txt = i1 <-> i2 <-> i3 <-> i4 <-> i5 <-> i6
  in setInfo "# ElasticConstants" txt

hklRender :: [(Int, Int, Int, Int, Double, Int)] -> Builder
hklRender hkl = let
  foo (h, k, l, a, b, c) = let
    txt = h <-> k <-> l <-> a <-> b <-> c
    in setInfo "# hklFamilies" txt
  in mconcat $ map foo hkl

categoryRender :: (Int,Int,Int,Int,Int) -> Builder
categoryRender (a,b,c,d,f) = let
  txt = a <-> b <-> c <-> d <-> f
  in setInfo "# Categories" txt

gridRender :: Gridinfo -> Builder
gridRender (Gridinfo (row, ceven, codd) (xstep, ystep) isHex) =
     setInfo "# GRID:"       (gridType isHex)
  <> setInfo "# XSTEP:"      xstep
  <> setInfo "# YSTEP:"      ystep
  <> setInfo "# NCOLS_ODD:"  codd
  <> setInfo "# NCOLS_EVEN:" ceven
  <> setInfo "# NROWS:"      row

dataRender :: Gridinfo -> Vector EBSDpoint -> Builder
dataRender g = V.foldl' (<>) (fromText "")  . V.map (pointRender g)

orientationRender :: Quaternion -> Builder
orientationRender q = let
  euler = fromQuaternion q
  in     renderF 5 (phi1 euler)
     <-> renderF 5 (phi  euler)
     <-> renderF 5 (phi2 euler)

pointRender :: Gridinfo -> EBSDpoint -> Builder
pointRender g EBSDpoint{..} = let
  (x, y) = getGridPoint g (xpos, ypos)
  in render ("\n  " :: Text)
     <>  orientationRender rotation
     <-> renderF 5 x
     <-> renderF 5 y
     <-> renderF 1 qi
     <-> renderF 3 ci
     <-> phaseNum
     <-> detecInt

-- Calculate colunm position (row,col) from ID sequence
-- Origin at (1,1)
getGridPoint :: Gridinfo -> (Int, Int) -> (Double, Double)
getGridPoint g (xi, yi) = let
  (xstep, ystep) = xystep g
  in (fromIntegral xi * xstep, fromIntegral yi * ystep)

class Render a where
  render :: a -> Builder

instance Render Builder where
  render = id

instance Render BS.ByteString where
  render = fromByteString

instance Render Text where
  render = fromText

instance Render String where
  render = fromString

instance Render Int where
  render = fromLazyText . toLazyText . decimal

instance Render Double where
  render = renderF 4

renderF :: (RealFloat a)=> Int -> a -> Builder
renderF n = fromLazyText . toLazyText . formatRealFloat Fixed (Just n)

(<->) :: (Render a, Render b)=> a -> b -> Builder
a <-> b = render a <> (fromChar ' ') <> render b

setInfo :: (Render a)=> Text -> a -> Builder
setInfo name value = let
  ntxt = render name
  txt  = render value
  in fromChar '\n' <> ntxt <> fromChar '\t' <> txt

gridType :: Bool -> Text
gridType g
  | g         = "HexGrid"
  | otherwise = "SqrGrid"
