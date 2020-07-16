{-# LANGUAGE
    RecordWildCards
  , OverloadedStrings
  , TypeSynonymInstances
  , FlexibleInstances
  #-}
module File.ANGWriter
  ( renderANGFile
  , renderANG
  ) where

import Blaze.ByteString.Builder
import Blaze.ByteString.Builder.Char8
import Data.Monoid ((<>))
import Data.Text (Text)
import Data.Text.Lazy.Builder (toLazyText)
import Data.Text.Lazy.Builder.Int
import Data.Text.Lazy.Builder.RealFloat
import Data.Vector (Vector)
import qualified Data.Vector          as V
import qualified Data.ByteString      as BS
import qualified Data.ByteString.Lazy as BSL

import File.ANGReader
import Texture.Orientation

renderANGFile :: String -> ANGdata -> IO ()
renderANGFile fileName = BSL.writeFile fileName . renderANG

renderANG :: ANGdata -> BSL.ByteString
renderANG = toLazyByteString . renderANGdata

renderANGdata :: ANGdata -> Builder
renderANGdata ANGdata{..} =
  headRenderA grid ebsdInfo
  <> mconcat (map ((comment <>) . phaseRender) phases)
  <> comment
  <> gridRender grid
  <> comment
  <> headRenderB ebsdInfo
  <> dataRender grid nodes

headRenderA :: ANGgrid -> ANGinfo -> Builder
headRenderA ANGgrid{..} ANGinfo{..} = let
  (xstart, ystart, zstart) = origin
  in (render ("# TEM_PIXperUM" :: Text) <> tab <> render pixperum)
  <> setInfo "# x-star"          xstart
  <> setInfo "# y-star"          ystart
  <> setInfo "# z-star"          zstart
  <> setInfo "# WorkingDistance" workDist

headRenderB :: ANGinfo -> Builder
headRenderB ANGinfo{..} =
     setInfo "# OPERATOR:" operator
  <> comment
  <> setInfo "# SAMPLEID:" sampleID
  <> comment
  <> setInfo "# SCANID:"   scanID
  <> comment

phaseRender :: ANGphase -> Builder
phaseRender ANGphase{..} =
     setInfo "# Phase"          phase
  <> setInfo "# MaterialName"   materialName
  <> setInfo "# Formula"        formula
  <> setInfo "# Info"           info
  <> setInfo "# Symmetry"       symmetry
  <> latticeRender              latticeCons
  <> setInfo "# NumberFamilies" numFamilies
  <> hklRender hklFamilies
  <> mconcat (map elasticRender elasticConst)
  <> categoryRender categories

latticeRender :: (Double, Double, Double, Double, Double, Double) -> Builder
latticeRender (a, b, c, w1, w2, w3) = let
  txt = renderF 3 a  <-> renderF 3 b  <-> renderF 3 c <->
        renderF 3 w1 <-> renderF 3 w2 <-> renderF 3 w3
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

gridRender :: ANGgrid -> Builder
gridRender (ANGgrid (row, ceven, codd) (xstep, ystep) _ isHex) =
     setInfo "# GRID:"       (gridType isHex)
  <> setInfo "# XSTEP:"      xstep
  <> setInfo "# YSTEP:"      ystep
  <> setInfo "# NCOLS_ODD:"  codd
  <> setInfo "# NCOLS_EVEN:" ceven
  <> setInfo "# NROWS:"      row

dataRender :: ANGgrid -> Vector ANGpoint -> Builder
dataRender g = V.foldl' (<>) (fromText "")  . V.map (pointRender g)

orientationRender :: Quaternion -> Builder
orientationRender q = let
  euler = fromQuaternion q
  toPlus x = if x < 0 then x + 2*pi else x
  in     renderF 5 (toPlus $ phi1 euler)
     <-> renderF 5 (toPlus $ phi  euler)
     <-> renderF 5 (toPlus $ phi2 euler)

pointRender :: ANGgrid -> ANGpoint -> Builder
pointRender g ANGpoint{..} = let
  (x, y) = getGridPoint g (xpos, ypos)
  in render ("\n  " :: Text)
     <>  orientationRender rotation
     <-> renderF 5 x
     <-> renderF 5 y
     <-> renderF 1 iq
     <-> renderF 3 ci
     <-> phaseNum
     <-> detecInt

-- Calculate colunm position (row,col) from ID sequence
-- Origin at (1,1)
getGridPoint :: ANGgrid -> (Int, Int) -> (Double, Double)
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
  render = renderF 6

renderF :: (RealFloat a)=> Int -> a -> Builder
renderF n = fromLazyText . toLazyText . formatRealFloat Fixed (Just n)

(<->) :: (Render a, Render b)=> a -> b -> Builder
a <-> b = render a <> spc <> render b

setInfo :: (Render a)=> Text -> a -> Builder
setInfo name value = let
  ntxt = render name
  txt  = render value
  in eol <> ntxt <> tab <> txt

tab :: Builder
tab = fromChar '\t'

eol :: Builder
eol = fromChar '\n'

spc :: Builder
spc = fromChar ' '

comment :: Builder
comment = render ("\n#" :: Text)

gridType :: Bool -> Text
gridType g
  | g         = "HexGrid"
  | otherwise = "SqrGrid"
