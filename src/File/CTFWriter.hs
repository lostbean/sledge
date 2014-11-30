{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TypeSynonymInstances #-}
{-# LANGUAGE FlexibleInstances #-}

module File.CTFWriter
       ( renderCTFFile
       ) where

import qualified Data.Vector          as V
import qualified Data.ByteString      as BS
import qualified Data.ByteString.Lazy as BSL
import qualified Data.Text.Lazy       as TL

import           Data.Monoid            ((<>), mconcat, mempty)
import           Data.Text              (Text)
import           Data.Text.Lazy.Builder (toLazyText)
import           Data.Vector            (Vector)

import           Blaze.ByteString.Builder
import           Blaze.ByteString.Builder.Char8
import           Data.Text.Lazy.Builder.RealFloat
import           Data.Text.Lazy.Builder.Int

import           File.CTFReader
import           Texture.Orientation

renderCTFFile :: String -> CTFdata -> IO ()
renderCTFFile fileName ang = let
  txt = toLazyByteString $ renderCTFdata ang
  in BSL.writeFile fileName txt

renderCTFdata :: CTFdata -> Builder
renderCTFdata CTFdata{..} = let
  CTFinfo{..} = ebsdInfo
  ps = map phaseRender phases
  in render ("Channel Text File" :: Text)
  <> setInfo "Prj"     project
  <> setInfo "Author"  author
  <> setInfo "JobMode" jobMode
  <> (gridRender grid <> eol)
  <> renderMany info
  <> mconcat ps
  <> dataHeader
  <> dataRender grid nodes

phaseRender :: CTFphase -> Builder
phaseRender CTFphase{..} = let
  is = map ((tab <>) . render) phaseInfo
  in setInfo "Phases" phaseID
  <> eol
  <> latticeRender latticeCons
  <-> render materialName
  <> mconcat is

dataHeader :: Builder
dataHeader = let
  hs :: [Text]
  hs = [ "Phase","X", "Y", "Bands", "Error", "Euler1"
       , "Euler2", "Euler3", "MAD", "BC", "BS" ]
  in eol <> renderMany hs

latticeRender :: (Double, Double, Double, Double, Double, Double) -> Builder
latticeRender (a, b, c, alpha1, alpha2, alpha3) = let
  sep = fromChar ';'
  lat = renderF 3 a <> sep <> renderF 3 b <> sep <> renderF 3 c
  ang = render alpha1 <> sep <> render alpha2 <> sep <> render alpha3
  in lat <-> ang

gridRender :: CTFgrid -> Builder
gridRender (CTFgrid (row, col) (xstep, ystep) (x0, y0, z0)) =
     setInfo "XCells" row
  <> setInfo "YCells" col
  <> setInfo "XStep"  xstep
  <> setInfo "YStep"  ystep
  <> setInfo "AcqE1"  x0
  <> setInfo "AcqE2"  y0
  <> setInfo "AcqE3"  z0

dataRender :: CTFgrid -> Vector CTFpoint -> Builder
dataRender g = V.foldl' (<>) (fromText "")  . V.map (pointRender g)

orientationRender :: Quaternion -> Builder
orientationRender q = let
  euler = fromQuaternion q
  in     renderF 4 (phi1 euler)
     <-> renderF 4 (phi  euler)
     <-> renderF 4 (phi2 euler)

pointRender :: CTFgrid -> CTFpoint -> Builder
pointRender g CTFpoint{..} = let
  (x, y) = getGridPoint g (xpos, ypos)
  in eol
     <>  render phase
     <-> render x
     <-> render y
     <-> render bands
     <-> render errorID
     <-> orientationRender rotation
     <-> render angulardev
     <-> render bandcontrast
     <-> render bandslope

-- Calculate colunm position (row,col) from ID sequence
-- Origin at (1,1)
getGridPoint :: CTFgrid -> (Int, Int) -> (Double, Double)
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

tab :: Builder
tab = fromChar '\t'

eol :: Builder
eol = render ("\n" :: Text)

renderF :: (RealFloat a)=> Int -> a -> Builder
renderF n = fromLazyText . TL.take 6 . toLazyText . formatRealFloat Fixed (Just n)

(<->) :: (Render a, Render b)=> a -> b -> Builder
a <-> b = render a <> tab <> render b

renderMany :: (Render a)=> [a] -> Builder
renderMany [] = mempty
renderMany (x : xs) = let
  rs = map ((tab <>) . render) xs
  in render x <> mconcat rs

setInfo :: (Render a)=> Text -> a -> Builder
setInfo name value = let
  ntxt = render name
  txt  = render value
  in eol <> ntxt <> tab <> txt
