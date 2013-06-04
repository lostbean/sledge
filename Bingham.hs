{-# LANGUAGE NamedFieldPuns             #-}
{-# LANGUAGE RecordWildCards            #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE BangPatterns               #-}
{-# LANGUAGE FlexibleInstances          #-}
{-# LANGUAGE MultiParamTypeClasses      #-}

module Hammer.Texture.Bingham
       ( Bingham
       , concentrations
       , directions
       , normalization
       , mkBingham
       , binghamPDF

       , computeF
       , renderBingham
       , renderFull
       , writeQuater
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Control.Applicative ((<$>))
import           Data.Vector         (Vector)
import           System.Random       (randomIO)
  
import           Math.Gamma

import           Hammer.Math.Algebra
import           Hammer.Texture.Orientation
import           Hammer.Texture.BinghamNormalization
import           Hammer.Render.VTK.VTKRender
  
--import           Debug.Trace
--dbg s x = trace (s L.++ show x) x

epsilon :: Double
epsilon = 1e-8

data DFDzx =
  DFDzx
  { dFdz1 :: Double
  , dFdz2 :: Double
  , dFdz3 :: Double
  } deriving (Show)

data DistDir =
  DistDir
  { v1 :: Vec4
  , v2 :: Vec4
  , v3 :: Vec4
  } deriving (Show)

data Bingham =
  Bingham
  { concentrations :: Vec3        -- ^ concentration values for each directions where
                                  -- z1 <= z2 <= z3 <= (z4 = 0) 
  , directions     :: DistDir     -- ^ unitary vectors Vx in rows
  , normalization  :: Double      -- ^ normalization factor F
  , bingMode       :: Quaternion  -- ^ mode of the distribution
  , partials       :: DFDzx       -- ^ partials derivates of F in each direction
  , scatter        :: Mat4        -- ^ scatter matrix
  } deriving (Show)

(|>) :: DistDir -> Vec4 -> Vec3
DistDir{..} |> v = Vec3 (v1 &. v) (v2 &. v) (v3 &. v)

mkBingham :: (Double, Quaternion) -> (Double, Quaternion) -> (Double, Quaternion) -> Bingham
mkBingham  x1 x2 x3 = let
  [(z1, d1), (z2, d2), (z3, d3)] = L.sortBy (\a b -> fst a `compare` fst b) [x1, x2, x3]
  -- add check orthogonality
  dir = DistDir (quater2vec d1) (quater2vec d2) (quater2vec d3)
  z   = Vec3 z1 z2 z3
  -- space of 3-sphere (S3 c R4)
  f   = computeF 3 z1 z2 z3
  dz1 = computedFdz1 3 z1 z2 z3
  dz2 = computedFdz2 3 z1 z2 z3
  dz3 = computedFdz3 3 z1 z2 z3
  dF  = DFDzx dz1 dz2 dz3
  s   = scatterMatrix f dir z dF mode
  mode = binghamMode dir
  in Bingham z dir f mode dF s

quater2vec :: Quaternion -> Vec4
quater2vec (Quaternion (q0, qv)) = extendHeadWith q0 qv

vec2quater :: Vec4 -> Quaternion
vec2quater (Vec4 q0 q1 q2 q3) = Quaternion (q0, Vec3 q1 q2 q3)
  
-- | Eval the Bingham distribution at specific point
binghamPDF :: Bingham -> Quaternion -> Double
binghamPDF Bingham{..} q = let
  v = quater2vec q
  k = concentrations &. (mapVec (\x->x*x) $ directions |> v)
  in exp k / normalization
     
-- =================================== Statistics ======================================

scatterMatrix :: Double -> DistDir -> Vec3 -> DFDzx -> Quaternion -> Mat4
scatterMatrix f DistDir{..} z DFDzx{..} qMode = let
  mode = quater2vec qMode
  km = 1 - (dFdz1 + dFdz2 + dFdz3) / f
  sm = km *& outer mode mode
  s1 = (dFdz1 / f) *& outer v1 v1
  s2 = (dFdz2 / f) *& outer v2 v2
  s3 = (dFdz3 / f) *& outer v3 v3
  in sm &+ s1 &+ s2 &+ s3

-- | Break the distribution directions (3x4 matrix) in
-- a 3x3 matrix by extracting the colunm n.
cutDistDir :: DistDir -> Int -> (Mat3, Vec3)
cutDistDir DistDir{..} n
  | n == 0    = (mkMat cut0, mkVec _1)
  | n == 1    = (mkMat cut1, mkVec _2)
  | n == 2    = (mkMat cut2, mkVec _3)
  | otherwise = (mkMat cut3, mkVec _4)
  where
    mkMat f = Mat3 (f v1) (f v2) (f v3) 
    mkVec f = Vec3 (f v1) (f v2) (f v3) 
    cut0 (Vec4 _ b c d) = Vec3 b c d
    cut1 (Vec4 a _ c d) = Vec3 a c d
    cut2 (Vec4 a b _ d) = Vec3 a b d
    cut3 (Vec4 a b c _) = Vec3 a b c 

-- | Find the mode of a bingham distribution. The normalized mode is
-- the quaternion orthogonal to all distbution directions @DistDir@.
binghamMode :: DistDir -> Quaternion
binghamMode directions = let
  -- find the sub matrix 3x3 with the maximum abs determinat
  dets = V.generate 4 (cutDistDir directions)
  cut  = V.maxIndex $ V.map (abs . det . fst) dets
  (m, c) = dets V.! cut
  -- solve the linear system A x = b
  solve a b = (inverse a) *. b
  Vec3 xa xb xc = solve m (neg c)
  -- re-insert the cut colunm with 1
  x | cut == 0  = Vec4 1 xa xb xc
    | cut == 1  = Vec4 xa 1 xb xc
    | cut == 2  = Vec4 xa xb 1 xc
    | otherwise = Vec4 xa xb xc 1
  Vec4 q0 q1 q2 q3 = normalize x
  in Quaternion (q0, Vec3 q1 q2 q3)

-- =========================================== Sampling ===========================================

-- | Metroplis-Hastings sampler for the Bingham distribution
sampleBingham :: Bingham -> Int -> IO [Quaternion]
sampleBingham dist@Bingham{..} n = func [] ti pi qi burn 0
  where
    burn = -20
    (v, conc) = symmEigen scatter
    z  = vecMap sqrt conc
  
    qi = bingMode
    ti = binghamPDF dist qi        -- target
    pi = centralGaussianPDF z v qi -- proposal

    func acc t0 p0 q0 count accepted = do
      q <- quaternionNormalRand z v q0
      let
        t = binghamPDF dist q        -- target
        p = centralGaussianPDF z v q -- proposal
        a = (p0 * t) / (p * t0)      -- candidate changes of approvel
        out rnd
          -- trace ("> ratio = " ++ show (fromIntegral accepted / fromIntegral count))
          | count > n = return $  acc
                        
          | count < 0 &&
            a > rnd   = func       acc  t  p  q (count+1)  accepted    -- burn-in (training samples)
          | count < 0 = func       acc  t0 p0 q0 (count+1) accepted    -- burn-in (training samples)
                        
          | a > rnd   = func (q  : acc) t  p  q  (count+1) (accepted+1)
          | otherwise = func (q0 : acc) t0 p0 q0 (count+1)  accepted
      random01 >>= out
  
-- | Sample a value between [0,1] from a uniform distribution.
random01 :: IO Double
random01 = abs <$> randomIO
           
measureMean :: [Quaternion] -> Vec4
measureMean [] = zero
measureMean l = let
  w = L.foldl' (\acc x -> let v = quater2vec x in acc &+ v) zero l
  n = fromIntegral $ length l
  in w &* (1 / n)

measureCovariance :: [Quaternion] -> Mat4
measureCovariance [] = zero
measureCovariance l = let
  --m = measureMean l
  w = L.foldl' (\acc x -> let v = quater2vec x in acc &+ (outer v v)) zero l
  n = fromIntegral $ length l
  in w &* (1 / n)

-- | Inverse of error function.
invERF :: Double -> Double 
invERF x
  | x < 0     = -invERF (-x)
  | otherwise = let 
    a = 0.147
    y1 = (2 / (pi * a) + log (1 - x * x) / 2)
    y2 = sqrt (y1 * y1 - ( 1 / a) * log (1 - x * x))
    in sqrt (y2 - y1)

-- | Compute the normal PDF.
normalPDF :: Double -> Double -> Double -> Double
normalPDF x mu sigma = let
  dx = x - mu
  s2 = 2 * sigma * sigma
  in exp(-(dx * dx) / s2) / (sqrt (2 * pi) * sigma)

-- | Generate a random sample from a univariate normal distribution.
-- The random input must range [0, 1]
normalSample :: Double -> Double -> Double -> Double
normalSample mu sigma rnd = mu + sigma * sqrt 2 * invERF (2 * rnd - 1)

-- | Sample from a multivariate normal in principal components form for
-- a zero-mean distribution. Uses the eigenvectors and eigenvalues to calculate
-- the inverse of the covariance matrix.
multiNormalRand :: Vec4 -> Mat4 -> IO Vec4
multiNormalRand z cv = randomIO >>= sample
  where
    -- randomIO gives [-1,1] therefore needs abs before call normalSample
    fVecRnd = vecZipWith (\s r -> normalSample 0 s (abs r)) z
    sample  = return . (.* (transpose cv)) . fVecRnd

-- | Sample from an angular central gaussian zero-mean distribution and add
-- it to a given quaternion. Used to create innovation in the Metropolis Hasting
-- sampler.
quaternionNormalRand :: Vec4 -> Mat4 -> Quaternion -> IO Quaternion
quaternionNormalRand z cv q = let
  foo = vec2quater .  normalize . (&+ quater2vec q)
  in foo <$> multiNormalRand z cv
  
-- | Compute an angular central gaussian pdf in principal components form for
-- a zero-mean distribution. Uses the eigenvectors and eigenvalues to calculate
-- of the inverse of the covariance matrix. 
centralGaussianPDF :: Vec4 -> Mat4 -> Quaternion -> Double
centralGaussianPDF z cv qx = let
  inte = vecFoldr (*) z
  d    = 4
  x    = quater2vec qx
  k    = inte * surface_area_sphere (d - 1)
  -- Is the same as: SUM ( <dx, V[i]> / z[i] )^2
  e    = normsqr $ (x .* (transpose cv)) &! (vecMap (1/) z)
  in 1 / (k * e ** (fromIntegral d / 2))

-- | Compute a multivariate normal pdf in principal components form.
-- Uses the eigenvectors and eigenvalues to calculate the inverse of
-- the covariance matrix. 
multiNormalPDF :: Vec4 -> Mat4 -> Quaternion -> Quaternion -> Double
multiNormalPDF z cv qmu qx = let
  inte = vecFoldr (*) z
  mu = quater2vec qmu
  x  = quater2vec qx
  d  = 4
  dx = x &- mu
  k  = ((2 * pi) ** (d/2)) * inte
  -- Is the same as: SUM ( <dx, V[i]> / z[i] )^2
  e  = normsqr $ (dx .* (transpose cv)) &! (vecMap (1/) z)
  in (exp $ (-0.5) * e) / k 

-- ====================================== Plot Space ===========================================

sphere :: Int -> Int -> (Vector Vec3, Vector (Int, Int, Int, Int)) 
sphere nPhi nTheta
  | nPhi   < 3 = sphere 3 nTheta
  | nTheta < 3 = sphere nPhi 3
  | otherwise  = (V.map toCart ql, mkMesh ll)
  where
    -- phi   -> azimuthal angle
    -- theta -> polar     angle
    step_phi   = 2 * pi / (fromIntegral nPhi)
    step_theta = pi / (fromIntegral nTheta)
    funcPHI   = (step_phi *) . fromIntegral
    funcTHETA = (step_theta *) . fromIntegral
    qlm = V.fromList [ (funcPHI np, funcTHETA nt) | nt <- [1..nTheta-1], np <- [0..nPhi-1] ]
    q0  = (0, 0)
    qn  = (2 * pi, pi)
    ql  = q0 `V.cons` qlm `V.snoc` qn
    toCart (phi, theta) = let
      x = sin theta * cos phi
      y = sin theta * sin phi
      z = cos theta
      in Vec3 x y z

    ll0 = V.replicate nPhi 0
    lln = V.replicate nPhi (nPhi * (nTheta-1) + 1)
    llm = V.fromList [ V.enumFromN ((n-1) * nPhi + 1) nPhi | n <- [1..nTheta-1] ]
    ll  = ll0 `V.cons` llm `V.snoc` lln

    mkMesh :: Vector (Vector Int) -> Vector (Int, Int, Int, Int)
    mkMesh vl = V.ifoldl' (\acc i l -> acc V.++ mkStrip l (vl V.! (i+1))) V.empty (V.init vl) 
    mkStrip l1 l2 = let
      func i a b = (b, l2 V.! (i+1), l1 V.! (i+1), a)
      merge = (V.head l2, V.head l1, V.last l1, V.last l2)
      in merge `V.cons` V.izipWith func (V.init l1) (V.init l2)

instance RenderCell (Int, Int, Int, Int) where
  makeCell (a,b,c,d) = U.fromList [a,b,c,d]
  getType _          = VTK_QUAD

axisPair2q :: Double -> Vec3 -> Quaternion
axisPair2q omega n
  | l == 0    = Quaternion (c, s *& (normalize n))
  | otherwise = Quaternion (c, s *& (normalize n))
  where
    s = sin (0.5*omega)
    c = cos (0.5*omega)
    l = norm n
    
q2axisPair :: Quaternion -> (Double, Vec3)
q2axisPair (Quaternion (q0, v1)) = (omega, s *& v1)
  where
    omega = 2 * acos q0
    s = 1 / sin (0.5 * omega)

renderBingham :: Bingham -> Int -> Double -> VTK Vec3
renderBingham dist n omega = let
  mesh@(ps, _) = sphere n n
  vtk          = renderSphereVTK mesh
  attr         = mkPointAttr "intensity" (\a _ -> intensity V.! a)
  intensity    = V.map (binghamPDF dist . axisPair2q omega) ps
  in addDataPoints vtk attr

renderSphereVTK :: (Vector Vec3, Vector (Int, Int, Int, Int)) -> VTK Vec3
renderSphereVTK (ps, quads) = mkUGVTK "hypersphere" (V.convert ps) quads

renderFull :: Bingham -> Int -> VTK Vec3
renderFull dist n = let
  step   = pi / fromIntegral n
  vtk    = renderSphereVTK mesh
  omegas = V.generate n ((step *) . fromIntegral)
  attr o = mkPointAttr ("intensity-" ++ show o) (\a _ -> (intensity o) V.! a)
  mesh@(ps, _) = sphere n n
  intensity o  = V.map (binghamPDF dist . axisPair2q o) ps
  in V.foldr (\o acc -> addDataPoints acc (attr o)) vtk omegas
     -- foldl crashs
     -- in V.foldr (\o acc -> addDataPoints acc (attr o)) vtk omegas

renderPoints :: [Quaternion] -> VTK Vec3
renderPoints lq = let
  (omega, pos) = V.unzip . V.map q2axisPair $ V.fromList lq
  pids = V.enumFromN (0 :: Int) (V.length pos)
  vtk  = mkUGVTK "samples" (V.convert pos) pids
  attr = mkPointAttr ("Omegas") (\a _ -> (omega) V.! a)
  in addDataPoints vtk attr

-- ====================================== Test ===========================================

writeQuater :: String -> VTK Vec3 -> IO ()
writeQuater name = writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vtu")

testSample :: Int -> IO ()
testSample n = let
  d1 = (10, Quaternion (0, Vec3 0 1 0))
  d2 = (60, Quaternion (0, Vec3 1 0 0))
  d3 = (60, Quaternion (0, Vec3 0 0 1))
  dist = mkBingham d1 d2 d3
  in do
     a <- sampleBingham dist n
     putStrLn $ show dist
     putStrLn $ showPretty $ scatter dist
     putStrLn $ showPretty $ measureCovariance a
     writeQuater "BingPDF" $ renderFull dist 20
     writeQuater "BingSamples" $ renderPoints a

testDist :: Bingham
testDist = let
  d1 = (30, Quaternion (0, Vec3 0 1 0))
  d2 = (2, Quaternion (0, Vec3 1 0 0))
  d3 = (1, Quaternion (1, Vec3 0 0 0))
  dist = mkBingham d1 d2 d3
  in dist
