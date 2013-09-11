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
       , bingMode
       , scatter

         -- * Main functions
       , mkBingham
       , binghamPDF
       , sampleBingham
       , fitBingham

         -- * Render Bingham distribution
       , renderBingham
       , renderBinghamOmega
       , renderBinghamToEuler

       , testSampleFit
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Control.Applicative ((<$>))
import           Data.Vector         (Vector)
import           System.Random       (randomIO, randoms, newStdGen)

import           Hammer.Math.Algebra
import           Hammer.Math.Optimum
import           Hammer.Texture.Orientation
import           Hammer.Texture.BinghamConstant
import           Hammer.Render.VTK.VTKRender

import           Debug.Trace
dbg s x = trace (s L.++ show x) x

data DFDzx =
  DFDzx
  { dFdz1 :: Double
  , dFdz2 :: Double
  , dFdz3 :: Double
  } deriving (Show)

-- | Variance directions of the Bingham distribution. The direction must
-- normalized and orthonogal to each other.
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
  x4 = (0, lastDir (snd x1) (snd x2) (snd x3))
  zd = L.sortBy (\a b -> fst a `compare` fst b) [x1, x2, x3, x4]

  -- find the highest concentration and set it as the mode and equal zero
  -- therefore all others concentrations are negative
  mode = snd $ last zd
  [(z1, d1), (z2, d2), (z3, d3)] = let
    zmax = fst $ last zd
    in take 3 $ map (\(z, d) -> (z - zmax, d)) zd

  -- add check orthogonality
  dir = DistDir (quaterVec d1) (quaterVec d2) (quaterVec d3)
  z   = Vec3 z1 z2 z3
  -- space of 3-sphere (S3 c R4)
  (f, (dz1, dz2, dz3, _)) = computeAllF z1 z2 z3 0
  dF  = DFDzx dz1 dz2 dz3
  s   = scatterMatrix f dir dF mode
  in Bingham z dir f mode dF s

-- | Eval the Bingham distribution at specific point
binghamPDF :: Bingham -> Quaternion -> Double
binghamPDF Bingham{..} q = let
  v = quaterVec q
  k = concentrations &. (mapVec (\x->x*x) $ directions |> v)
  in exp k / normalization

-- =================================== Statistics ======================================

scatterMatrix :: Double -> DistDir -> DFDzx -> Quaternion -> Mat4
scatterMatrix f DistDir{..} DFDzx{..} qMode = let
  mode = quaterVec qMode
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

-- | Find the 4th direction normalized and orthogonal to all other distbution directions.
lastDir :: Quaternion -> Quaternion -> Quaternion -> Quaternion
lastDir q1 q2 q3 = let
  dir  = DistDir (quaterVec q1) (quaterVec q2) (quaterVec q3)
  -- find the sub matrix 3x3 with the maximum abs determinat
  dets = V.generate 4 (cutDistDir dir)
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
  in mkQuaternion x

-- ======================================== Sampling =====================================

-- | Metroplis-Hastings sampler for the Bingham distribution
sampleBingham :: Bingham -> Int -> IO [Quaternion]
sampleBingham dist@Bingham{..} n = func [] ti wi qi burn (0 :: Int)
  where
    (v, conc) = symmEigen scatter
    burn      = -20
    z         = vecMap sqrt conc

    qi = bingMode
    ti = binghamPDF dist qi        -- target
    wi = centralGaussianPDF z v qi -- proposal

    func acc t0 p0 q0 count accepted = do
      q <- quaternionNormalRand z v q0
      let
        t = binghamPDF dist q        -- target
        p = centralGaussianPDF z v q -- proposal
        a = (p0 * t) / (p * t0)      -- candidate changes of approvel
        out rnd
          -- trace ("> ratio = " ++ show (fromIntegral accepted / fromIntegral count))
          | count > n = return $  acc
          -- burn-in (training samples)
          | count < 0 &&
            a > rnd   = func       acc  t  p  q  (count+1) accepted
          | count < 0 = func       acc  t0 p0 q0 (count+1) accepted

          | a > rnd   = func (q  : acc) t  p  q  (count+1) (accepted+1)
          | otherwise = func (q0 : acc) t0 p0 q0 (count+1)  accepted
      random01 >>= out

-- | Sample a value between [0,1] from a uniform distribution.
random01 :: IO Double
random01 = abs <$> randomIO

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
  foo = mkQuaternion . (&+ quaterVec q)
  in foo <$> multiNormalRand z cv

-- | Compute an angular central gaussian pdf in principal components form for
-- a zero-mean distribution. Uses the eigenvectors and eigenvalues to calculate
-- of the inverse of the covariance matrix.
centralGaussianPDF :: Vec4 -> Mat4 -> Quaternion -> Double
centralGaussianPDF z cv qx = let
  inte = vecFoldr (*) z
  d    = 4
  x    = quaterVec qx
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
  mu = quaterVec qmu
  x  = quaterVec qx
  d  = 4
  dx = x &- mu
  k  = ((2 * pi) ** (d/2)) * inte
  -- Is the same as: SUM ( <dx, V[i]> / z[i] )^2
  e  = normsqr $ (dx .* (transpose cv)) &! (vecMap (1/) z)
  in (exp $ (-0.5) * e) / k

-- ===================================== Bingham Fit =====================================

newtype DFF  = DFF {unDFF :: Vec4} deriving (Show)

-- | Fit a Bingham distribution on a given set of quaternions using MLE (Maximum Likelihhod
-- Estimation) to determine the concentration values.
fitBingham :: Vector Quaternion -> Bingham
fitBingham qs = let
  scatter = calcInertiaMatrix qs

  (ev, ex) = symmEigen scatter
  -- sort decomposition by eigenvalues (highest eigenvalue -> highest concentration value)
  (sortex, sortei) = let
    xs = zip [_1 ex, _2 ex, _3 ex, _4 ex] [1 :: Int ..]
    in unzip $ L.sortBy (\a b -> compare (fst a) (fst b)) xs

  -- v4 x4 will be axis with the highest eigenvalue (the highest concentration or mode)
  ((axis1, axis2, axis3, _), dff) = let
    v                = transpose ev
    [v1, v2, v3, v4] = map getVec sortei
    [x1, x2, x3, x4] = sortex
    getVec i
      | i == 1    = _1 v
      | i == 2    = _2 v
      | i == 3    = _3 v
      | otherwise = _4 v
    in ((v1, v2, v3, v4), DFF $ Vec4 x1 x2 x3 x4)

  -- z4 by definition is equal zero
  (z1, z2, z3) = lookupConcentration dff

  in mkBingham
     -- for positive concentrations
     (z1, mkQuaternion axis1)
     (z2, mkQuaternion axis2)
     (z3, mkQuaternion axis3)

-- | Calculates the inertia matrix that describes the distribution.
calcInertiaMatrix :: Vector Quaternion -> Mat4
calcInertiaMatrix qs
  | n > 0     = total &* (1/n)
  | otherwise = zero
  where
    n     = fromIntegral (V.length qs)
    total = V.foldl' func zero qs
    func acc q = let
      v = quaterVec q
      in acc &+ (outer v v)

-- | Find the concentration values using the gradient descent method.
lookupConcentration :: DFF -> (Double, Double, Double)
lookupConcentration dFF = unVec3 $ bfgs defaultBFGS (errorFunc dFF) zero

evaldFF :: Vec3 -> DFF
evaldFF (Vec3 z1 z2 z3) = let
  (f, (dFdz1, dFdz2, dFdz3, dFdz4)) = computeAllF z1 z2 z3 0
  in DFF $ Vec4 (dFdz1 / f) (dFdz2 / f) (dFdz3 / f) (dFdz4 / f)

-- | Evaluates the error function.
evalError :: DFF -> Vec3 -> Double
evalError dY z = normsqr $ (unDFF $ evaldFF z) &- (unDFF dY)

-- | Calculates the partial derivatives of the error function.
errorFunc :: DFF -> Vec3 -> (Double, Vec3)
errorFunc dY (Vec3 z1 z2 z3) = let
  k = 0.001
  fdz z
    | z == 0    = k * k
    | otherwise = abs $ z * k
  dz1 = fdz z1
  dz2 = fdz z2
  dz3 = fdz z3

  g   = evalError dY (Vec3  z1      z2      z3      )
  g1  = evalError dY (Vec3 (z1+dz1) z2      z3      )
  g2  = evalError dY (Vec3  z1     (z2+dz2) z3      )
  g3  = evalError dY (Vec3  z1      z2      (z3+dz3))

  dg1 = (g1 - g) / dz1
  dg2 = (g2 - g) / dz2
  dg3 = (g3 - g) / dz3
  in (g, (Vec3 dg1 dg2 dg3))

-- ====================================== Plot Space =====================================

sphere :: Int -> Int -> (Vector Vec3, Vector SphereCell)
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

    mkMesh :: Vector (Vector Int) -> Vector SphereCell
    mkMesh vl = V.ifoldl' (\acc i l -> acc V.++ mkStrip l (vl V.! (i+1))) V.empty (V.init vl)
    mkStrip l1 l2 = let
      func i a b = SphereCell (b, l2 V.! (i+1), l1 V.! (i+1), a)
      merge      = SphereCell (V.head l2, V.head l1, V.last l1, V.last l2)
      in merge `V.cons` V.izipWith func (V.init l1) (V.init l2)

-- local instance to avoid conflict when exported.
newtype SphereCell = SphereCell (Int, Int, Int, Int)

instance RenderCell SphereCell where
  makeCell (SphereCell (a, b, c, d)) = U.fromList [a, b, c, d]
  getType _                          = VTK_QUAD

renderSphereVTK :: (Vector Vec3, Vector SphereCell) -> VTK Vec3
renderSphereVTK (ps, quads) = mkUGVTK "hypersphere" (V.convert ps) quads

renderPoints :: [Quaternion] -> VTK Vec3
renderPoints lq = let
  (pos, omega) = V.unzip . V.map (axisAngle . fromQuaternion) $ V.fromList lq
  pids = V.enumFromN (0 :: Int) (V.length pos)
  vtk  = mkUGVTK "samples" (V.convert pos) pids
  attr = mkPointAttr ("Omegas") (\a _ -> (omega) V.! a)
  in addDataPoints vtk attr

-- | Render one sphere with n divisions of a the Bingham distribution at
-- a fixed rotation angle (Omega ~ [0 .. 2*pi]).
renderBinghamOmega :: Bingham -> Int -> Double -> VTK Vec3
renderBinghamOmega dist n omega = let
  mesh@(ps, _) = sphere n n
  vtk          = renderSphereVTK mesh
  attr         = mkPointAttr "intensity" (\a _ -> intensity V.! a)
  axis2q x     = toQuaternion $ mkAxisPair x (Rad omega)
  intensity    = V.map (binghamPDF dist . axis2q) ps
  in addDataPoints vtk attr

-- | Render one sphere with n divisions of a the Bingham distribution.
renderBingham :: Bingham -> Int -> VTK Vec3
renderBingham dist n = let
  step   = pi / fromIntegral n
  vtk    = renderSphereVTK mesh
  omegas = V.generate n ((step *) . fromIntegral)
  attr o = mkPointAttr ("intensity-" ++ show o) (\a _ -> (intensity o) V.! a)
  mesh@(ps, _) = sphere n n
  axis2q o x   = toQuaternion $ mkAxisPair x (Rad o)
  intensity o  = V.map (binghamPDF dist . axis2q o) ps
  in V.foldr (\o acc -> addDataPoints acc (attr o)) vtk omegas
     -- foldl crashs
     -- in V.foldr (\o acc -> addDataPoints acc (attr o)) vtk omegas

-- | Render the PDF of a Bingham distribution in the Euler space.
renderBinghamToEuler :: (Angle a, Num a)=> a -> (Int, Int, Int) -> Bingham -> VTK Double
renderBinghamToEuler step (np1, np, np2) dist = let
  start = toAngle 0
  angs = V.enumFromStepN start step (max np1 (max np np2))
  eu   = [ mkEuler (angs V.! i) (angs V.! j) (angs V.! k) | i <- [0 .. np1-1]
                                                          , j <- [0 .. np -1]
                                                          , k <- [0 .. np2-1]]
  inte = V.map (binghamPDF dist . toQuaternion) (V.fromList eu)
  orig = (fromAngle start, fromAngle start, fromAngle start)
  spc  = (fromAngle step,  fromAngle step,  fromAngle step)
  vtk  = mkSPVTK "Euler_ODF" (np1, np, np2) orig spc
  attr = mkPointAttr "PDF" (\i _ -> inte V.! i)
  in addDataPoints vtk attr

-- ====================================== Test ===========================================

writeQuater :: (RenderElemVTK a)=> String -> VTK a -> IO ()
writeQuater name = writeUniVTKfile ("/home/edgar/Desktop/" ++ name ++ ".vtu") True

testDist :: Bingham
testDist = let
  d1 = (30, mkQuaternion (Vec4 0 0 1 0))
  d2 = (2,  mkQuaternion (Vec4 0 1 0 0))
  d3 = (1,  mkQuaternion (Vec4 1 0 0 0))
  dist = mkBingham d1 d2 d3
  in dist

testSample :: Int -> IO ()
testSample n = let
  d1 = (5, mkQuaternion (Vec4 0 0 1   0 ))
  d2 = (2, mkQuaternion (Vec4 0 0 0 (-1)))
  d3 = (1, mkQuaternion (Vec4 0 1 0   0 ))
  dist = mkBingham d1 d2 d3
  in do
     a <- sampleBingham dist n
     putStrLn $ show dist
     putStrLn $ showPretty $ scatter dist
     putStrLn $ showPretty $ calcInertiaMatrix $ V.fromList a
     writeQuater "Bing-PDF-testSample" $ renderBingham dist 20
     writeQuater "Bing-Samples-testSample" $ renderPoints a
     writeQuater "Euler-PDF-testSample" $ renderBinghamToEuler (Deg 5) (72, 36, 72) dist

-- | Sample and fit the sampled quaternions.
testSampleFit :: Double ->  Double ->  Double -> IO ()
testSampleFit z1 z2 z3 = let
  d1 = (z1, mkQuaternion (Vec4 0 0 1 0))
  d2 = (z2, mkQuaternion (Vec4 0 1 0 0))
  d3 = (z3, mkQuaternion (Vec4 0 0 0 1))
  dist = mkBingham d1 d2 d3
  in do
    a <- sampleBingham dist 100000
    let
      av    = V.fromList a
      dist2 = fitBingham av
    putStrLn ">>> Initial distribution:"
    putStrLn $ show dist
    putStrLn $ showPretty $ scatter dist
    putStrLn "\n>>> Sampled inertia matrix:"
    putStrLn $ showPretty $ calcInertiaMatrix av
    putStrLn "\n>>> Final distribution:"
    putStrLn $ show dist2
    putStrLn $ showPretty $ scatter dist2
    writeQuater "Bing-PDF-Intial" $ renderBingham dist 20
    writeQuater "Bing-PDF-Final" $ renderBingham dist2 20
    writeQuater "Bing-Samples-testSampleFit" $ renderPoints a

-- | Use Monte Carlo integration to check to normalization of the distribution.
testNormalization :: Double -> Double -> Double -> IO Double
testNormalization z1 z2 z3 = let
  d1 = (z1, mkQuaternion (Vec4 0 0 1 0))
  d2 = (z2, mkQuaternion (Vec4 0 1 0 0))
  d3 = (z3, mkQuaternion (Vec4 0 0 0 1))
  dist = mkBingham d1 d2 d3
  in do
    gen <- newStdGen
    let
      n  = 1000000
      qs = take n $ randoms gen
      s  = sum $ map (binghamPDF dist) qs
      v  = surface_area_sphere 3
    return $ (v / (fromIntegral n)) * s
