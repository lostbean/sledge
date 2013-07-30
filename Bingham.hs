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
       , binghamMode
       , sampleBingham
       , fitBingham

         -- * Render Bingham distribution
       , renderBingham
       , renderBinghamOmega
       , renderBinghamToEuler
       , testFit
       , writeAllTables
       ) where

import qualified Data.List           as L
import qualified Data.Vector         as V
import qualified Data.Vector.Unboxed as U

import           Control.Applicative ((<$>))
import           Data.Vector         (Vector)
import           System.Random       (randomIO)

import           Hammer.Math.Algebra
import           Hammer.Texture.Orientation
import           Hammer.Texture.BinghamNormalization
import           Hammer.Texture.BinghamTable
import           Hammer.Render.VTK.VTKRender

import           Debug.Trace
dbg s x = trace (s L.++ show x) x

epsilon :: Double
epsilon = 1e-8

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
  [(z1, d1), (z2, d2), (z3, d3)] = L.sortBy (\a b -> fst a `compare` fst b) [x1, x2, x3]
  -- add check orthogonality
  dir = DistDir (quaterVec d1) (quaterVec d2) (quaterVec d3)
  z   = Vec3 z1 z2 z3
  -- space of 3-sphere (S3 c R4)
  f   = computeF 3 z1 z2 z3
  dz1 = computedFdz1 3 z1 z2 z3
  dz2 = computedFdz2 3 z1 z2 z3
  dz3 = computedFdz3 3 z1 z2 z3
  dF  = DFDzx dz1 dz2 dz3
  s   = scatterMatrix f dir dF mode
  mode = binghamMode dir
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

fitBingham :: Vector Quaternion -> Bingham
fitBingham qs = let
  scatter = calcInertiaMatrix qs

  v  = transpose . fst $ symmEigen scatter
  v1 = _2 v
  v2 = _3 v
  v3 = _4 v

  t1 = (v1 .* scatter) &. v1
  t2 = (v2 .* scatter) &. v2
  t3 = (v3 .* scatter) &. v3

  (z1, z2, z3) = lookupConcentration (Vec3 t1 t2 t3)

  in mkBingham
     (z1, mkQuaternion v1)
     (z2, mkQuaternion v2)
     (z3, mkQuaternion v3)

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

lookupConcentration :: Vec3 -> (Double, Double, Double)
lookupConcentration dY = let
  (iz1, iz2, iz3) = fst $ findNearestNaive dY
  initz1 = tableIndex U.! iz1
  initz2 = tableIndex U.! iz2
  initz3 = tableIndex U.! iz3
  initDelta = 10
  numIter   = 50 :: Int
  go 0 step = step
  go n step = case stepDown dY step of
    Just next -> go (n-1) next
    _         -> step
  in fst $ go numIter ((initz1, initz2, initz3), initDelta)

evalLookupError :: Vec3 -> (Double, Double, Double) -> Double
evalLookupError dY (z1, z2, z3) = let
  f     = interpolateF3D     z1 z2 z3
  dFdz1 = interpolatedFdZ13D z1 z2 z3
  dFdz2 = interpolatedFdZ13D z1 z2 z3
  dFdz3 = interpolatedFdZ13D z1 z2 z3
  dF    = Vec3 dFdz1 dFdz2 dFdz3
  in dbg ("err> " ++ show (z1, z2, z3)) $ normsqr $ dF &* (1/f) &- dY

evalErrGradient :: Vec3 -> (Double, Double, Double) -> (Double, Double, Double)
evalErrGradient dY (z1, z2, z3) = let
  dz = 0.01
  g  = evalLookupError dY (z1   , z2   , z3   )
  g1 = evalLookupError dY (z1+dz, z2   , z3   )
  g2 = evalLookupError dY (z1   , z2+dz, z3   )
  g3 = evalLookupError dY (z1   , z2   , z3+dz)
  dgdz1 = (g1 - g) / dz
  dgdz2 = (g2 - g) / dz
  dgdz3 = (g3 - g) / dz
  in (dgdz1, dgdz2, dgdz3)

type Step = ((Double, Double, Double), Double)

stepDown :: Vec3 -> Step -> Maybe Step
stepDown dY (con@(z1, z2, z3), delta) = next
  where
    refErr = evalLookupError dY (z1, z2, z3)
    scales = V.fromList [2.5, 1.5, 1.0, 0.5, 0.25]
    tries  = V.map foo scales
    next   = V.find ((refErr >) . evalLookupError dY . fst) tries
    (d1, d2, d3) = evalErrGradient dY (z1, z2, z3)
    foo scale = let
      k = scale * delta
      nz1 = z1 - k * d1
      nz2 = z2 - k * d2
      nz3 = z3 - k * d3
      in ((nz1, nz2, nz3), k)

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

testSample :: Int -> IO ()
testSample n = let
  d1 = (1, mkQuaternion (Vec4 0 0 1 0))
  d2 = (3, mkQuaternion (Vec4 0 1 0 0))
  d3 = (6, mkQuaternion (Vec4 0 0 0 1))
  dist = mkBingham d1 d2 d3
  in do
     a <- sampleBingham dist n
     putStrLn $ show dist
     putStrLn $ showPretty $ scatter dist
     putStrLn $ showPretty $ calcInertiaMatrix $ V.fromList a
     writeQuater "BingPDF" $ renderBingham dist 20
     writeQuater "BingSamples" $ renderPoints a
     writeQuater "EulerPDF" $ renderBinghamToEuler (Deg 5) (72, 36, 72) dist

testDist :: Bingham
testDist = let
  d1 = (30, mkQuaternion (Vec4 0 0 1 0))
  d2 = (2,  mkQuaternion (Vec4 0 1 0 0))
  d3 = (1,  mkQuaternion (Vec4 1 0 0 0))
  dist = mkBingham d1 d2 d3
  in dist

testFit :: Bingham
testFit = fitBingham qs
  where
    qs = V.fromList [
      mkQuaternion $ Vec4 0.8776         0         0   (-0.4794),
      mkQuaternion $ Vec4 0.8752         0         0   (-0.4838),
      mkQuaternion $ Vec4 0.8727         0         0   (-0.4882),
      mkQuaternion $ Vec4 0.8703         0         0   (-0.4925),
      mkQuaternion $ Vec4 0.8678         0         0   (-0.4969),
      mkQuaternion $ Vec4 0.8653         0         0   (-0.5012),
      mkQuaternion $ Vec4 0.8628         0         0   (-0.5055),
      mkQuaternion $ Vec4 0.8603         0         0   (-0.5098),
      mkQuaternion $ Vec4 0.8577         0         0   (-0.5141),
      mkQuaternion $ Vec4 0.8551         0         0   (-0.5184),
      mkQuaternion $ Vec4 0.8525         0         0   (-0.5227),
      mkQuaternion $ Vec4 0.8499         0         0   (-0.5269),
      mkQuaternion $ Vec4 0.8473         0         0   (-0.5312),
      mkQuaternion $ Vec4 0.8446         0         0   (-0.5354),
      mkQuaternion $ Vec4 0.8419         0         0   (-0.5396),
      mkQuaternion $ Vec4 0.8392         0         0   (-0.5438),
      mkQuaternion $ Vec4 0.8365         0         0   (-0.5480),
      mkQuaternion $ Vec4 0.8337         0         0   (-0.5522),
      mkQuaternion $ Vec4 0.8309         0         0   (-0.5564),
      mkQuaternion $ Vec4 0.8281         0         0   (-0.5605),
      mkQuaternion $ Vec4 0.8253         0         0   (-0.5646),
      mkQuaternion $ Vec4 0.8225         0         0   (-0.5688),
      mkQuaternion $ Vec4 0.8196         0         0   (-0.5729),
      mkQuaternion $ Vec4 0.8168         0         0   (-0.5770),
      mkQuaternion $ Vec4 0.8139         0         0   (-0.5810),
      mkQuaternion $ Vec4 0.8110         0         0   (-0.5851),
      mkQuaternion $ Vec4 0.8080         0         0   (-0.5891),
      mkQuaternion $ Vec4 0.8051         0         0   (-0.5932),
      mkQuaternion $ Vec4 0.8021         0         0   (-0.5972),
      mkQuaternion $ Vec4 0.7991         0         0   (-0.6012),
      mkQuaternion $ Vec4 0.7961         0         0   (-0.6052),
      mkQuaternion $ Vec4 0.7930         0         0   (-0.6092),
      mkQuaternion $ Vec4 0.7900         0         0   (-0.6131),
      mkQuaternion $ Vec4 0.7869         0         0   (-0.6171),
      mkQuaternion $ Vec4 0.7838         0         0   (-0.6210),
      mkQuaternion $ Vec4 0.7807         0         0   (-0.6249),
      mkQuaternion $ Vec4 0.7776         0         0   (-0.6288),
      mkQuaternion $ Vec4 0.7744         0         0   (-0.6327),
      mkQuaternion $ Vec4 0.7712         0         0   (-0.6365),
      mkQuaternion $ Vec4 0.7681         0         0   (-0.6404),
      mkQuaternion $ Vec4 0.7648         0         0   (-0.6442),
      mkQuaternion $ Vec4 0.7616         0         0   (-0.6480),
      mkQuaternion $ Vec4 0.7584         0         0   (-0.6518),
      mkQuaternion $ Vec4 0.7551         0         0   (-0.6556),
      mkQuaternion $ Vec4 0.7518         0         0   (-0.6594),
      mkQuaternion $ Vec4 0.7485         0         0   (-0.6631),
      mkQuaternion $ Vec4 0.7452         0         0   (-0.6669),
      mkQuaternion $ Vec4 0.7418         0         0   (-0.6706),
      mkQuaternion $ Vec4 0.7385         0         0   (-0.6743),
      mkQuaternion $ Vec4 0.7351         0         0   (-0.6780),
      mkQuaternion $ Vec4 0.7317         0         0   (-0.6816),
      mkQuaternion $ Vec4 0.7283         0         0   (-0.6853),
      mkQuaternion $ Vec4 0.7248         0         0   (-0.6889),
      mkQuaternion $ Vec4 0.7214         0         0   (-0.6925),
      mkQuaternion $ Vec4 0.7179         0         0   (-0.6961),
      mkQuaternion $ Vec4 0.7144         0         0   (-0.6997),
      mkQuaternion $ Vec4 0.7109         0         0   (-0.7033),
      mkQuaternion $ Vec4 0.7074         0         0   (-0.7068),
      mkQuaternion $ Vec4 0.7038         0         0   (-0.7104),
      mkQuaternion $ Vec4 0.7003         0         0   (-0.7139),
      mkQuaternion $ Vec4 0.6967         0         0   (-0.7174),
      mkQuaternion $ Vec4 0.6931         0         0   (-0.7208),
      mkQuaternion $ Vec4 0.6895         0         0   (-0.7243),
      mkQuaternion $ Vec4 0.6859         0         0   (-0.7277),
      mkQuaternion $ Vec4 0.6822         0         0   (-0.7311),
      mkQuaternion $ Vec4 0.6786         0         0   (-0.7345),
      mkQuaternion $ Vec4 0.6749         0         0   (-0.7379),
      mkQuaternion $ Vec4 0.6712         0         0   (-0.7413),
      mkQuaternion $ Vec4 0.6675         0         0   (-0.7446),
      mkQuaternion $ Vec4 0.6637         0         0   (-0.7480),
      mkQuaternion $ Vec4 0.6600         0         0   (-0.7513),
      mkQuaternion $ Vec4 0.6562         0         0   (-0.7546),
      mkQuaternion $ Vec4 0.6524         0         0   (-0.7578),
      mkQuaternion $ Vec4 0.6486         0         0   (-0.7611),
      mkQuaternion $ Vec4 0.6448         0         0   (-0.7643),
      mkQuaternion $ Vec4 0.6410         0         0   (-0.7675),
      mkQuaternion $ Vec4 0.6372         0         0   (-0.7707),
      mkQuaternion $ Vec4 0.6333         0         0   (-0.7739),
      mkQuaternion $ Vec4 0.6294         0         0   (-0.7771),
      mkQuaternion $ Vec4 0.6255         0         0   (-0.7802),
      mkQuaternion $ Vec4 0.6216         0         0   (-0.7833),
      mkQuaternion $ Vec4 0.6177         0         0   (-0.7864),
      mkQuaternion $ Vec4 0.6137         0         0   (-0.7895),
      mkQuaternion $ Vec4 0.6098         0         0   (-0.7926),
      mkQuaternion $ Vec4 0.6058         0         0   (-0.7956),
      mkQuaternion $ Vec4 0.6018         0         0   (-0.7986),
      mkQuaternion $ Vec4 0.5978         0         0   (-0.8016),
      mkQuaternion $ Vec4 0.5938         0         0   (-0.8046),
      mkQuaternion $ Vec4 0.5898         0         0   (-0.8076),
      mkQuaternion $ Vec4 0.5857         0         0   (-0.8105),
      mkQuaternion $ Vec4 0.5817         0         0   (-0.8134),
      mkQuaternion $ Vec4 0.5776         0         0   (-0.8163),
      mkQuaternion $ Vec4 0.5735         0         0   (-0.8192),
      mkQuaternion $ Vec4 0.5694         0         0   (-0.8220),
      mkQuaternion $ Vec4 0.5653         0         0   (-0.8249),
      mkQuaternion $ Vec4 0.5612         0         0   (-0.8277),
      mkQuaternion $ Vec4 0.5570         0         0   (-0.8305),
      mkQuaternion $ Vec4 0.5529         0         0   (-0.8333),
      mkQuaternion $ Vec4 0.5487         0         0   (-0.8360),
      mkQuaternion $ Vec4 0.5445         0         0   (-0.8388),
      mkQuaternion $ Vec4 0.5403         0         0   (-0.8415)
      ]
