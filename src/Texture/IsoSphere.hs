{-# LANGUAGE
    RecordWildCards
  , BangPatterns
  #-}
module Texture.IsoSphere
  ( IsoSphere (subLevel, vertices, faces, centers)
  , isoSphereZero
  , isoSphere
  , refineIsoSphere
  , scaleIsoSphere
  , angularDistance
  , nearestPoint
  , genIsoSphereSO3Grid
  , getOptSubDivisionLevel
  , renderIsoSphereFaces
  , renderQuaternions
  ) where

import qualified Data.List                   as L
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as U
import qualified Data.Vector.Unboxed.Mutable as UM
import qualified Data.Vector.Mutable         as VM

import Data.Function    (on)
import Control.Monad.ST (runST, ST)
import Data.STRef

import Linear.Vect
import Hammer.VTK

import Texture.Orientation
import Texture.Symmetry

-- | IsoSphere is an geodesic grid on sphere that is based on the subdivision of
-- isodecahedron. It has an uniform and homogeneous distribution of points.
data IsoSphere
  = IsoSphere
  { subLevel :: !Int
  , vertices :: U.Vector Vec3D
  , faces    :: V.Vector (U.Vector (Int, Int, Int))
  , centers  :: V.Vector (U.Vector Vec3D)
  } deriving (Show)

-- | Creates an IsoShpere with given subdivision level.
isoSphere :: Int -> IsoSphere
isoSphere n = U.foldl' (\acc _ -> refineIsoSphere acc) isoSphereZero (U.generate (abs n) id)

-- | Creates an IsoShpere with no subdivision (isodecahedron).
isoSphereZero :: IsoSphere
isoSphereZero = let
  t  = 0.5 * (1 + sqrt 5)
  vs = U.map normalize $ U.fromList
       [ Vec3 (-1) t  0, Vec3 1  t  0, Vec3 (-1) (-t) 0, Vec3 1 (-t)  0
       , Vec3 0 (-1)  t, Vec3 0  1  t, Vec3 0 (-1) (-t), Vec3 0 1  (-t)
       , Vec3 t  0 (-1), Vec3 t  0  1, Vec3 (-t) 0 (-1), Vec3 (-t) 0  1 ]
  fs = U.fromList
       -- 5 faces around point 0
       [ (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11)
       -- 5 adjacent faces
       , (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8)
       -- 5 faces around point 3
       , (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9)
       -- 5 adjacent faces
       , (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1) ]
  in IsoSphere
     { subLevel = 0
     , vertices = vs
     , faces    = V.singleton fs
     , centers  = V.singleton (U.map (getFaceCenter vs) fs)
     }

getFaceCenter :: (DotProd a v, LinearMap a v, Norm a v, UM.Unbox (v a)) => U.Vector (v a) -> (Int, Int, Int) -> v a
getFaceCenter vs (i1, i2, i3) = normalize $ (v1 &+ v3 &+ v2) &* (1/3)
  where
    v1 = vs U.! i1
    v2 = vs U.! i2
    v3 = vs U.! i3

-- | Refine IsoSphere by subdividing each triangular face in four new triangles.
refineIsoSphere :: IsoSphere -> IsoSphere
refineIsoSphere ms@IsoSphere {..} = let
  fs    = if V.null faces then U.empty else V.last faces
  vsize = U.length vertices
  fsize = U.length fs

  -- Each face gives rise to 4 new faces
  newfsize = 4 * fsize

  addF mf poll fid (v1, v2, v3) = do
      let foffset = 4 * fid
      v12 <- addAddGetNewPos poll (v1, v2)
      v23 <- addAddGetNewPos poll (v2, v3)
      v31 <- addAddGetNewPos poll (v3, v1)
      UM.write mf foffset       (v12, v23, v31)
      UM.write mf (foffset + 1) (v31, v1,  v12)
      UM.write mf (foffset + 2) (v12, v2,  v23)
      UM.write mf (foffset + 3) (v23, v3,  v31)

  in runST $ do
    -- new list of faces
    mf <- UM.replicate newfsize (-1,-1,-1)
    -- initialize table of edges
    me <- VM.replicate vsize U.empty
    -- initialize counter to the next position after 'vertices'
    mk <- newSTRef vsize
    -- run subdivision (index only)
    U.zipWithM_ (addF mf (mk, me)) (U.enumFromN 0 fsize) fs
    -- retrieve all mutable data
    ff <- U.unsafeFreeze mf
    ef <- V.unsafeFreeze me
    kf <- readSTRef mk
    -- calculate subdivision points
    ps <- fillVertices ef kf vertices
    -- new list of centers
    let cs = U.map (getFaceCenter ps) ff
    return $ ms
      { subLevel = subLevel + 1
      , faces    = V.snoc faces ff
      , vertices = ps
      , centers  = V.snoc centers cs
      }

-- | Calculate a new points for each existing edge and add it to its correspondent position
-- in a new vertex list.
fillVertices :: V.Vector (U.Vector (Int, Int)) -> Int -> U.Vector Vec3D -> ST s (U.Vector Vec3D)
fillVertices edges vmax points = do
  mv <- UM.new vmax
  let
    func i = do
      -- outer Vector of edges and points must have the same size
      -- copy original point
      UM.write mv i (points U.! i)
      let es = edges V.! i
      -- add new points (one per edge)
      U.mapM_ (\(j, vid) -> UM.write mv vid (getV i j)) es
  mapM_ func [0..n - 1]
  U.unsafeFreeze mv
  where
    n = U.length points
    -- subdivision rule for new points (one per edge)
    getV ia ib = let
      va = points U.! ia
      vb = points U.! ib
      in normalize $ 0.5 *& (va &+ vb)

-- | Verify if an edge (Int, Int) between two vertices was already assigned then retrieves
-- its correspondent position otherwise register the edge to the next available position.
addAddGetNewPos :: (UM.Unbox a, Num a)=> (STRef s a, VM.MVector s (U.Vector (Int, a)))
                -> (Int, Int) -> ST s a
addAddGetNewPos (kref, vs) (i, j) = do
  -- the edges are stored in Vector (Vector (Int, a)) with the lowest value of edge in the
  -- outer Vector and the higher value in the first element of the inner Vector. The second
  -- element of the inner Vector stores the position assigned to this edge. This position
  -- will be used to store the new point (correspondent to each edge) in the new list of
  -- points.
  xs <- VM.read vs a
  case U.find ((==b) . fst) xs of
   Just o -> return (snd o) -- edge already exist, return its correspondent position
   _      -> do
     -- get current free position
     !k <- readSTRef kref
     -- assigned edge to the current free position
     VM.write vs a (U.cons (b, k) xs)
     -- set counter to the next free position
     writeSTRef kref (k+1)
     return k
  where
    a = min i j
    b = max i j

scaleIsoSphere :: Double -> IsoSphere -> IsoSphere
scaleIsoSphere k iso = iso {vertices = U.map (k *&) (vertices iso)}

-- | Angular distance (Radians) between neighboring points at a given subdivision level.
angularDistance :: Int -> Double
angularDistance = (k /) . fromIntegral . ((2 :: Int)^) . abs
  where k = pi * (72 / 180)

genIsoSphereSO3Grid :: (Angle a)=> Symm -> a -> U.Vector Quaternion
genIsoSphereSO3Grid symm a = vs
  where
    (n, w) = getOptSubDivisionLevel a
    step   = fromAngle (w :: Rad)
    nstep  = floor (pi / step)
    iso    = isoSphere n
    minw   = getMinDistFZPlanes symm
    vs     = U.cons mempty (U.concatMap getLayer (U.enumFromStepN step step nstep))
    toQ t  = U.map (\v -> toQuaternion (mkAxisPair v (Rad t)))
    getLayer t
      | t <= minw = toQ t (vertices iso)
      | otherwise = U.filter (isInRodriFZ symm) (toQ t (vertices iso))

-- | Find the minimum subdivision level that provides the given angular step size.
getOptSubDivisionLevel :: (Angle a0, Angle a1)=> a0 -> (Int, a1)
getOptSubDivisionLevel a = go 0
  where
    w = fromAngle a
    go !n
      | o > w     = go (n + 1)
      | otherwise = (n, toAngle o)
      where o = angularDistance n

-- | Fast query to the nearest point. N-aray tree search with expected time complexity
-- O(l) where l is the subdivision level.
nearestPoint :: IsoSphere -> Vec3D -> (Int, Vec3D, Double)
nearestPoint IsoSphere{..} q = getClosest 0 fid0
  where
    lmax = V.length faces - 1
    fid0 = U.minIndex $ U.map getD (centers V.! 0)
    qn   = normalize q
    getD = acos . (qn &.)
    -- find the face position the closest child
    getNextChild !l !fid = let
      cs   = centers V.! (l+1)
      fo   = 4 * fid
      fids = U.enumFromStepN fo 1 4
      i    = U.minIndex $ U.map (getD . (cs U.!)) fids
      in fids U.! i
    -- recursively find the closest face
    getClosest !l !fid
      | l >= lmax = getClosestOnFace face
      | otherwise = getClosest (l+1) next
      where
        fs   = faces V.! l
        face = fs    U.! fid
        next = getNextChild l fid
    -- find the closest vertex on the closest face
    getClosestOnFace (ia, ib, ic) = let
      xs = map (\i -> let v = vertices U.! i in (i, v, getD v)) [ia, ib, ic]
      in L.minimumBy (compare `on` (\(_, _, d) -> d)) xs

-- | Render IsoSphere faces.
renderIsoSphereFaces :: IsoSphere -> [VTKAttrPointValue Vec3D] -> VTK Vec3D
renderIsoSphereFaces IsoSphere{..} attrs = mkUGVTK "IsoSphere" vertices fs attrs []
  where fs = if V.null faces then U.empty else V.last faces

-- | Render quaternion points.
renderQuaternions :: U.Vector Quaternion -> [VTKAttrPointValue Vec3D] -> VTK Vec3D
renderQuaternions qs attrs = mkUGVTK "IsoSpace" (U.map quaternionToVec3 qs) is attrs []
  where is = U.generate (U.length qs) id

quaternionToVec3 :: Quaternion -> Vec3D
quaternionToVec3 q = let
  (q0, qv) = splitQuaternion q
  omega = 2 * acos q0
  in omega *& normalize qv
