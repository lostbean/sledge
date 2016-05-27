{-# LANGUAGE
    FlexibleContexts
  , FlexibleInstances
  , RecordWildCards
  #-}
module Texture.Sampler
  ( hitAndRunSlice
  , HRSCfg (..)
  , defaultCfg
  ) where

import Data.IORef
import Linear.Vect
import System.Random

import Texture.Orientation

-- | Function to obtain a linear shoter.
class HasRandomDir a where
  linearShoter :: a -> IO (Double -> a)

data HRSCfg
  = HRSCfg
  { initShotDist     :: Double  -- ^ Initial absolute distance
  , mixingFraction   :: Double  -- ^ Mixing factor. 0 means always use 'initShotDist' and
                                -- 1 means always use last range
  , longShotFraction :: Double  -- ^ Fraction of frequency of long shot ('maxShotDist') tries.
  , maxShotDist      :: Double  -- ^ Maximum absolute distance
  } deriving (Show)

-- | Default configuration for 'hitAndRunSlice' with 'initShotDist' = pi/18,
-- 'mixingFraction' = 0.01, 'longShotFraction' = 0.05 and 'maxShotDist' = pi.
defaultCfg :: HRSCfg
defaultCfg = HRSCfg
  { initShotDist     = pi/18
  , mixingFraction   = 0.05
  , longShotFraction = 0.05
  , maxShotDist      = pi
  }

hitAndRunSlice :: (HasRandomDir a)=> HRSCfg -> (a -> Double) -> a -> Int -> IO [a]
hitAndRunSlice cfg@HRSCfg{..} p x0 n = do
  -- initialize counter and dynamic range
  count <- newIORef (0 :: Int)
  kdyn  <- newIORef initShotDist
  -- run sampler
  xs <- sampler x0 0 (count, kdyn)
  -- print info
  readIORef count >>= print
  readIORef kdyn  >>= print
  return xs
  where
    -- shrink linear range with k
    shrinkrange (kmin, kmax) k
      | k > 0     = (kmin, k)
      | otherwise = (k, kmax)

    sampler xi i (count, kdyn)
      | i >= n    = return []
      | otherwise = do
        -- get level of slice
        u      <- randomRIO (0, p xi)
        -- get a linear shoter (a infinite line on the sliced plane)
        shoter <- linearShoter xi

        -- get a successful shot
        let
          scanline range = do
            -- random shot in the range
            k <- randomRIO range
            let xj = shoter k
            let fj = p xj
            modifyIORef count (+1)
            if fj >= u

              then do
              -- Good shot, then update dynamic range
              rangeMixer cfg kdyn range
              -- get rest and add value to result list
              xs <- sampler xj (i+1) (count, kdyn)
              return (xj : xs)

              -- Bad shot, shrink range and try again
              else scanline (shrinkrange range k)

        -- get initial range with both extremes lying outside and using dynamic range
        r <- getInitRange cfg count p shoter u =<< readIORef kdyn
        --readIORef kdyn >>= putStrLn . ("Kdyn: " ++) . show
        -- get a successful shot
        scanline r

-- | Mix successful (accepted) range and average distance value.
rangeMixer :: HRSCfg -> IORef Double -> (Double, Double) -> IO ()
rangeMixer HRSCfg{..} km (rl, ru) = let
  m    = min 1 (abs mixingFraction)
  avgk = 0.5 * (abs rl + abs ru)
  in modifyIORef km $ \k -> (1 - m) * k + m * avgk

-- | Grow range progressively until both extremes are out (p(x) > level). Frequently tries
-- to explore the maximum range.
getInitRange :: HRSCfg -> IORef Int -> (a -> Double) -> (Double -> a) ->
                Double -> Double -> IO (Double, Double)
getInitRange HRSCfg{..} count p shoter y k0 = do
  x <- randomRIO (0, 1)
  if x < longShotFraction
    then maxShot
    else do
    kl <- go (-k0) (-k0)
    ku <- go k0 k0
    return (kl, ku)
  where
    maxShot = return (-maxShotDist, maxShotDist)
    go step k
      | abs k > maxShotDist = return (maxShotDist * signum k)
      | p (shoter k) < y    = modifyIORef count (+1) >> return k
      | otherwise           = modifyIORef count (+1) >> go step (k + step)

instance HasRandomDir Double where
  linearShoter x0 = do
    t <- func <$> randomRIO (0,1)
    return (\k -> x0 + k * t)
    where
      func :: Double -> Double
      func x = if x >= 0.5 then 1 else (-1)

instance HasRandomDir (Vec2 Double) where
  linearShoter x0 = do
    t <- sampleOne >>= func
    return (\k -> x0 &+ k *& t)
    where
      func v
        | l > 1     = sampleOne >>= func
        | l == 0    = return v
        | otherwise = return (v &* (1/l))
        where l = norm v
      sampleOne = do
        x <- randomRIO (-1,1)
        y <- randomRIO (-1,1)
        return (Vec2 x y)

instance HasRandomDir Quaternion where
  linearShoter x0 = do
    t <- sampleOne >>= func
    return (\k -> x0 #<= (toQuaternion . mkAxisPair t $ Rad k))
    where
      func v
        | l > 1     = sampleOne >>= func
        | l == 0    = return v
        | otherwise = return (v &* (1/l))
        where l = norm v
      sampleOne = do
        x <- randomRIO (-1,1)
        y <- randomRIO (-1,1)
        z <- randomRIO (-1,1)
        return (Vec3 x y z)
