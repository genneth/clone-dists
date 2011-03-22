{-# LANGUAGE TypeFamilies, EmptyDataDecls, ViewPatterns, FlexibleContexts #-}

import MC0DFramework
import Control.Monad
import Control.Monad.Random
import Data.List (intersperse)
import Control.Monad.State.Strict
import Control.Monad.Writer

data BrdU = BrdU -- unit division and stratification rate

instance DivisionProcess BrdU where
  data CellType BrdU = I Double | A0 Double | B0 Double | A1 Double | B1 Double
    deriving Show
  data MeanFields BrdU = BrdUMF
  timeToLive (I ttl) = ttl
  timeToLive (A0 ttl) = ttl
  timeToLive (B0 ttl) = ttl
  timeToLive (A1 ttl) = ttl
  timeToLive (B1 ttl) = ttl
  progeny _ _ (I _) = do
    a <- liftM A0 (exponentialVariable 1.0)
    b <- liftM B0 (exponentialVariable 1.0)
    return [a,b]
  progeny _ _ (A0 _) = do
      fate <- getRandomR (0.0, 1.0)
      if fate < r
        then pair cA cA
        else if fate < (1.0 - r)
          then pair cA cB
          else pair cB cB
    where
      r = 0.25 :: Double
      cA = liftM A1 (exponentialVariable 1.0)
      cB = liftM B1 (exponentialVariable 1.0)
      pair c1 c2 = do
        a <- c1
        b <- c2
        return [a, b]
  progeny a b (A1 t) = progeny a b (A0 t)
  progeny _ _ (B0 _) = return []
  progeny _ _ (B1 _) = return []
  cellType (I _) = 4
  cellType (A0 _) = 0
  cellType (B0 _) = 1
  cellType (A1 _) = 2
  cellType (B1 _) = 3
  evalMeanFields _ = BrdUMF

main = do
    (_, _, log) <- evalRandStateIO 
      (runPopulationWithRecording 5.0 stats)
      (0.0, BrdU, addCells emptyCulture 0.0 $ replicate 100 (I 0))
    forM log $ \(t,h) -> do
      putStr $ (show t) ++ " "
      sequence_ $ intersperse (putStr " ") (map (putStr . show) $ h)
      putStrLn ""
  where
    stats = do
      (t,_,cells) <- get
      tell $ [(t, [count 0 cells, count 1 cells, count 2 cells, count 3 cells])]
      return ()
