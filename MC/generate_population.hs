import Control.Monad.Random (evalRandIO)
import Control.Monad.Random.Class (MonadRandom, getRandomR)
import Data.List (sort, group, intersperse, foldl')
import Control.Monad (liftM)
import System.Environment (getArgs)

sampleCDF invcdf = do
	p <- getRandomR (0.0, 1.0)
	return $ invcdf p

-- exponential:
-- p(l, x) = l exp(-l x)
-- c(l, x) = 1 - exp(-l x)
-- c^-1(l, p) = ln(1-p) / (-l)
expInvCDF lambda p = (log (1 - p)) / (-lambda)

data CellType = A | B | C deriving Show
data (Real a) => Cell a = Cell {cellType:: ! CellType, cellLifetime:: ! a}
	deriving Show

data Clone a = Clone ! (Cell a) [Clone a]
	deriving Show

data Parameters = Parameters {lambda:: ! Double, r:: ! Double, rho:: ! Double, time::Double}

generateClone _  C = return $ Clone (Cell C (1.0/0.0)) []
generateClone ps B = do
	l <- sampleCDF $ expInvCDF gamma
	c <- generateClone ps C
	return $ Clone (Cell B l) [c]
	where
		gamma = l * rh / (1-rh)
		l = lambda ps
		rh = rho ps
generateClone ps A = do
	l <- sampleCDF $ expInvCDF (lambda ps)
	r' <- getRandomR (0.0, 1.0)
	cs <- if r' < (r ps) then
		do c1 <- generateClone ps A
		   c2 <- generateClone ps A
		   return $ [c1, c2]
	     else if r' > 1-(r ps) then
		do c1 <- generateClone ps B
		   c2 <- generateClone ps B
		   return $ [c1, c2]
	     else
		do c1 <- generateClone ps A
		   c2 <- generateClone ps B
		   return $ [c1, c2]
	return $ Clone (Cell A l) cs

countBS t (Clone c cs) =
	if l > t then
		case cellType c of
			A -> (1,0)
			B -> (1,0)
			C -> (0,1)
	else
		foldl' sumPair (0,0) $ map (countBS (t - l)) cs
	where
		l = cellLifetime c
		sumPair (a,b) (a',b') = (a+a',b+b')

generateBS ps = do
	c <- generateClone ps A
	return $ countBS (time ps) c

aggregate :: (Num t, Ord t) => [(t,t)] -> [(t,t,Int)]
aggregate = map (\rs@((b,s):_) -> (b,s,length rs)) . group . sort . filter (\(b,s) -> (b+s)>1)

generateData ps n = liftM aggregate $ sequence $ replicate n $ generateBS ps

showTriple (b,s,n) = (show b) ++ " " ++ (show s) ++ " " ++ (show n)

main = do
	args <- getArgs
	let n = read (args !! 0)
	let ps = Parameters 0.25 0.25 0.6 6
	agg <- evalRandIO $ generateData ps n
	--print $! agg
	putStrLn $ "[" ++ concat (intersperse "; " $ map showTriple $! agg) ++ "]"
