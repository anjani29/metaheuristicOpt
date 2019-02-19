

sphere <- function(X){
return(sum(X^2))
 }

Vmax <- 2
ci <- 1.5
cg <- 1.5
w <- 0.7
 numVar <- 5
rangeVar <- matrix(c(-10,10), nrow=2)
 
resultPSO <- PSO(sphere, optimType="MIN", numVar, numPopulation=20, 
                 maxIter=100, rangeVar, Vmax, ci, cg, w)

 optimum.value <- sphere(resultPSO)




PSO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, Vmax=2, ci=1.49445, cg=1.49445, w=0.729){
	# calculate the dimension of problem if not specified by user
	dimension <- ncol(rangeVar)

	# parsing rangeVar to lowerBound and upperBound
	lowerBound <- rangeVar[1,]
	upperBound <- rangeVar[2,]
	
	# if user define the same upper bound and lower bound for each dimension
	if(dimension==1){
		dimension <- numVar
	}

	## convert optimType to numerical form
	## 1 for minimization and -1 for maximization
	if(optimType == "MAX") optimType <- -1 else optimType <- 1

	# generate initial population of particle
	particles <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# calculate the initial local best
	# local best for this step is a global best
	Gbest <- calcBest(FUN, optimType, particles)
	Lbest <- particles

	# initial velocity of each particle
	velocity <- generateRandom(numPopulation, dimension, -Vmax, Vmax)

	# find the best particle position
	bestParticle <- engine.PSO(FUN, optimType, maxIter, lowerBound, upperBound, Vmax, ci, cg, w, Gbest, Lbest, particles, velocity)

	return(bestParticle)
}



engine.PSO <- function(FUN, optimType, maxIter, lowerBound, upperBound, Vmax, ci, cg, w, Gbest, Lbest, particles, velocity){
	FLbest <- calcFitness(FUN, optimType, Lbest)
	FGbest <- optimType*FUN(Gbest)
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
	for (t in 1:maxIter){
		for (i in 1:nrow(particles)){
			for (d in 1:ncol(particles)){
				# pick random rumber
				ri <- runif(1)
				rg <- runif(1)

				# update the particle velocity
				newV <- w * velocity[i,d] + ci*ri*(Lbest[i,d]-particles[i,d]) + cg*rg*(Gbest[d]-particles[i,d])
				
				# check range velocity
				if(newV < -Vmax) newV <- -Vmax
				if(newV > Vmax) newV <- Vmax
				velocity[i,d] <- newV

				newPos <- particles[i,d] + velocity[i,d]
				# check range search space
				if(length(lowerBound)==1){
					if(newPos < lowerBound) newPos <- lowerBound
					if(newPos > upperBound) newPos <- upperBound
				}else{
					if(newPos < lowerBound[d]) newPos <- lowerBound[d]
					if(newPos > upperBound[d]) newPos <- upperBound[d]
				}
				particles[i,d] <- newPos

				# check Local best and Global best
				F <- optimType*FUN(particles[i,])
				if(F < FLbest[i]){
					Lbest[i,] <- particles[i,]
					FLbest[i] <- F
					if(FLbest[i] < FGbest){
						Gbest <- Lbest[i,]
						FGbest <- FLbest[i]
					}
				}
			}
		}
		curve[t] <- FGbest
		setTxtProgressBar(progressbar, t)
	}
	close(progressbar)
	curve <- curve*optimType
	## plot(c(1:maxIter), curve, type="l", main="PSO", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  ## ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(Gbest)
}
