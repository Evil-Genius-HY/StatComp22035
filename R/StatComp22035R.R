#' @title Example data of an adjacent matrix 
#' @name MA
#' @description Example data used to illustrate the performance of \code{PageRank}.
#' @examples
#' \dontrun{
#' data(MA)
#' attach(MA)
#' NumNode <- 50
#' nrow <- ncol <- NumNode
#' IniVec <- runif(NumNode, 0, 1)
#' Score <- PageRank(MA, 0.85, 100000, IniVec)
#' }
NULL

#' @title Hybrid Particle Swarm Optimization Algorithm
#' @description Improved PSO Algorithm to handle complex optimization problem
#' @param solveFun Target function
#' @param limitX The limit of the independent variable
#' @param particleNum Number of particles
#' @param w Weight of inertia
#' @param c1 Constant of acceleration
#' @param c2 Constant of acceleration
#' @param iters Maximum number of iterations
#' @param alpha Value increment threshold of the optimal fitness value
#' @return The optimal (maximum) solution to the target function
#' @examples
#' \dontrun{
#' solveFun <- function(x){
#'   log(sin(pi*x)+2)
#' }
#' limitX <- c(-1, 2) 
#' particleNum <- 20 
#' w <- 1          
#' c1 <- c2 <- 2   
#' iters <- 10000  
#' alpha <- 0.0002
#' solution <- Hybrid.PSO(solveFun,particleNum,limitX,w,c1,c2,iters,alpha)
#' }
#' @importFrom stats runif
#' @export
Hybrid.PSO <- function(solveFun,particleNum,limitX,w,c1,c2,iters,alpha){
  pbest <- NULL
  gbest <- NULL
  gbestAdd <- NULL
  vmax <- 0.15 * (limitX[2] - limitX[1])
  xMat <- matrix(c(x = runif(particleNum, limitX[1], limitX[2])), dimnames = list(NULL, c("x")))
  vMat <- matrix(c(x = runif(particleNum, -vmax, vmax)),  dimnames = list(NULL, c("x")))
  adjusts <- apply(xMat, 1, solveFun)
  pbest <- cbind(xMat, adjusts)
  idxAdjusts <- ncol(pbest)
  gbest <- pbest[which.max(pbest[, idxAdjusts]),]
  for (i in 1:iters){
    mapply(function(no, adj){
      if(adj > pbest[no, idxAdjusts]){
        pbest[no, ] <<- c(xMat[no, ], adj)
      }
    }, 1:length(adjusts), adjusts)
    
    if (max(pbest[, idxAdjusts]) > gbest[idxAdjusts]) {
      gbestAdd <- max(pbest[, idxAdjusts]) - gbest[idxAdjusts]
      gbest <- pbest[which.max(pbest[, idxAdjusts]), ]
    }
    xMatOld <- xMat
    xMat <- xMat + vMat
    vMat <- w*vMat +  
      c1 * runif(1, 0, 1) * (pbest[, 1:(idxAdjusts - 1), drop=F] - xMatOld) + 
      c2 * runif(1, 0, 1) * (matrix(rep(gbest[1:(idxAdjusts - 1)], particleNum), ncol = idxAdjusts - 1 , byrow = T)-xMatOld) 
    vMat[which(vMat < (-vmax))] <- (-vmax)
    vMat[which(vMat > (vmax))]  <- (vmax)
    xMat[which(xMat < (limitX[1]))] <- (limitX[1])
    xMat[which(xMat > (limitX[2]))] <- (limitX[2])
    adjusts <- apply(xMat, 1, solveFun)
    if (!is.null(gbestAdd) && gbestAdd < alpha) {
      return(gbest)
    }
  }
  cat("The iteration limit was reached and the result could not be solved")
}

#' @title Import packages for homework review file.
#' @name HomeworkPackage
#' @description Import packages for homework review file
#' @examples
#' \dontrun{
#' library(Rcpp)
#' library(microbenchmark)
#' library(knitr)
#' library(rmarkdown)
#' library(boot)
#' library(bootstrap)
#' library(DAAG) 
#' library(MASS)
#' library(ggplot2)
#' }
#' @import microbenchmark
#' @import Rcpp
#' @import knitr
#' @import rmarkdown
#' @import boot
#' @import bootstrap
#' @import DAAG
#' @import MASS
#' @import ggplot2
NULL
