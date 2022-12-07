## ----eval=FALSE---------------------------------------------------------------
#  Hybrid.PSO <- function(solveFun,particleNum,limitX,w,c1,c2,iters,alpha){
#    pbest <- NULL
#    gbest <- NULL
#    gbestAdd <- NULL
#    vmax <- 0.15 * (limitX[2] - limitX[1])
#    xMat <- matrix(c(x = runif(particleNum, limitX[1], limitX[2])), dimnames = list(NULL, c("x")))
#    vMat <- matrix(c(x = runif(particleNum, -vmax, vmax)),  dimnames = list(NULL, c("x")))
#    adjusts <- apply(xMat, 1, solveFun)
#    pbest <- cbind(xMat, adjusts)
#    idxAdjusts <- ncol(pbest)
#    gbest <- pbest[which.max(pbest[, idxAdjusts]),]
#    for (i in 1:iters){
#      mapply(function(no, adj){
#        if(adj > pbest[no, idxAdjusts]){
#          pbest[no, ] <<- c(xMat[no, ], adj)
#        }
#      }, 1:length(adjusts), adjusts)
#  
#      if (max(pbest[, idxAdjusts]) > gbest[idxAdjusts]) {
#        gbestAdd <- max(pbest[, idxAdjusts]) - gbest[idxAdjusts]
#        gbest <- pbest[which.max(pbest[, idxAdjusts]), ]
#      }
#      xMatOld <- xMat
#      xMat <- xMat + vMat
#      vMat <- w*vMat +
#        c1 * runif(1, 0, 1) * (pbest[, 1:(idxAdjusts - 1), drop=F] - xMatOld) +
#        c2 * runif(1, 0, 1) * (matrix(rep(gbest[1:(idxAdjusts - 1)], particleNum), ncol = idxAdjusts - 1 , byrow = T)-xMatOld)
#      vMat[which(vMat < (-vmax))] <- (-vmax)
#      vMat[which(vMat > (vmax))]  <- (vmax)
#      xMat[which(xMat < (limitX[1]))] <- (limitX[1])
#      xMat[which(xMat > (limitX[2]))] <- (limitX[2])
#      adjusts <- apply(xMat, 1, solveFun)
#      if (!is.null(gbestAdd) && gbestAdd < alpha) {
#        return(gbest)
#      }
#    }
#    cat("The iteration limit was reached and the result could not be solved")
#  }
#  

## ---- fig.align='center'------------------------------------------------------
solveFun <- function(x){
  x*(sin(10*pi*x)+2)
}
limitX <- c(-1, 2) 
particleNum <- 20 
w <- 1          
c1 <- c2 <- 2   
iters <- 10000  
alpha <- 0.0002
solution <- Hybrid.PSO(solveFun,particleNum,limitX,w,c1,c2,iters,alpha)
cat("The optimal solution of the target function is x=", solution[1], ", f(x)=",solution[2], ".")
x <- seq(-1, 2, 0.01)
plot(x, solveFun(x), type = "l")
grid(col='grey60')
points(solution[1],solution[2],col="red")

## ----eval=FALSE---------------------------------------------------------------
#  NumericVector PageRank(NumericMatrix MA, double d, int MaxIter, NumericVector IniVec){
#  
#    int N = MA.nrow();
#    NumericMatrix MT(N, N);
#    NumericVector UpdVec(N);
#  
#    for(int i = 0; i < N; i++){
#      double sum = 0;
#      for(int j = 0; j < N; j++){
#        sum = sum + MA(j,i);
#      }
#      if(sum>0){
#        for(int k = 0; k < N; k++){
#          MT(k,i) = MA(k,i)/sum;
#        }
#      }
#    }
#  
#    for(int i = 0; i < MaxIter; i++){
#      for(int j = 0; j < N; j++){
#        UpdVec[j] = 0;
#        for(int k = 0; k < N; k++){
#          UpdVec[j] = UpdVec[j] + d*MT(j,k)*IniVec[k];
#        }
#        UpdVec[j] = UpdVec[j] + (1-d)/N;
#      }
#      double Diff = 0;
#      for(int l = 0; l < N; l++){
#        Diff = Diff + abs(UpdVec[l]-IniVec[l]);
#      }
#      if(Diff < 1e-10){
#        break;
#      }
#      for(int m = 0; m < N; m++){
#        IniVec[m] = UpdVec[m];
#      }
#    }
#  
#    return(UpdVec);
#  }

## ----fig.align='center'-------------------------------------------------------
library(Rcpp)
library(ggplot2)
library(StatComp22035)

data(MA)
NumNode <- 50
nrow <- ncol <- NumNode
IniVec <- runif(NumNode, 0, 1)
Score <- PageRank(MA, 0.85, 100000, IniVec)
data <- data.frame(
  index = 1:NumNode,  
  value = Score
)
ggplot(data, aes(x=index, y=value)) + 
  geom_bar(stat = "identity", width=0.5)


