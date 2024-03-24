#' Fit Non-linear Partial Least Squares Regression (NPLSR) Model for Binary
#' Response
#'
#' This function fits a Non-linear Partial Least Squares Regression
#' (NPLSR) model for binary response variables.
#'
#' @param Y Response matrix (binary) of dimension I x L.
#' @param X Matrix of predictors (explanatory variables) of dimension
#' I x J*K.
#' @param I Number of observations.
#' @param J Number of columns in the first dimension of the original
#' matrix of predictors.
#' @param K Number of columns in the second dimension of the original
#' matrix of predictors.
#' @param L Number of response variables.
#' @param S Number of latent variables (default is 2).
#' @param tolerance Tolerance value for convergence (default is 0.000005).
#' @param maxiter Maximum number of iterations (default is 100).
#' @param penalization Penalization parameter (default is 0.1).
#' @param OptimMethod Optimization method (default is "CG").
#'
#' @return A list containing various components of the NPLSR model including
#' latent variables, loadings, predictions, and percent correct classifications.
#'
#' @export
#'
#' @examples
#' pls <-NPLSRBinFit(Y,X,8,7,6,9)
#' @import grLogBiplotRegBRec
#' @import grLogBiplotRegARec
#' @import JLogBiplotRegBRec
#' @import JLogBiplotRegARec
#' @import RidgeBinaryLogistic
#'
NPLSRBinFit <- function(Y, X, I, J, K, L, S=2, tolerance=0.000005, maxiter=100,
                       penalization=0.1, OptimMethod="CG")
{

TT = matrix(0, I, S)
P = matrix(0,J*K,S)
WJ = matrix(0,J,S)
WK = matrix(0,K,S)
Q = matrix(0, L, S)
q0 = matrix(0, nrow = L, ncol = 1)
for (l in 1:L) q0[l] = RidgeBinaryLogistic(y = Y[, l], matrix(1,I, 1), penalization = 0)$beta
Q = q0
U = matrix(1, I, 1)
# We suppose thet the X variables are at least centered

for (s in 1:S){
  us=matrix(X[,1],I,1)
  U=cbind(U,us)
  parQ=rep(1,L)/sqrt(L)
  Q=cbind(Q,parQ)
  error=1
  iter=0
  ts = us
  while ((error>tolerance) & (iter<maxiter)){
    iter=iter+1
    told=ts
    Zs <- crossprod(us, X)
    Zs <- matrix(Zs, nrow = J, ncol = K)
    svd.z <- svd(Zs)
    wjs <- svd.z$u[, 1]
    wks <- svd.z$v[, 1]
    ws = kronecker(wks, wjs)
    # Actualiza ts
    ts <- X %*% ws
    us=ts
    U[,s+1]=us
    # Update Q
    resbipQ <- optim(parQ, fn=JLogBiplotRegBRec, gr=grLogBiplotRegBRec, method=OptimMethod, X=Y, A=U, B=Q,lambda=penalization)
    parQ=resbipQ$par
    Q[,s+1]=parQ
    # Update U
    resbipU <- optim(us, fn=JLogBiplotRegARec, gr=grLogBiplotRegARec, method=OptimMethod, X=Y, A=U, B=Q, lambda=penalization)
    us=resbipU$par
    U[,s+1]=us
    #ts=us
    error=sum((told-ts)^2)
  }
  TT[,s]=ts
  WJ[,s]=wjs
  WK[,s]=wks
  P[,s]=ws
  X = X - ts %*% t(ws)
}
U=U[,-1]
Q=Q[,-1]

Lin= cbind(1, TT) %*% t(cbind(q0,Q))
Expected=exp(Lin)/(1+exp(Lin))
C=solve(t(TT)%*%TT)%*%t(TT)%*%U
B= P%*%C%*%t(Q)
Pred=matrix(as.numeric(Expected>0.5), nrow=I)
Right=(Y==Pred)
PercentCorrect=sum(Right)/(I*L)
PercentCorrectCols=apply(Right, 2, sum)/I
result=list(TT=TT, P=P, U=U, Q=Q, q0=q0, B=B, WJ= WJ, WK=WK, Linterm=Lin, Expected=Expected, Predictions=Pred, PercentCorrect=PercentCorrect,
            PercentCorrectCols=PercentCorrectCols)
return(result)
}

