#' Perform Non-linear Partial Least Squares Regression (NPLSR) for
#' Binary Responses
#'
#' This function conducts Non-linear Partial Least Squares Regression
#' (NPLSR) for binary response variables.
#'
#' @param Y Response matrix (binary) of dimension I x L.
#' @param X Matrix of predictors (explanatory variables) of dimension
#'  I x J*K.
#' @param I Number of observations.
#' @param J Number of columns in the first dimension of the original
#' matrix of predictors.
#' @param K Number of columns in the second dimension of the original
#' matrix of predictors.
#' @param L Number of response variables.
#' @param inames Names for the observations (optional, default is NULL).
#' @param xnames Names for the predictors (optional, default is NULL).
#' @param ynames Names for the response variables (optional, default is NULL).
#' @param S Number of latent variables (default is 2).
#' @param tolerance Tolerance value for convergence (default is 0.00005).
#' @param maxiter Maximum number of iterations (default is 100).
#' @param penalization Penalization parameter (default is 0.1).
#' @param OptimMethod Optimization method (default is "CG").
#'
#' @return An object of class 'PLSRBin' containing various components of the
#' NPLSR model including scores, loadings, coefficients, predictions, and percent correct classifications.
#'
#' @export
#'
#' @examples
#' pls <-NPLSRBin(Y,X,8,7,6,9)
#' plsbip=Biplot.PLSRBIN(pls)
#' plot(plsbip)
#'
#' @import NPLSRBinFit

NPLSRBin <- function(Y, X, I, J, K ,L, inames=NULL, xnames=NULL, ynames=NULL,
                    S=2, tolerance=0.00005,maxiter=100, penalization=0.1,
                    OptimMethod="CG")
{


  if (!CheckBinaryVector(Y)) stop("The response must be binary (0 or 1)")

  dimnames=paste("Comp.", 1:S)
  if (is.null(inames)){
    inames = paste("i",1:I,sep="")
  }
  if (is.null(xnames)){
    for (k in 1:K)
      for (j in 1:J){
        xnames=c(xnames, paste(k, j,sep="."))
      }
  }
  if (is.null(ynames)){
    ynames = paste("Y",1:L,sep="")
  }

  result=list()
  result$Method="PLSR for binary responses"
  result$X=X
  result$Y=Y
  result$tolerance=tolerance
  result$maxiter=maxiter
  result$penalization=penalization
  myfit=NPLSRBinFit(Y=Y, X=X, I, J, K, L, S=S, tolerance=tolerance, maxiter=maxiter, penalization=penalization)

  rownames(myfit$TT)=inames
  colnames(myfit$TT)=dimnames
  rownames(myfit$B)=xnames
  colnames(myfit$B)=ynames
  rownames(myfit$P)=xnames
  colnames(myfit$P)=dimnames

  result$XScores=myfit$TT
  result$XLoadings=myfit$P
  result$YScores=myfit$U
  result$YLoadings=myfit$Q
  rownames(result$YLoadings)=ynames
  colnames(result$YLoadings)=paste("Dim", 1:S)
  result$Coefficients=myfit$B
  result$XStructure=cor(result$X,myfit$TT)
  result$BinaryFits=myfit$fit
  result$Intercepts=myfit$q0
  result$LinTerm=myfit$Linterm
  result$Expected=myfit$Expected
  result$Predictions=myfit$Predictions
  rownames(result$Predictions)=inames
  colnames(result$Predictions)=ynames
  result$PercentCorrect=myfit$PercentCorrect
  result$PercentCorrectCols=myfit$PercentCorrectCols
  result$Initial_Transformation = 4
  result$ScaledX = X
  class(result)="PLSRBin"
  return(result)
}


