# function to lexicographically sort the points in the point-cloud to
# deal with permutation-invariance in the ensemble model which predicts the
# amount of generated power. The ensemble model has been trained on
# lexicographically sorted point-clouds, but it is not permutation-invariant
# with respect to its input. Thus, point-clouds must be soerted before calling
# the ensemble model

permutation_invariant_objective <- function(x) {
  X <- data.frame(x=matrix(x,ncol=2))
  ixs <- order(X$x.1,X$x.2)
  x <- as.vector(as.matrix(X[ixs,]))
  return( objective(x) )
}
