#' Relevant Component Analysis
#'
#' @description
#' \code{rca} performs relevant component analysis (RCA) for the given data.
#' It takes a data set and a set of positive constraints as arguments
#' and returns a linear transformation of the data space into better
#' representation, alternatively, a Mahalanobis metric over the data space.
#'
#' The new representation is known to be optimal in an information
#' theoretic sense under a constraint of keeping equivalent data
#' points close to each other.
#'
#' @param x \code{n * d} matrix or data frame of original data.
#'
#' @param chunks a vector of size \code{N} describing the chunklets:
#' \code{-1} in the \code{i}-th place says that point \code{i} does not
#' belong to any chunklet; integer \code{j} in place \code{i} says
#' that point \code{i} belongs to chunklet \code{j};
#' The chunklets indexes should be \code{1:number-of-chunklets}.
#'
#' @param useD optional. When not given, RCA is done in the
#' original dimension and \code{B} is full rank. When \code{useD} is given,
#' RCA is preceded by constraints based LDA which reduces
#' the dimension to \code{useD}. \code{B} in this case is of rank \code{useD}.
#'
#' @return
#' A list of the RCA results:
#' \itemize{
#' \item \code{B}: The RCA suggested Mahalanobis matrix.
#' Distances between data points \code{x1}, \code{x2} should be
#' computed by \code{(x2 - x1)^T * B * (x2 - x1)}
#' \item \code{RCA}: The RCA suggested transformation of the data.
#' The data should be transformed by \code{RCA * data}
#' \item \code{newX}: The data after the RCA transformation.
#' \code{newX = data * RCA}
#' }
#'
#' @details
#' The three returned objects are just different forms of the same output.
#' If one is interested in a Mahalanobis metric over the original data space,
#' the first argument is all she/he needs. If a transformation into another
#' space (where one can use the Euclidean metric) is preferred, the second
#' returned argument is sufficient. Using \code{A} and \code{B} are equivalent
#' in the following sense:
#'
#' if \code{y1 = A * x1}, \code{y2 = A * y2}  then
#'
#' \code{(x2 - x1)^T * B * (x2 - x1) = (y2 - y1)^T * (y2 - y1)}
#'
#' @note
#' Note that any different sets of instances (chunklets),
#' e.g. \code{{1, 3, 7}} and \code{{4, 6}}, might belong to the
#' same class and might belong to different classes.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom stats cov
#'
#' @export rca
#'
#' @references
#' Aharon Bar-Hillel, Tomer Hertz, Noam Shental, and Daphna Weinshall (2003).
#' Learning Distance Functions using Equivalence Relations.
#' \emph{Proceedings of 20th International Conference on
#' Machine Learning (ICML2003)}
#'
#' @examples
#' library("MASS") # generate synthetic multivariate normal data
#' set.seed(42)
#' k <- 100L # sample size of each class
#' n <- 3L # specify how many classes
#' N <- k * n # total sample size
#' x1 <- mvrnorm(k, mu = c(-16, 8), matrix(c(15, 1, 2, 10), ncol = 2))
#' x2 <- mvrnorm(k, mu = c(0, 0), matrix(c(15, 1, 2, 10), ncol = 2))
#' x3 <- mvrnorm(k, mu = c(16, -8), matrix(c(15, 1, 2, 10), ncol = 2))
#' x <- as.data.frame(rbind(x1, x2, x3)) # predictors
#' y <- gl(n, k) # response
#'
#' # fully labeled data set with 3 classes
#' # need to use a line in 2D to classify
#' plot(x[, 1L], x[, 2L],
#'   bg = c("#E41A1C", "#377EB8", "#4DAF4A")[y],
#'   pch = rep(c(22, 21, 25), each = k)
#' )
#'
#' # generate synthetic chunklets
#' chunks <- vector("list", 300)
#' for (i in 1:100) chunks[[i]] <- sample(1L:100L, 10L)
#' for (i in 101:200) chunks[[i]] <- sample(101L:200L, 10L)
#' for (i in 201:300) chunks[[i]] <- sample(201L:300L, 10L)
#'
#' chks <- x[unlist(chunks), ]
#'
#' # make "chunklet" vector to feed the chunks argument
#' chunksvec <- rep(-1L, nrow(x))
#' for (i in 1L:length(chunks)) {
#'   for (j in 1L:length(chunks[[i]])) {
#'     chunksvec[chunks[[i]][j]] <- i
#'   }
#' }
#'
#' # relevant component analysis
#' rcs <- rca(x, chunksvec)
#'
#' # learned transformation of the data
#' rcs$RCA
#'
#' # learned Mahalanobis distance metric
#' rcs$B
#'
#' # whitening transformation applied to the chunklets
#' chkTransformed <- as.matrix(chks) %*% rcs$RCA
#'
#' # original data after applying RCA transformation
#' # easier to classify - using only horizontal lines
#' xnew <- rcs$newX
#' plot(xnew[, 1L], xnew[, 2L],
#'   bg = c("#E41A1C", "#377EB8", "#4DAF4A")[gl(n, k)],
#'   pch = c(rep(22, k), rep(21, k), rep(25, k))
#' )

rca = function(x, chunks, useD = NULL) {

  n = nrow(x)
  d = ncol(x)
  if(is.null(useD)) useD = d

  # subtract the mean
  TM = colMeans(x)
  x  = x - (matrix(1L, n, 1L) %*% TM)

  # compute chunklet means and center data
  S = max(chunks)

  Cdata = matrix()
  AllInds = vector()
  M = matrix(NA, nrow = S, ncol = d)
  tmp = vector('list', S)

  for (i in 1L:S) {
    inds = which(chunks == i)
    M[i, ] = colMeans(x[inds, ])
    tmp[[i]] = x[inds, ] - matrix(1L, length(inds), 1L) %*% M[i, ]
    Cdata = do.call(rbind, tmp)
    AllInds = c(AllInds, inds)
  }

  # Compute inner covariance matrix
  InnerCov = cov(Cdata) * ((nrow(Cdata) - 1L) / nrow(Cdata))

  # Optional cFLD: find optimal projection: min | A S_w A^t | / | A S_t A^t |

  if (useD < d) {
    # Compute total covariance using only chunkleted points
    TotalCov = cov(x[AllInds, ])
    # TotalCov = cov(x)  # Compute total covariance using all the data.
    # More accurate but may lead to PSD problems
    tmp = eigen(solve(TotalCov) %*% InnerCov)
    D = tmp$values
    V = tmp$vectors
    # reorder the vectors in descending order
    V = V[ , ncol(V):1L]
    # A is the cFLD transformation
    # Acts on row data by multiplication from the right
    A = V[ , 1L:useD]
    InnerCov = t(A) %*% InnerCov %*% A
  } else {
    A = diag(d)
  }

  # RCA: whiten the data w.r.t the inner covariance matrix
  tmp = svd(InnerCov)
  U1  = tmp$u
  S1  = diag(tmp$d)
  V1  = tmp$v

  # Operate from the right on row vectors
  RCA  = as.matrix(A %*% (U1 %*% S1^(0.5)))

  # The total mean subtracted is re-added before transformation.
  # This operation is not required for distance computations,
  # as it only adds a constant to all points
  newX = as.matrix(x + matrix(1L, n, 1L) %*% TM) %*% RCA

  B = RCA %*% t(RCA)

  res = list('B' = B, 'RCA' = RCA, 'newX' = newX)
  res

}
