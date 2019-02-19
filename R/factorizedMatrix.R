methods::setClass("FactorizedMatrix", slots = list(A = "matrix", B = "matrix", scale = "numeric", offset = "matrix"))

FactorizedMatrix <- function(A, B) {
  if (!is.matrix(A) || !is.matrix(B)) stop("A and B must be matrices.")
  if (ncol(A) != nrow(B)) stop("Non-conformable arguments.")
  new("FactorizedMatrix", A=A, B=B, scale=1)
}

methods::setMethod(
  '[',
  c('FactorizedMatrix', 'numeric', 'missing'),
  function(x, i) {
    result = x@scale * as.vector(x@A[i,] %*% x@B)
    if (length(x@offset)>0) result = result + x@offset[i,]
    result
  }
)

methods::setMethod(
  '[',
  c('FactorizedMatrix', 'missing', 'numeric'),
  function(x, i, j) {
    result = x@scale * as.vector(x@A %*% x@B[,j])
    if (length(x@offset)>0) result = result + x@offset[,j]
    result
  }
)

methods::setMethod(
  '[',
  c('FactorizedMatrix', 'numeric', 'numeric'),
  function(x, i, j) {
    result = x@scale * as.numeric(x@A[i,] %*% x@B[,j])
    if (length(x@offset)>0) result = result + x@offset[i, j]
    result
  }
)

methods::setMethod(
  't',
  c('FactorizedMatrix'),
  function(x) {
    new("FactorizedMatrix", A=t(x@B), B=t(x@A), scale=x@scale, offset=t(x@offset))
  }
)

methods::setMethod(
  '-',
  c('matrix', 'FactorizedMatrix'),
  function(e1, e2) {
    if (!all(dim(e1)==dim(e2))) {
      stop("Error: non-conformable arrays")
    }
    if (length(e2@offset)>0) {
      offset = e2@offset + e1
    } else{
      offset = e1
    }
    new("FactorizedMatrix", A=e2@A, B=e2@B, scale=-e2@scale, offset=offset)
  }
)

methods::setMethod(
  'nrow',
  c('FactorizedMatrix'),
  function(x) {
    nrow(x@A)
  }
)

methods::setMethod(
  'ncol',
  c('FactorizedMatrix'),
  function(x) {
    ncol(x@B)
  }
)

methods::setMethod(
  'dim',
  c('FactorizedMatrix'),
  function(x) {
    c(nrow(x@A), ncol(x@B))
  }
)

methods::setMethod(
  'rowMeans',
  c('FactorizedMatrix'),
  function(x) {
    apply(apply(x@A, 1, function(row) row %*% x@B), 2, function(col) mean(col))
  }
)

#' @exportMethod as.matrix
methods::setMethod(
  'as.matrix',
  c('FactorizedMatrix'),
  function(x, typeof, ...) {
    if (missing(typeof)) typeof <- "numeric"
    result = x@scale * (x@A %*% x@B)
    if (length(x@offset)>0) result = result + x@offset
    result
  }
)

methods::setMethod(
  'crossprod',
  c('FactorizedMatrix'),
  function(x) {
    # TODO: Do this iteratively?
    crossprod(as.matrix(x))
  }
)
