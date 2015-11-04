nm <- function(x) {
   # Make sure input is an integer matrix or vector.
   if (is.data.frame(x)) x <- as.matrix(x)
   stopifnot(is.numeric(x))
   if (any(x < 0)) stop("x contains negative values")

   if (is.vector(x)) {
      ix <- as.integer(x)
      d <- 1L
      n <- length(ix)
   }
   else if (is.matrix(x)) {
      # The data must be transposed.
      ix <- apply(t(x), c(1,2), as.integer)
      d <- nrow(ix)
      n <- ncol(ix)
   }
   else {
      stop("x must be a matrix or a vector")
   }
   if (any((t(x)-ix) != 0 )) {
      warning("non-integer x")
   }
   
   # Allocate return value and call C function.
   return(.Call("R_call_mle_nm", ix, d, n))

}
