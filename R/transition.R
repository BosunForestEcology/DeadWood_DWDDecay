# Draw next-year DC for each element of DC_vec using a Markov transition matrix.
# transMatrix rows are current DC, columns are next-year DC.
# Row probabilities need not sum to 1; rmultinom normalises.
applyTransition <- function(DC_vec, transMatrix) {
  if (length(DC_vec) == 0L) return(integer(0L))
  stopifnot(
    all(!is.na(DC_vec)),
    all(DC_vec >= 1L & DC_vec <= nrow(transMatrix))
  )
  vapply(DC_vec, function(dc) {
    probs <- transMatrix[dc, ]
    if (all(probs == 0)) stop(sprintf("applyTransition: transition row %d is all zeros", dc))
    which(stats::rmultinom(1L, 1L, probs) == 1L)
  }, integer(1L))
}
