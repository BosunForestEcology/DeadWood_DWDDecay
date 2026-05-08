# Logistic survival function: probability of a piece remaining in its current DC at age t.
logisticPr <- function(t, a, b) 1 / (1 + exp(a + b * t))

# Vectorised: conditional probability of transitioning during [age_A, age_A + step]
# for each piece, given it has not yet transitioned by age_A.
#
# dc_vec          : integer vector of current DC (1-4) — selects logistic parameter row
# age_A_vec       : integer vector of ageSinceEntry at start of interval
# diameter_vec    : numeric vector of piece diameter_cm (must be >= 7.5 cm)
# logistic_params : data.table with columns (a0, a1, b0, b1), rows indexed by DC transition 1-4
# step            : interval length in years (default 5)
computeDWDTransProb <- function(dc_vec, age_A_vec, diameter_vec, logistic_params, step = 5L) {
  a   <- logistic_params$a0[dc_vec] + logistic_params$a1[dc_vec] * log(diameter_vec)
  b   <- logistic_params$b0[dc_vec] + logistic_params$b1[dc_vec] * log(diameter_vec)
  ab_A  <- a + b * age_A_vec
  ab_B  <- a + b * (age_A_vec + step)
  prc_A <- 1 - 1 / (1 + exp(ab_A))
  prc_B <- 1 - 1 / (1 + exp(ab_B))
  p     <- (prc_B - prc_A) / (1 - prc_A)
  p[prc_A >= 1 - 1e-10] <- 1.0
  pmax(0, pmin(1, p))
}
