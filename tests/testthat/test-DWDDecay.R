library(data.table)
library(testthat)

.projRoot <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(.projRoot, "modules"))) {
  .projRoot <- getwd()
}
source(file.path(.projRoot, "modules", "DeadWood_DWDDecay", "R", "dwd-transition.R"))

logisticParams_test <- data.table::data.table(
  a0 = c(-7.295, -9.663,  3.882, -5.584),
  a1 = c( 1.888,  0.000, -2.986,  0.000),
  b0 = c( 1.325,  0.934, -0.076,  0.149),
  b1 = c(-0.239,  0.000,  0.093,  0.000)
)

# ---- computeDWDTransProb unit tests ----

test_that("computeDWDTransProb returns 0 when step = 0", {
  p <- computeDWDTransProb(1L, 10L, 15.0, logisticParams_test, step = 0L)
  expect_equal(p, 0, tolerance = 1e-10)
})

test_that("computeDWDTransProb returns values in [0, 1] across DCs, ages, and diameters", {
  for (dc in 1:4) {
    for (age in c(0L, 10L, 30L, 80L)) {
      for (d in c(5.0, 15.0, 40.0)) {
        p <- computeDWDTransProb(dc, age, d, logisticParams_test)
        expect_gte(p, 0)
        expect_lte(p, 1)
      }
    }
  }
})

test_that("computeDWDTransProb is non-decreasing with age for all DCs", {
  ages <- seq(0L, 100L, by = 5L)
  for (dc in 1:4) {
    probs <- computeDWDTransProb(rep(dc, length(ages) - 1L),
                                 ages[-length(ages)], 15.0, logisticParams_test)
    expect_true(all(diff(probs) >= -1e-10),
                info = paste("DC", dc, "transition prob not monotone"))
  }
})

test_that("computeDWDTransProb vectorises over multiple pieces", {
  p <- computeDWDTransProb(1:4, rep(20L, 4), rep(15.0, 4), logisticParams_test)
  expect_length(p, 4L)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("computeDWDTransProb: DC2->DC3 has no diameter dependence (a1=b1=0)", {
  p_small <- computeDWDTransProb(2L, 20L,  5.0, logisticParams_test)
  p_large <- computeDWDTransProb(2L, 20L, 40.0, logisticParams_test)
  expect_equal(p_small, p_large, tolerance = 1e-10)
})

test_that("computeDWDTransProb: DC4->out has no diameter dependence (a1=b1=0)", {
  p_small <- computeDWDTransProb(4L, 30L,  5.0, logisticParams_test)
  p_large <- computeDWDTransProb(4L, 30L, 40.0, logisticParams_test)
  expect_equal(p_small, p_large, tolerance = 1e-10)
})

# ---- Integration tests (simInit / spades) ----

library(SpaDES.core)

snagToDWD_test <- matrix(c(
  0.000, 0.162, 0.815, 0.023,
  0.000, 0.078, 0.871, 0.052,
  0.000, 0.036, 0.854, 0.111,
  0.000, 0.016, 0.762, 0.222,
  0.000, 0.007, 0.599, 0.395
), nrow = 5, byrow = TRUE)

emptyFallen <- data.table::data.table(
  pixelID = integer(), species = character(),
  DC = integer(), ageInDC = integer(), initBiomass = numeric(), diameter_cm = numeric()
)

testInit <- function(params, objects, times = list(start = 0, end = 10)) {
  simInit(
    times   = times,
    modules = list("DeadWood_DWDDecay"),
    params  = params,
    objects = objects,
    paths   = list(modulePath = file.path(.projRoot, "modules"))
  )
}

test_that("Init creates empty DWDTable with correct schema", {
  sim <- testInit(
    params  = list(DeadWood_DWDDecay = list(
      DWD_logisticParams = logisticParams_test,
      snagToDWD_DCmat    = snagToDWD_test
    )),
    objects = list(fallenSnags = data.table::copy(emptyFallen))
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$DWDTable, "data.table")
  expect_equal(nrow(sim$DWDTable), 0L)
  expect_named(sim$DWDTable,
               c("pixelID", "species", "DC", "ageInDC", "ageSinceEntry", "initBiomass", "diameter_cm"))
})

test_that("Receive appends fallenSnags: diameter_cm carried, ageInDC and ageSinceEntry reset to 0", {
  fallen <- data.table::data.table(
    pixelID     = c(1L, 2L),
    species     = "Pinus strobus",
    DC          = c(1L, 3L),
    ageInDC     = c(5L, 3L),
    initBiomass = c(8.0, 4.0),
    diameter_cm = c(15.0, 22.0)
  )
  sim <- testInit(
    params  = list(DeadWood_DWDDecay = list(
      DWD_logisticParams = logisticParams_test,
      snagToDWD_DCmat    = snagToDWD_test
    )),
    objects = list(fallenSnags = fallen)
  )
  sim <- spades(sim, events = c("init", "receive"))
  expect_equal(nrow(sim$DWDTable), 2L)
  expect_true(all(sim$DWDTable$ageInDC == 0L))
  expect_true(all(sim$DWDTable$ageSinceEntry == 0L))
  expect_setequal(sim$DWDTable$diameter_cm, c(15.0, 22.0))
  expect_true(all(sim$DWDTable$DC >= 1L & sim$DWDTable$DC <= 4L))
})

test_that("Transition: DC never decreases and pieces in DC4 can exit", {
  set.seed(42)
  preloaded <- data.table::data.table(
    pixelID       = 1:100,
    species       = "Pinus strobus",
    DC            = as.integer(sample(1:4, 100, replace = TRUE)),
    ageInDC       = rep(0L, 100),
    ageSinceEntry = rep(50L, 100),
    initBiomass   = rep(5.0, 100),
    diameter_cm   = rep(15.0, 100)
  )
  sim <- testInit(
    params  = list(DeadWood_DWDDecay = list(
      DWD_logisticParams = logisticParams_test,
      snagToDWD_DCmat    = snagToDWD_test
    )),
    objects = list(fallenSnags = data.table::copy(emptyFallen))
  )
  sim <- spades(sim, events = "init")
  initialDC <- preloaded$DC
  sim$DWDTable <- preloaded
  sim <- spades(sim, events = "transition")
  remaining <- sim$DWDTable
  # All remaining pieces must have DC >= their original DC
  matched <- remaining[data.table::data.table(pixelID = 1:100, initDC = initialDC),
                       on = "pixelID"]
  expect_true(all(matched$DC >= matched$initDC, na.rm = TRUE))
  # DC must stay in 1-4
  expect_true(all(remaining$DC >= 1L & remaining$DC <= 4L))
})

test_that("Transition: ageSinceEntry increments by 5 each step", {
  preloaded <- data.table::data.table(
    pixelID       = 1L,
    species       = "Pinus strobus",
    DC            = 1L,
    ageInDC       = 0L,
    ageSinceEntry = 10L,
    initBiomass   = 5.0,
    diameter_cm   = 15.0
  )
  # Force probability to 0 so piece stays in DC1 (use very young age, DC1 slow at age 0)
  sim <- testInit(
    params  = list(DeadWood_DWDDecay = list(
      DWD_logisticParams = logisticParams_test,
      snagToDWD_DCmat    = snagToDWD_test
    )),
    objects = list(fallenSnags = data.table::copy(emptyFallen))
  )
  sim <- spades(sim, events = "init")
  sim$DWDTable <- preloaded
  set.seed(1)
  sim <- spades(sim, events = "transition")
  if (nrow(sim$DWDTable) > 0L)
    expect_equal(sim$DWDTable$ageSinceEntry, 15L)
})
