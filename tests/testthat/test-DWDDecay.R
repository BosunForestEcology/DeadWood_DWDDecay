library(SpaDES.core)
library(data.table)
library(testthat)

DWDTransMat_test <- matrix(
  c(0.55, 0.35, 0.06, 0.00, 0.00,
    0.00, 0.50, 0.37, 0.08, 0.00,
    0.00, 0.00, 0.48, 0.38, 0.08,
    0.00, 0.00, 0.00, 0.50, 0.36,
    0.00, 0.00, 0.00, 0.00, 0.72),
  nrow = 5, byrow = TRUE
)
snagToDWD_DCmap_test <- c(DC1 = 1L, DC2 = 2L, DC3 = 2L, DC4 = 3L, DC5 = 4L)
DWD_lossProb_test    <- c(DC1 = 0.00, DC2 = 0.05, DC3 = 0.06, DC4 = 0.14, DC5 = 0.28)

emptySnagTable <- data.table(
  pixelID = integer(), species = character(),
  DC = integer(), ageInDC = integer(), initBiomass = numeric()
)

# testInit is not available in SpaDES.core >= 3.x; use simInit directly.
# This test file lives at: modules/DeadWood_DWDDecay/tests/testthat/
# The project root (containing modules/) is 4 levels up.
.projRoot <- normalizePath(file.path(getwd(), "..", "..", "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(.projRoot, "modules"))) {
  .projRoot <- getwd()
}

testInit <- function(moduleName, params, objects, times = list(start = 0, end = 10)) {
  simInit(
    times   = times,
    modules = list(moduleName),
    params  = params,
    objects = objects,
    paths   = list(modulePath = file.path(.projRoot, "modules"))
  )
}

test_that("DWDDecay init creates empty DWDTable with correct schema", {
  sim <- testInit(
    "DeadWood_DWDDecay",
    params  = list(DWDDecay = list(
      DWDTransMat     = DWDTransMat_test,
      snagToDWD_DCmap = snagToDWD_DCmap_test,
      DWD_lossProb    = DWD_lossProb_test
    )),
    objects = list(fallenSnags = data.table::copy(emptySnagTable))
  )
  sim <- spades(sim, events = "init")
  expect_s3_class(sim$DWDTable, "data.table")
  expect_equal(nrow(sim$DWDTable), 0L)
  expect_named(sim$DWDTable, c("pixelID", "species", "DC", "ageInDC", "initBiomass"))
})

test_that("DWDDecay receive appends fallenSnags with DC remapped", {
  fallen <- data.table(
    pixelID     = c(1L, 2L),
    species     = "Pinus strobus",
    DC          = c(1L, 3L),
    ageInDC     = c(2L, 1L),
    initBiomass = c(8.0, 4.0)
  )
  sim <- testInit(
    "DeadWood_DWDDecay",
    times   = list(start = 0, end = 1),
    params  = list(DWDDecay = list(
      DWDTransMat     = DWDTransMat_test,
      snagToDWD_DCmap = snagToDWD_DCmap_test,
      DWD_lossProb    = DWD_lossProb_test
    )),
    objects = list(fallenSnags = fallen)
  )
  sim <- spades(sim, events = c("init", "receive"))
  expect_equal(nrow(sim$DWDTable), 2L)
  # DC1 snag -> DC1 DWD; DC3 snag -> DC2 DWD (per snagToDWD_DCmap)
  expect_equal(sim$DWDTable[pixelID == 1L, DC], 1L)
  expect_equal(sim$DWDTable[pixelID == 2L, DC], 2L)
  # ageInDC resets to 0 upon entering DWD pool
  expect_true(all(sim$DWDTable$ageInDC == 0L))
})

test_that("DWDDecay annual DC never decreases (upper-triangular matrix)", {
  set.seed(77)
  initial_DCs <- as.integer(sample(1:5, 50, replace = TRUE))
  sim <- testInit(
    "DeadWood_DWDDecay",
    times   = list(start = 0, end = 20),
    params  = list(DWDDecay = list(
      DWDTransMat     = DWDTransMat_test,
      snagToDWD_DCmap = snagToDWD_DCmap_test,
      DWD_lossProb    = rep(0, 5)
    )),
    objects = list(fallenSnags = data.table::copy(emptySnagTable))
  )
  sim <- spades(sim, events = "init")
  sim$DWDTable <- data.table(
    pixelID     = 1:50,
    species     = "Pinus strobus",
    DC          = initial_DCs,
    ageInDC     = rep(0L, 50),
    initBiomass = rep(5.0, 50)
  )
  sim <- spades(sim)
  final_DCs <- sim$DWDTable[order(pixelID), DC]
  expect_true(all(final_DCs >= initial_DCs))
})

test_that("DWDDecay annual with 100% loss probability removes all DC5 records", {
  preloaded <- data.table(
    pixelID     = 1:10,
    species     = "Pinus strobus",
    DC          = rep(5L, 10),
    ageInDC     = rep(0L, 10),
    initBiomass = rep(3.0, 10)
  )
  sim <- testInit(
    "DeadWood_DWDDecay",
    times   = list(start = 0, end = 1),
    params  = list(DWDDecay = list(
      DWDTransMat     = DWDTransMat_test,
      snagToDWD_DCmap = snagToDWD_DCmap_test,
      DWD_lossProb    = c(0, 0, 0, 0, 1)
    )),
    objects = list(fallenSnags = data.table::copy(emptySnagTable))
  )
  sim <- spades(sim, events = "init")
  sim$DWDTable <- preloaded
  set.seed(3)
  sim <- spades(sim, events = "annual")
  expect_equal(nrow(sim$DWDTable), 0L)
})
