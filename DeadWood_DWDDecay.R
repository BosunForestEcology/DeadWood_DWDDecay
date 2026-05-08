# spatialExtent field omitted: removed from SpaDES.core API in version >= 3.0
defineModule(sim, list(
  name        = "DeadWood_DWDDecay",
  description = "Manages the downed woody debris pool. Receives fallen snags from DeadWood_snagDecay,
                 maps their decay class to the DWD scale (DC1-4; DC4 combines original DC4 and DC5),
                 and advances each piece through decay using a time- and diameter-dependent logistic
                 transition model. Pieces exiting DC4 are removed as fully decomposed.",
  keywords    = c("dead wood", "DWD", "downed woody debris", "decay class", "logistic"),
  authors     = structure(list(list(given = "First", family = "Last",
                                    role = c("aut", "cre"),
                                    email = "email@example.com", comment = NULL)),
                           class = "person"),
  childModules = character(0),
  version     = list(DeadWood_DWDDecay = "0.0.1"),
  timeframe   = as.POSIXlt(c(NA, NA)),
  timeunit    = "year",
  citation    = list(),
  documentation = list(),
  reqdPkgs    = list("data.table", "SpaDES.core (>= 3.0.0)"),
  parameters  = bindrows(
    defineParameter("DWD_logisticParams", "data.table",
                    data.table::data.table(
                      a0 = c(-7.295, -9.663,  3.882, -5.584),
                      a1 = c( 1.888,  0.000, -2.986,  0.000),
                      b0 = c( 1.325,  0.934, -0.076,  0.149),
                      b1 = c(-0.239,  0.000,  0.093,  0.000)
                    ),
                    NA, NA,
                    desc = "4-row data.table of logistic transition parameters (a0, a1, b0, b1).
                            Row 1: DC1->DC2+; Row 2: DC2->DC3+; Row 3: DC3->DC4+; Row 4: DC4->out.
                            a = a0 + a1*ln(D), b = b0 + b1*ln(D) where D is piece diameter in cm."),
    defineParameter("snagToDWD_DCmat", "matrix",
                    matrix(c(
                      0.000, 0.162, 0.815, 0.023,  # snag DC1
                      0.000, 0.078, 0.871, 0.052,  # snag DC2
                      0.000, 0.036, 0.854, 0.111,  # snag DC3
                      0.000, 0.016, 0.762, 0.222,  # snag DC4
                      0.000, 0.007, 0.599, 0.395   # snag DC5
                    ), nrow = 5, byrow = TRUE,
                    dimnames = list(paste0("snagDC", 1:5), paste0("DWDDC", 1:4))),
                    NA, NA,
                    desc = "5x4 probability matrix: rows = snag DC at fall, columns = DWD DC entered
                            (DC1-4; DC4 combines original DC4 and DC5). Source: Vanderwel et al. 2006 Table 4.")
  ),
  inputObjects = bindrows(
    expectsInput("fallenSnags", "data.table",
                 desc = "Snags that fell in the current timestep from snagDecay:
                         pixelID, species, DC, ageInDC, initBiomass, diameter_cm.")
  ),
  outputObjects = bindrows(
    createsOutput("DWDTable", "data.table",
                  desc = "Current DWD inventory: pixelID, species, DC (1-4), ageInDC,
                          ageSinceEntry, initBiomass, diameter_cm.")
  )
))

doEvent.DeadWood_DWDDecay <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim) + 5, "DeadWood_DWDDecay", "receive",    eventPriority = 2)
      sim <- scheduleEvent(sim, start(sim) + 5, "DeadWood_DWDDecay", "transition", eventPriority = 3)
    },
    receive = {
      sim <- Receive(sim)
      sim <- scheduleEvent(sim, time(sim) + 5, "DeadWood_DWDDecay", "receive",    eventPriority = 2)
    },
    transition = {
      sim <- Transition(sim)
      sim <- scheduleEvent(sim, time(sim) + 5, "DeadWood_DWDDecay", "transition", eventPriority = 3)
    },
    warning(paste("Undefined event type:", eventType, "in module DWDDecay"))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  if (!data.table::is.data.table(P(sim)$DWD_logisticParams) || nrow(P(sim)$DWD_logisticParams) != 4L)
    stop("DWD_logisticParams must be a data.table with 4 rows (one per DC transition).")
  if (!is.matrix(P(sim)$snagToDWD_DCmat) || !identical(dim(P(sim)$snagToDWD_DCmat), c(5L, 4L)))
    stop("snagToDWD_DCmat must be a 5x4 matrix (rows = snag DC, columns = DWD DC 1-4).")

  sim$DWDTable <- data.table::data.table(
    pixelID       = integer(),
    species       = character(),
    DC            = integer(),
    ageInDC       = integer(),
    ageSinceEntry = integer(),
    initBiomass   = numeric(),
    diameter_cm   = numeric()
  )
  return(invisible(sim))
}

Receive <- function(sim) {
  if (is.null(sim$fallenSnags) || nrow(sim$fallenSnags) == 0L) {
    return(invisible(sim))
  }
  incoming <- data.table::copy(sim$fallenSnags)
  mat <- P(sim)$snagToDWD_DCmat
  incoming[, DC            := vapply(DC, function(d) sample(4L, 1L, prob = mat[d, ]), integer(1L))]
  incoming[, ageInDC       := 0L]
  incoming[, ageSinceEntry := 0L]
  sim$DWDTable <- data.table::rbindlist(list(sim$DWDTable, incoming))
  return(invisible(sim))
}

Transition <- function(sim) {
  if (nrow(sim$DWDTable) == 0L) return(invisible(sim))

  if (anyNA(sim$DWDTable$diameter_cm))
    stop("DWDTable$diameter_cm contains NA. Ensure cohortData includes a non-NA diameter_cm column.")
  if (anyNA(sim$DWDTable$ageSinceEntry))
    stop("DWDTable$ageSinceEntry contains NA.")

  oldDC <- sim$DWDTable$DC

  transProb <- computeDWDTransProb(
    dc_vec          = oldDC,
    age_A_vec       = sim$DWDTable$ageSinceEntry,
    diameter_vec    = sim$DWDTable$diameter_cm,
    logistic_params = P(sim)$DWD_logisticParams
  )

  if (anyNA(transProb))
    stop("NA transition probabilities computed — check DWD_logisticParams and diameter_cm values.")
  didTransition <- stats::rbinom(length(transProb), 1L, transProb) == 1L

  sim$DWDTable[, ageSinceEntry := ageSinceEntry + 5L]
  sim$DWDTable[didTransition, DC := DC + 1L]
  sim$DWDTable[, ageInDC := data.table::fifelse(DC == oldDC, ageInDC + 5L, 0L)]
  sim$DWDTable <- sim$DWDTable[DC <= 4L]

  return(invisible(sim))
}
