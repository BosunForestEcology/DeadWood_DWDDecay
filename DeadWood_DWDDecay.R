# spatialExtent field omitted: removed from SpaDES.core API in version >= 3.0
defineModule(sim, list(
  name        = "DeadWood_DWDDecay",
  description = "Manages the downed woody debris pool. Receives fallen snags from DeadWood_snagDecay,
                 maps their decay class to the DWD scale, advances DC annually, and removes
                 fully decomposed records.",
  keywords    = c("dead wood", "DWD", "downed woody debris", "decay class", "Markov"),
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
    defineParameter("DWDTransMat", "matrix",
                    matrix(c(
                      0.55, 0.35, 0.06, 0.00, 0.00,  # from DC1
                      0.00, 0.50, 0.37, 0.08, 0.00,  # from DC2
                      0.00, 0.00, 0.48, 0.38, 0.08,  # from DC3
                      0.00, 0.00, 0.00, 0.50, 0.36,  # from DC4
                      0.00, 0.00, 0.00, 0.00, 0.72   # from DC5
                    ), nrow = 5, byrow = TRUE,
                    dimnames = list(paste0("from_DC", 1:5), paste0("to_DC", 1:5))),
                    NA, NA,
                    desc = "5x5 annual DC transition probability matrix for DWD (Pinus strobus). Replace with values from Table 2 of the integrated snag/DWD decay model paper."),
    defineParameter("snagToDWD_DCmat", "matrix",
                    matrix(c(
                      0.000, 0.162, 0.815, 0.016, 0.007,  # snag DC1
                      0.000, 0.078, 0.871, 0.035, 0.017,  # snag DC2
                      0.000, 0.036, 0.854, 0.074, 0.037,  # snag DC3
                      0.000, 0.016, 0.762, 0.140, 0.082,  # snag DC4
                      0.000, 0.007, 0.599, 0.226, 0.169   # snag DC5
                    ), nrow = 5, byrow = TRUE,
                    dimnames = list(paste0("snagDC", 1:5), paste0("DWDDC", 1:5))),
                    NA, NA,
                    desc = "5x5 probability matrix: rows = snag DC at fall, columns = DWD DC entered. DWD DC1 is never entered (all zeros). Source: Vanderwel et al. 2006 Table 4."),
    defineParameter("DWD_lossProb", "numeric",
                    c(DC1 = 0.00, DC2 = 0.05, DC3 = 0.06, DC4 = 0.14, DC5 = 0.28),
                    0, 1,
                    desc = "Annual complete-loss probability by DC (Pinus strobus). Replace with species-specific values from the integrated snag/DWD decay model paper.")
  ),
  inputObjects = bindrows(
    expectsInput("fallenSnags", "data.table",
                 desc = "Snags that fell in the current timestep from snagDecay.")
  ),
  outputObjects = bindrows(
    createsOutput("DWDTable", "data.table",
                  desc = "Current DWD inventory: pixelID, species, DC, ageInDC, initBiomass.")
  )
))

doEvent.DeadWood_DWDDecay <- function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, start(sim) + 5, "DeadWood_DWDDecay", "receive",     eventPriority = 2)
      sim <- scheduleEvent(sim, start(sim) + 5, "DeadWood_DWDDecay", "transition",  eventPriority = 3)
    },
    receive = {
      sim <- Receive(sim)
      sim <- scheduleEvent(sim, time(sim) + 5, "DeadWood_DWDDecay", "receive",     eventPriority = 2)
    },
    transition = {
      sim <- Transition(sim)
      sim <- scheduleEvent(sim, time(sim) + 5, "DeadWood_DWDDecay", "transition",  eventPriority = 3)
    },
    warning(paste("Undefined event type:", eventType, "in module DWDDecay"))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  if (all(P(sim)$DWDTransMat == 0))
    stop("DWDTransMat is the zero matrix — provide a real transition matrix in params.")
  if (length(P(sim)$DWD_lossProb) != 5L)
    stop("DWD_lossProb must have length 5 (one probability per decay class).")
  if (!is.matrix(P(sim)$snagToDWD_DCmat) || !identical(dim(P(sim)$snagToDWD_DCmat), c(5L, 5L)))
    stop("snagToDWD_DCmat must be a 5x5 matrix (rows = snag DC, columns = DWD DC).")

  sim$DWDTable <- data.table::data.table(
    pixelID     = integer(),
    species     = character(),
    DC          = integer(),
    ageInDC     = integer(),
    initBiomass = numeric()
  )
  return(invisible(sim))
}

Receive <- function(sim) {
  if (is.null(sim$fallenSnags) || nrow(sim$fallenSnags) == 0L) {
    return(invisible(sim))
  }
  incoming <- data.table::copy(sim$fallenSnags)
  mat <- P(sim)$snagToDWD_DCmat
  incoming[, DC      := vapply(DC, function(d) sample(5L, 1L, prob = mat[d, ]), integer(1L))]
  incoming[, ageInDC := 0L]
  sim$DWDTable <- data.table::rbindlist(list(sim$DWDTable, incoming))
  return(invisible(sim))
}

Transition <- function(sim) {
  if (nrow(sim$DWDTable) == 0L) return(invisible(sim))

  oldDC <- sim$DWDTable$DC
  sim$DWDTable[, DC := applyTransition(DC, P(sim)$DWDTransMat)]
  sim$DWDTable[, ageInDC := data.table::fifelse(DC == oldDC, ageInDC + 5L, 0L)]

  lossIdx <- sim$DWDTable[, stats::rbinom(.N, 1L, P(sim)$DWD_lossProb[DC]) == 1L]
  sim$DWDTable <- sim$DWDTable[!lossIdx]

  return(invisible(sim))
}
