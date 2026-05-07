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
    defineParameter("DWDTransMat", "matrix", matrix(0, 5, 5), NA, NA,
                    desc = "5x5 annual DC transition probability matrix for DWD."),
    defineParameter("snagToDWD_DCmap", "integer", 1:5, 1L, 5L,
                    desc = "Maps snag DC at fall to starting DWD DC (length 5)."),
    defineParameter("DWD_lossProb", "numeric", rep(0.05, 5), 0, 1,
                    desc = "Annual complete-loss probability by DC (length 5).")
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
      sim <- scheduleEvent(sim, start(sim) + 1, "DeadWood_DWDDecay", "receive", eventPriority = 2)
      sim <- scheduleEvent(sim, start(sim) + 1, "DeadWood_DWDDecay", "annual",  eventPriority = 3)
    },
    receive = {
      sim <- Receive(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "DeadWood_DWDDecay", "receive", eventPriority = 2)
    },
    annual = {
      sim <- Annual(sim)
      sim <- scheduleEvent(sim, time(sim) + 1, "DeadWood_DWDDecay", "annual",  eventPriority = 3)
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
  if (length(P(sim)$snagToDWD_DCmap) != 5L)
    stop("snagToDWD_DCmap must have length 5.")

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
  incoming[, DC      := P(sim)$snagToDWD_DCmap[DC]]
  incoming[, ageInDC := 0L]
  sim$DWDTable <- data.table::rbindlist(list(sim$DWDTable, incoming))
  return(invisible(sim))
}

Annual <- function(sim) {
  if (nrow(sim$DWDTable) == 0L) return(invisible(sim))

  oldDC <- sim$DWDTable$DC
  sim$DWDTable[, DC := applyTransition(DC, P(sim)$DWDTransMat)]
  sim$DWDTable[, ageInDC := data.table::fifelse(DC == oldDC, ageInDC + 1L, 0L)]

  lossIdx <- sim$DWDTable[, stats::rbinom(.N, 1L, P(sim)$DWD_lossProb[DC]) == 1L]
  sim$DWDTable <- sim$DWDTable[!lossIdx]

  return(invisible(sim))
}
