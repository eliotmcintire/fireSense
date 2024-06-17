# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense",
  description = "A landscape fire model, sensitive to environmental changes (e.g.
                 weather and land-cover).",
  keywords = c("fire", "percolation", "environmental control", "feedback",
               "weather", "vegetation", "land-cover"),
  authors = c(
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = numeric_version("1.0.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "fireSense.Rmd"), ## same file
  reqdPkgs = list("data.table", "ggplot2", "ggspatial", "raster", "terra"),
  parameters = rbind(
    defineParameter("plotIgnitions", "logical", FALSE, NA, NA,
                    "whether to plot ignitions, escapes, and burns"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "optional. When to start plotting."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "optional. Interval between plot events."),
    defineParameter(".runInitialTime", "numeric", start(sim), NA, NA,
                    "time to simulate initial fire"),
    defineParameter(".runInterval", "numeric", NA, NA, NA,
                    paste("optional. Interval between two runs of this module,",
                          "expressed in units of simulation time. By default, 1 year.")),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "optional. When to start saving output to a file."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "optional. Interval between save events."),
    defineParameter("whichModulesToPrepare", "character",
                    default = c("fireSense_SpreadPredict", "fireSense_IgnitionPredict", "fireSense_EscapePredict"),
                    NA, NA,
                    paste("Which fireSense predict modules to prep? Defaults to all 3.",
                          "Must include `fireSense_IgnitionPredict`."))
  ),
  inputObjects = rbind(
    expectsInput("fireSense_IgnitionPredicted", "data.frame",
                 "A SpatRaster of ignition probabilities."),
    expectsInput("fireSense_EscapePredicted", "SpatRaster",
                 "A SpatRaster of escape probabilities."),
    expectsInput("fireSense_SpreadPredicted", "SpatRaster",
                 "A SpatRaster of spread probabilities.")
  ),
  outputObjects = rbind(
    createsOutput("burnDT", "data.table",
                  "Data table with pixel IDs of most recent burn."),
    createsOutput("burnMap", "SpatRaster",
                  "A raster of cumulative burns"),
    createsOutput("burnSummary", "data.table",
                  "Describes details of all burned pixels."),
    createsOutput("rstAnnualBurnID", "SpatRaster",
                  "annual raster whose values distinguish individual fires"),
    createsOutput("rstCurrentBurn", "SpatRaster",
                  "A binary raster with 1 values representing burned pixels.")
  )
))

## event types
#   - type `init` is required for initialization

doEvent.fireSense = function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {
      ## bail early if there's a problem with ignition, escape, or spread rasters
      if (!is.null(sim$fireSense_IgnitionPredicted))
        stopifnot(length(na.omit(sim$fireSense_IgnitionPredicted[])) > 0)

      if (!is.null(sim$fireSense_EscapePredicted))
        stopifnot(length(na.omit(sim$fireSense_EscapePredicted[])) > 0)

      if (!is.null(sim$fireSense_SpreadPredicted))
        stopifnot(length(na.omit(sim$fireSense_SpreadPredicted[])) > 0)

      ## trying to avoid the raster warning no non-missing arguments to max
      sim$burnMap <- rast(sim$rasterToMatch)
      burnVals <- as.vector(sim$flammableRTM)
      sim$burnMap[burnVals == 0] <- NA
      sim$burnMap[burnVals == 1] <- 0
      rm(burnVals)
      gc()

      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "burn", eventPriority = 5.13)

      if (!is.na(P(sim)$.plotInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, moduleName, "plot", eventPriority = .last())
    },
    burn = {
      sim <- burn(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "burn", eventPriority = 5.13)
    },
    plot = {
      sim <- plot(sim)

      if (!is.na(P(sim)$.plotInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, moduleName, "plot", eventPriority = .last())
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  invisible(sim)
}

burn <- function(sim) {
  moduleName <- current(sim)$moduleName

  ## Ignite
  notNA <- which(!is.na(sim$fireSense_IgnitionPredicted[]))
  ignitionProbs <- sim$fireSense_IgnitionPredicted[notNA]

  ignited <- notNA[which(
    rbinom(n = length(ignitionProbs),
           size = 1,
           prob = pmin(ignitionProbs, 1)
    ) > 0
  )]

  ignited <- sample(ignited) # Randomize order

  rm(ignitionProbs)
  if (!is.na(P(sim)$.plotInitialTime)) {
    mod$ignitions <- ignited
  }

  if (length(ignited) > 0L) {
    if ("fireSense_EscapePredict" %in% P(sim)$whichModulesToPrepare) {
      ## Escape
      adjacent <- SpaDES.tools::adj(
        x = sim$fireSense_EscapePredicted,
        cells = ignited,
        directions = 8,
        returnDT = TRUE
      )

      if (is.matrix(adjacent))
        adjacent <- as.data.table(adjacent)

      from <- unique(adjacent, by = "from")
      pix <- from$from
      from[, `:=`(probEscape = as.vector(sim$fireSense_EscapePredicted)[pix])]

      # Update probEscape to get p0
      from <- from[!is.na(probEscape),]
      p0 <- with(
        data = from[adjacent, on = "from"][, probEscape := (1 - (1 - probEscape)^(1 / .N)), by = "from"],
        expr = {
          p0 <- sim$fireSense_EscapePredicted
          p0[to] <- probEscape
          p0
        }
      )

      mod$spreadState <- SpaDES.tools::spread2(
        landscape = sim$fireSense_EscapePredicted,
        start = ignited,
        iterations = 1,
        spreadProb = p0,
        directions = 8L,
        asRaster = FALSE
      )

      if (!is.na(P(sim)$.plotInitialTime)) {
        mod$escapes <- unique(mod$spreadState[state == "activeSource", ]$initialPixels)
      }
    }
    if ("fireSense_SpreadPredict" %in% P(sim)$whichModulesToPrepare) {
      ## Spread
      # Note: if none of the cells are active SpaDES.tools::spread2() returns spreadState unchanged

      mod$spreadState <- SpaDES.tools::spread2(
        landscape = sim$fireSense_SpreadPredicted,
        spreadProb = sim$fireSense_SpreadPredicted,
        directions = 8L,
        start = mod$spreadState,
        asRaster = FALSE)

      mod$spreadState[ , fire_id := .GRP, by = "initialPixels"] # Add an fire_id column

      sim$rstAnnualBurnID <- rast(sim$fireSense_SpreadPredicted)
      sim$rstCurrentBurn <- rast(sim$fireSense_SpreadPredicted)

      sim$rstAnnualBurnID[mod$spreadState$pixels] <- mod$spreadState$fire_id
      sim$rstCurrentBurn[mod$spreadState$pixels] <- 1
      sim$burnMap[mod$spreadState$pixels] <- sim$burnMap[mod$spreadState$pixels] + 1

      #get fire year, pixels burned, area burned, poly ID of all burned pixels
      # Make burnSummary --> similar to SCFM
      sim$burnDT <- mod$spreadState

      tempDT <- sim$burnDT[, .(.N), by = "initialPixels"]
      tempDT$year <- time(sim)
      tempDT$areaBurnedHa <- tempDT$N * prod(res(sim$fireSense_SpreadPredicted)) * 1e-4
      setnames(tempDT, c("initialPixels"), c("igLoc"))
      sim$burnSummary <- rbind(sim$burnSummary, tempDT)
    }
  }

  invisible(sim)
}

plot <- function(sim) {
  if (P(sim)$plotIgnitions) {
    ## this plot treats escapes and ignitions as points, but burns as rasters
    ## it is impossible to show escapes as a raster, and burns do not plot well as points
    ## this requires some ggplot hacks

    flam <- as.data.frame(as(sim$flammableRTM, "SpatialPixelsDataFrame"))
    names(flam) <- c("value", "x", "y")
    flam[flam$value == 0,]$value <- "unburnable"
    flam[flam$value == 1,]$value <- "burnable"

    escapes <- raster::xyFromCell(sim$flammableRTM, cell = mod$escapes) %>%
      as.data.table(.)
    ignitions <- raster::xyFromCell(sim$flammableRTM, cell = mod$ignitions) %>%
      as.data.table(.)

    ## there should be an easier anti-join in data.table
    both <- rbind(escapes, ignitions)
    both <- both[, .N, .(x, y)] #want only xy of points that did not escape
    ignitions <- both[N == 1, .(x, y)]

    ignitions <- SpatialPointsDataFrame(coords = ignitions,
                                        proj4string = crs(sim$flammableRTM),
                                        data = data.frame(stat = as.factor(rep(x = "ignited",
                                                                               times = nrow(ignitions))))
    )

    escapes <- SpatialPointsDataFrame(coords = escapes,
                                      proj4string = crs(sim$flammableRTM),
                                      data = data.frame(stat = as.factor(rep(x = "escaped",
                                                                             times = nrow(escapes))))
    )

    ignitions <- ggspatial::df_spatial(ignitions)
    escapes <- ggspatial::df_spatial(escapes)
    burns <- as.data.frame(as(sim$rstCurrentBurn, "SpatialPixelsDataFrame"))
    names(burns) <- c("value", "x", "y")
    burns$value <- "burned"
    burns$value <- factor(burns$value, levels = c("ignited", "escaped", "burned"))
    g <- ggplot() +
      geom_raster(data = flam,
                  aes(x = x, y = y, fill = value),
                  show.legend = TRUE) +
      geom_raster(data = burns,
                  aes(x = x, y = y, fill = value),
                  show.legend = FALSE) +
      geom_point(data = ignitions,
                 aes(x = x, y = y),
                 color = "#EFFD5F",
                 show.legend = FALSE,
                 cex = 0.7) +
      geom_point(data = escapes,
                 aes(x = x, y = y),
                 color = "#EC9706",
                 show.legend = FALSE,
                 cex = 0.7) +
      theme_minimal() +
      scale_fill_manual(name = "fire status",
                        values = c("burnable" = "#028A0F", #green = burnable
                                   "burned" = "#D0312D",  #red = burned
                                   "escaped" = "#EC9706", #orange = escaped
                                   "ignited" = "#EFFD5F", #yellow = #ignited
                                   "unburnable" = "#C5C6D0"), #grey = unburnable
                        drop = FALSE)
    #there is a warning about geom_tile, but it can't be used with geom_point
    ggsave(plot = g, filename = paste0('firePlotGG', time(sim), ".png"),
           device = "png", path = file.path(outputPath(sim), "figures"))

    mod$ignitions <- NULL
    mod$escapes <- NULL
  }

  invisible(sim)
}
