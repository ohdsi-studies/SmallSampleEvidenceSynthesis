# Copyright 2019 Observational Health Data Sciences and Informatics
#
# This file is part of SmallSampleEvidenceSynthesis
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

testRealWorldExample <- function(localFolder = "C:/home/Research/SmallCountMetaAnalysis/Example") {
  # library(survival)
  examples <- read.csv(file.path(localFolder, "Examples.csv"), stringsAsFactors = FALSE)

  databases <- unique(examples$database)
  examples <- unique(examples[, c("exampleId", "targetName", "comparatorName", "outcomeName")])
  
  loadStrataFile <- function(database, exampleId) {
    return(readRDS(file.path(localFolder, sprintf("StratPop_%s_%s.rds", database, exampleId))))
  }
  
  loadStrataFiles <- function(exampleId) {
    populations <- lapply(databases, loadStrataFile, exampleId = exampleId)
    names(populations) <- databases
    return(populations)
  }
  allPopulations <- lapply(examples$exampleId, loadStrataFiles) 
  names(allPopulations) <- examples$exampleId
  
  x <- seq(log(0.1), log(10), length.out = 100)
  
  # Generate likelihood data:
  likelihoodData <- data.frame()
  maEstimates <- data.frame()
  for (i in 1:nrow(examples)) {
    writeLines(paste("Processing example", examples$exampleId[i]))
    exampleId <- examples$exampleId[i]
    populations <- allPopulations[[i]]
    estimates <- data.frame()
    flexFits <- data.frame()
    grids <- data.frame()
    for (database in databases) {
      writeLines(paste("Processing database", database))
      population <- populations[[database]]
      data1 <- data.frame(exampleId = exampleId,
                         database = database,
                         x = x,
                         y = getLikelihoodData(population, x),
                         type = "Log likelihood")
      data1$y <- data1$y - max(data1$y)
      likelihoodData <- rbind(likelihoodData, data1)  
      
      model <- fitIndividualModel(population, useCyclops = FALSE)    
      if (!is.na(model$seLogRr) && !is.infinite(model$seLogRr) && model$seLogRr < 10) {
        
        
        data2 <- data.frame(exampleId = exampleId,
                            database = database,
                            x = x,
                            y = dnorm(x, model$logRr, model$seLogRr, log = TRUE),
                            type = "Normal approximation")
        
        data2$y <- data2$y - max(data2$y)
        # data2$y <- data2$y / normFactor
        likelihoodData <- rbind(likelihoodData, data2)
        estimates <- rbind(estimates, model)
      }
      # wideX <- seq(from = -10, to = 10, by = 0.2)
      # wideY <- getLikelihoodData(population, wideX)
      # fit <- fitFlexFun(wideX, wideY)
      fit <- fitFlexFun(x, data1$y)
      data3 <- data.frame(exampleId = exampleId,
                          database = database,
                          x = x,
                          y = flexFun(x, fit$mu, fit$sigma, fit$gamma),
                          type = "Flex fit")
      data3$y <- data3$y - max(data3$y)
      likelihoodData <- rbind(likelihoodData, data3) 
      flexFits <- rbind(flexFits, fit)
      
      grid <- getLikelihoodData(population, log(seq(from = 0.1, to = 10, by = 0.01)))
      grids <- rbind(grids, grid)
    }
    writeLines("Processing combined estimates")
    # Summary across DBs
    pooledPop <- combinePopulations(populations)
    data1 <- data.frame(exampleId = exampleId,
                        database = "Combined",
                        x = x,
                        y = getLikelihoodData(pooledPop, x),
                        type = "Log likelihood")
    data1$y <- data1$y - max(data1$y)
    likelihoodData <- rbind(likelihoodData, data1)
    fit <- fitIndividualModel(pooledPop, useCyclops = TRUE)
    maEstimates <- rbind(maEstimates, data.frame(exampleId = examples$exampleId[i],
                                                 type = "Pooled Cox",
                                                 hr = fit$rr,
                                                 lb = fit$ci95Lb,
                                                 ub = fit$ci95Ub))
    
    meta <- meta::metagen(TE = estimates$logRr, seTE = estimates$seLogRr)
    # model <- summary(meta)$random
    model <- summary(meta)$fixed
    data2 <- data.frame(exampleId = exampleId,
                        database = "Combined",
                        x = x,
                        y = dnorm(x, model$TE, model$seTE, log = TRUE),
                        type = "Normal approximation")
    
    data2$y <- data2$y - max(data2$y)
    likelihoodData <- rbind(likelihoodData, data2)
    maEstimates <- rbind(maEstimates, data.frame(exampleId = examples$exampleId[i],
                                                 type = "Traditional meta-analysis",
                                                 hr = exp(model$TE),
                                                 lb = exp(model$lower),
                                                 ub = exp(model$upper)))
    
    data3 <- data.frame(exampleId = exampleId,
                        database = "Combined",
                        x = x,
                        y = combiFlexFun(x, flexFits),
                        type = "Flex fit")
    
    data3$y <- data3$y - max(data3$y)
    likelihoodData <- rbind(likelihoodData, data3)
    
    model <- profileCombiFlexFun(flexFits)
    maEstimates <- rbind(maEstimates, data.frame(exampleId = examples$exampleId[i],
                                                 type = "Flex fit meta-analysis",
                                                 hr = model$hr,
                                                 lb = model$lb,
                                                 ub = model$ub))
    
    # Save for David:
    write.csv(estimates, sprintf("c:/temp/normal_example_%s.csv", examples$exampleId[i]), row.names = FALSE)
    
    
    write.csv(flexFits, sprintf("c:/temp/parametric_example_%s.csv", examples$exampleId[i]), row.names = FALSE)
    
    names(grids) <- log(seq(from = 0.1, to = 10, by = 0.01))
    write.csv(grids, sprintf("c:/temp/grids_example_%s.csv", examples$exampleId[i]), row.names = FALSE)
    
    
  }
  likelihoodData$database <- factor(likelihoodData$database, levels = c(databases[order(databases)], "Combined"))
  
  bgColor <- data.frame(xmin = min(likelihoodData$x), 
                        xmax = max(likelihoodData$x), 
                        ymin = -20, #min(likelihoodData$y), 
                        ymax = 0,
                        x = 0,
                        y = 0,
                        exampleId = examples$exampleId,
                        database = "Combined")
  
  maEstimates$label <- sprintf("%s: %0.2f (%0.2f-%0.2f)", maEstimates$type, maEstimates$hr, maEstimates$lb, maEstimates$ub)
  # labels <- aggregate(label ~ exampleId, maEstimates, paste) 
  maEstimates$database <- "Combined"
  maEstimates$x <- log(9)
  maEstimates$y <- rep(3:5 * -3.7, nrow(examples))
  maEstimates$database <- factor(maEstimates$database, levels = c(databases[order(databases)], "Combined"))
  
  exampleLabels <- sprintf("%s vs. %s\nfor %s", examples$targetName, examples$comparatorName, examples$outcomeName)
  names(exampleLabels) <- examples$exampleId
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  ggplot2::ggplot(likelihoodData, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), data = bgColor, fill = rgb(0.95, 0.95, 1.0)) +
    ggplot2::geom_vline(xintercept = log(breaks), color = "#AAAAAA", lty = 1, size = 0.5) +
    ggplot2::geom_line(ggplot2::aes(group = type, color = type, linetype = type), size = 1, alpha = 0.7) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
    ggplot2::geom_label(ggplot2::aes(label = label), hjust = 1, data = maEstimates) +
    ggplot2::scale_y_continuous("Log likelihood", limits = c(-20, 0)) +
    ggplot2::scale_x_continuous("Hazard Ratio", limits = log(c(0.1, 10)), breaks = log(breaks), labels = breaks) +
    ggplot2::facet_grid(database~exampleId, labeller = ggplot2::labeller(exampleId = exampleLabels)) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), 
                   panel.background = ggplot2::element_blank(), 
                   panel.grid.major = ggplot2::element_blank(), 
                   legend.title = ggplot2::element_blank(), 
                   axis.ticks = ggplot2::element_blank(), 
                   axis.text.y = ggplot2::element_blank(), 
                   strip.background = ggplot2::element_blank(),
                   legend.position = "top")
  
  ggplot2::ggsave("c:/temp/plotLogistic.png", width = 11, height = 9)
  
  # temp <- likelihoodData
  likelihoodData <- temp
  
  data1 <- aggregate(y ~ x + exampleId, data = likelihoodData[likelihoodData$database != "Combined" & likelihoodData$type == "Flex fit", ], sum)
  data1$y <- data1$y - max(data1$y)
  data1$database <- "Combined"
  data1$type = "Flex fit"
  
  likelihoodData <- rbind(likelihoodData, data1)
  
  likelihoodData[likelihoodData$database != "MDCD" & likelihoodData$type == "Normal approximation" & likelihoodData$exampleId == 2, ]
  
  data <- likelihoodData[likelihoodData$database == "CCAE" & likelihoodData$type == "Log likelihood" & likelihoodData$exampleId == 3, ]
  data <- likelihoodData[likelihoodData$database == "Optum" & likelihoodData$type == "Log likelihood" & likelihoodData$exampleId == 2, ]
  data <- likelihoodData[likelihoodData$database == "MDCD" & likelihoodData$type == "Log likelihood" & likelihoodData$exampleId == 2, ]
  data <- likelihoodData[likelihoodData$database == "MDCR" & likelihoodData$type == "Log likelihood" & likelihoodData$exampleId == 2, ]
  
  fitFunction <- function(data) {
    coxProfile <- data.frame(beta = data$x, ll = data$y)
    beta <- data$x
    ll <- data$y
    fit <- fitFlexFun(coxProfile)
    plotLikelihoodFit(data.frame(beta = data$x, ll = data$y), fit)

  }
  
  
  
  
  fittedFunction <- lapply()
  
  
  # Fit individual models, then perform meta-analysis:
  estimates <- lapply(populations, fitIndividualModel)
  estimates <- do.call("rbind", estimates)
  estimates <- estimates[!is.na(estimates$seLogRr) & estimates$seLogRr < 100, ]
  meta <- meta::metagen(TE = estimates$logRr, seTE = estimates$seLogRr)
  # ma <- summary(meta)$random
  ma <- summary(meta)$random
  maEstimate <- data.frame(rr = exp(ma$TE),
                           ci95Lb = ma$lower,
                           ci95Ub = ma$upper,
                           logRr = ma$TE,
                           seLogRr = ma$seTE)
  
  # Pool data, then fit single model:
  pooledPop <- combinePopulations(populations)
  fit <- fitIndividualModel(pooledPop)
  print(maEstimate)
  print(fit)
  
  plotLikelihood(populations[[2]])
}

fitIndividualModel <- function(population, useCyclops = FALSE) {
  if (useCyclops) {
    # Using Cyclops (doesn't really assume normality)
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
    fit <- Cyclops::fitCyclopsModel(cyclopsData)
    mode <- coef(fit)
    ci95 <- confint(fit, 1, level = .95)
    return(data.frame(rr = exp(mode),
                      ci95Lb = exp(ci95[2]),
                      ci95Ub = exp(ci95[3]),
                      logRr = mode,
                      seLogRr = (ci95[3] - ci95[2])/(2 * qnorm(0.975))))
  } else {
    # Using vanilla Cox regression:
    fit <- coxph(Surv(time, y) ~ x + strata(stratumId), population)
    ci95 <- confint(fit)
    return(data.frame(rr = exp(fit$coefficients[1]),
                      ci95Lb = ci95[1],
                      ci95Ub = ci95[2],
                      logRr = fit$coefficients[1],
                      seLogRr = sqrt(fit$var[1,1])))
  }
}

getLikelihoodData <- function(population, x) {
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
  cyclopsData$rowNames <- NULL
  getLikelihood <- function(x) {
    Cyclops::fitCyclopsModel(cyclopsData, startingCoefficients = x, fixedCoefficients = 1)$log_likelihood
  }
  y <- sapply(x, getLikelihood)
  return(y)
}


plotLikelihood <- function(population) {
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
  fit <- Cyclops::fitCyclopsModel(cyclopsData)
  mode <- coef(fit)
  ci95 <- confint(fit, 1, level = .95)
  
  # remove this
  # threshold <- qchisq(0.95, df = 1) / 2
  # temp <- Cyclops::fitCyclopsModel(cyclopsData, startingCoefficients = ci95[3], fixedCoefficients = 1)
  # y <- temp$log_likelihood
  # fit$log_likelihood - y
  ci99 <- c(0,-5, 5)
  mode = 0
  
  # Plot likelihood distribution ----------------------------------------------------
  ci99 <- confint(fit, 1, level = .99)
  if (is.na(ci99[2])) {
    ci99[2] <- mode - (ci99[3] - mode) * 2
  }
  if (is.na(ci99[3])) {
    ci99[3] <- mode + (mode - ci99[2]) * 2
  }
  x <- seq(from = ci99[2], to = ci99[3], length.out = 100)
  y <- rep(0, 100)
  for (i in 1:100) {
    # Set starting coefficient, then tell Cyclops that it is fixed:
    temp <- Cyclops::fitCyclopsModel(cyclopsData, startingCoefficients = x[i], fixedCoefficients = 1)
    y[i] <- temp$log_likelihood
  }
  plot <- ggplot2::ggplot(data.frame(x = x, y = y), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = mode), data = data.frame(mode = mode)) +
    ggplot2::geom_vline(linetype = "dotted", ggplot2::aes(xintercept = ci), data = data.frame(ci = ci95[2:3])) +
    ggplot2::xlab("Beta") +
    ggplot2::ylab("Log likelihood")
  return(plot)
}

combinePopulations <- function(populations) {
  highestId <- 0
  for (i in 1:length(populations)) {
    # Making sure stratum IDs are unique:
    populations[[i]]$stratumId <- populations[[i]]$stratumId + highestId
    highestId <- max(populations[[i]]$stratumId)
  }
  pooledPop <- do.call("rbind", populations) 
  return(pooledPop)
}
