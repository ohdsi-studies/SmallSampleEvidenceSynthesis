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
  
  files <- list.files(localFolder)
  populations <- lapply(files, function(x) readRDS(file.path(localFolder, x)))
  
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
  
  plotLikelihood(populations[[1]])
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

plotLikelihood <- function(population) {
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
  fit <- Cyclops::fitCyclopsModel(cyclopsData)
  mode <- coef(fit)
  ci95 <- confint(fit, 1, level = .95)
  
  
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