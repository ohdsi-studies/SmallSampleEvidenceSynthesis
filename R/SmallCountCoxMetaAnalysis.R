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

# profileCoxLikelihood <- function(population) {
#   cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), 
#                                             data = population, 
#                                             modelType = "cox")
#   fit <- Cyclops::fitCyclopsModel(cyclopsData)
#   mode <- coef(fit)
#   ci95 <- confint(fit, 1, level = .95)
#   ci99 <- confint(fit, 1, level = .99)
#   if (is.na(ci99[2])) {
#     ci99[2] <- mode - (ci99[3] - mode) * 2
#   }
#   if (is.na(ci99[3])) {
#     ci99[3] <- mode + (mode - ci99[2]) * 2
#   }
#   beta <- seq(from = ci99[2], to = ci99[3], length.out = 100)
#   ll <- rep(0, 100)
#   for (i in 1:100) {
#     # Set starting coefficient, then tell Cyclops that it is fixed:
#     temp <- Cyclops::fitCyclopsModel(cyclopsData, startingCoefficients = beta[i], fixedCoefficients = 1)
#     ll[i] <- temp$log_likelihood
#   }
#   result <- data.frame(beta = beta, ll = ll)
#   attr(result, "mode") <- mode
#   attr(result, "se") <- (ci95[3] - ci95[2]) / (2*qnorm(0.975))
#   return(result)
# }

combiFlexFun <- function(x, fits) {
  if (is.data.frame(fits)) {
    ll <- lapply(split(fits, 1:nrow(fits)), function(fit) flexFun(x, fit$mu, fit$sigma, fit$gamma))
    ll <- do.call(cbind, ll)
    ll <- apply(ll, 1, sum)
    return(ll)
  } else {
    
  }
}

profileCombiFlexFun <- function(fits, alpha = 0.05) {
  # fit <- nlm(function(x) {y <- -combiFlexFun(x, fits); print(y); return(y)}, 0)
  # fit <- optim(0, function(x) {y <- -combiFlexFun(x, fits); print(y); return(y)}, method = "Brent", lower = -10, upper = 10)
  fit <- suppressWarnings(optim(0, function(x) -combiFlexFun(x, fits)))
  # logRr <- fit$estimate           
  # threshold <- -fit$minimum - qchisq(1 - alpha, df = 1) / 2
  logRr <- fit$par           
  threshold <- -fit$value - qchisq(1 - alpha, df = 1) / 2
  
  
  precision = 1e-07
  
  # Binary search for upper bound
  L <- logRr
  H <- 10
  while (H >= L) {
    M <- L + (H - L)/2
    llM <- combiFlexFun(M, fits)
    metric <- threshold - llM
    # writeLines(paste('M =', M, 'Metric = ',metric))
    if (metric > precision) {
      H <- M
    } else if (-metric > precision) {
      L <- M
    } else {
      ub <- M
     # print(M)
      break
    }
    if (M == logRr || M == 10)
      print("Error")
  }
  
  # Binary search for lower bound  
  L <- -10
  H <- logRr
  while (H >= L) {
    M <- L + (H - L)/2
    llM <- combiFlexFun(M, fits)
    metric <- threshold - llM
    # writeLines(paste('M =', M, 'Metric = ',metric))
    if (metric > precision) {
      L <- M
    } else if (-metric > precision) {
      H <- M
    } else {
      lb <- M
      # print(M)
      break
    }
    if (M == -10 || M == logRr)
      print("Error")
  }
  result <- data.frame(hr = exp(logRr),
                       lb = exp(lb),
                       ub = exp(ub),
                       logRr = logRr,
                       seLogRr = (ub - lb)/(2*qnorm(0.975)))
  return(result)
}

flexFun <- function(x, mu, sigma, gamma) {
  # Cubic:
  # return((gamma*(x - mu))^3 + ((-(x - mu)^2)/(2*sigma^2)))
  
  # Logistic:
  # return((1/(1 + exp(gamma*(x - mu)))) * ((-(x - mu)^2)/(sigma^2)))
  
  # Skew normal:
  # return(dsn(x, xi = mu, omega = sigma, alpha = gamma, log = T))
  
  # Exp:
  return(((exp(gamma*(x - mu)))) * ((-(x - mu)^2)/(2*sigma^2)))
}

profileCombiGrids <- function(grids, alpha = 0.05) {
  grid <- apply(grids, 2, sum)
  maxIdx <- which(grid == max(grid))[1]
  logRr <- as.numeric(names(grid)[maxIdx])
  threshold <- grid[maxIdx] - qchisq(1 - alpha, df = 1) / 2
  lbIdx <- min(which(grid[1:maxIdx] > threshold))
  if (lbIdx == 1) {
    warning("Lower bound of confidence interval out of range")
  }
  lb <- as.numeric(names(grid)[lbIdx])
  ubIdx <- maxIdx + max(which(grid[(maxIdx + 1):length(grid)] > threshold))
  if (lbIdx == length(grid)) {
    warning("Upper bound of confidence interval out of range")
  }
  ub <- as.numeric(names(grid)[ubIdx])
  result <- data.frame(hr = exp(logRr),
                       lb = exp(lb),
                       ub = exp(ub),
                       logRr = logRr,
                       seLogRr = (ub - lb)/(2*qnorm(0.975)))
  return(result)
}

fitFlexFun <- function(beta, ll) {
  # Scale to standard (translate in log space so intersects at (0,0)):
    
  sumSquares <- function(p, maxAnchor = TRUE) {
    approxLl <- flexFun(beta, p[1], p[2], p[3])
    if (maxAnchor) {
      approxLl <- approxLl - max(approxLl)
    } else {
      approxLl <- approxLl - approxLl[idx]
    }
    result <- sum((approxLl - ll)^2)
    # print(paste(paste(p, collapse = ","), result))
    return(result)
  }
  
  ll <- ll - max(ll)
  fit <- tryCatch({
    suppressWarnings(nlm(sumSquares, c(0,1,0), maxAnchor = TRUE))
  }, error = function(e){
    list(estimate = c(0,0,0), minimum = Inf)
  })
  result <- data.frame(mu = fit$estimate[1],
                       sigma  = fit$estimate[2],
                       gamma = fit$estimate[3])
  minimum <- fit$minimum
  
  fit <- tryCatch({
    suppressWarnings(optim(c(0,1,0), sumSquares, maxAnchor = TRUE))
  }, error = function(e){
    list(par = c(0,0,0), value = Inf)
  })
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1],
                         sigma  = fit$par[2],
                         gamma = fit$par[3])
    minimum <- fit$value
  }
  
  idx <- which(abs(beta) == min(abs(beta)))
  ll <- ll - ll[idx]
  fit <- tryCatch({
    suppressWarnings(nlm(sumSquares, c(0,1,0), maxAnchor = FALSE))
  }, error = function(e){
    list(estimate = c(0,0,0), minimum = Inf)
  })
  if (fit$minimum < minimum) {
    result <- data.frame(mu = fit$estimate[1],
                         sigma  = fit$estimate[2],
                         gamma = fit$estimate[3])
    minimum <- fit$minimum
  }
  fit <- tryCatch({
    suppressWarnings(optim(c(0,1,0), sumSquares, maxAnchor = FALSE))
  }, error = function(e){
    list(par = c(0,0,0), value = Inf)
  })
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1],
                         sigma  = fit$par[2],
                         gamma = fit$par[3])
    minimum <- fit$value
  }
  threshold <- 0.05
  if (minimum / length(beta) > threshold) {
    warning("Mean squared error greater than ", threshold, ". Probably bad fit.")
  }
  result$sigma <- abs(result$sigma)
  return(result)
}

# Remove this
plotLikelihoodFit <- function(beta, ll) {
  pseudoCoxFit <- result
  # pseudoCoxFit <- fitFlexFun(beta, ll)
  pseudoCoxLl <- flexFun(x = beta, 
                         mu = pseudoCoxFit$mu,
                         sigma = pseudoCoxFit$sigma,
                         gamma = pseudoCoxFit$gamma)
  coxLl <- ll
  coxLl <- coxLl - max(coxLl)
  pseudoCoxLl <- pseudoCoxLl - max(pseudoCoxLl)
  plotData <- rbind(data.frame(beta = beta,
                               ll = coxLl,
                               label = "Cox"),
                    data.frame(beta = beta,
                               ll = pseudoCoxLl,
                               label = "Pseudo Cox"))
  ggplot2::ggplot(plotData, ggplot2::aes(x = beta, y = ll, group = label, color = label)) +
    ggplot2::geom_line(size = 1)
}