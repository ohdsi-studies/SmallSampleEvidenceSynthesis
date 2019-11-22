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

profileCoxLikelihood <- function(population) {
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), 
                                            data = population, 
                                            modelType = "cox")
  fit <- Cyclops::fitCyclopsModel(cyclopsData)
  mode <- coef(fit)
  ci95 <- confint(fit, 1, level = .95)
  ci99 <- confint(fit, 1, level = .99)
  if (is.na(ci99[2])) {
    ci99[2] <- mode - (ci99[3] - mode) * 2
  }
  if (is.na(ci99[3])) {
    ci99[3] <- mode + (mode - ci99[2]) * 2
  }
  beta <- seq(from = ci99[2], to = ci99[3], length.out = 100)
  ll <- rep(0, 100)
  for (i in 1:100) {
    # Set starting coefficient, then tell Cyclops that it is fixed:
    temp <- Cyclops::fitCyclopsModel(cyclopsData, startingCoefficients = beta[i], fixedCoefficients = 1)
    ll[i] <- temp$log_likelihood
  }
  result <- data.frame(beta = beta, ll = ll)
  attr(result, "mode") <- mode
  attr(result, "se") <- (ci95[3] - ci95[2]) / (2*qnorm(0.975))
  return(result)
}

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
  fit <- optim(0, function(x) -combiFlexFun(x, fits))
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
  return((gamma*(x - mu))^3 + ((-(x - mu)^2)/sigma))
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
  fit <- nlm(sumSquares, c(0,1,0), maxAnchor = TRUE)
  result <- data.frame(mu = fit$estimate[1],
                       sigma  = fit$estimate[2],
                       gamma = fit$estimate[3])
  minimum <- fit$minimum
  
  fit <- optim(c(0,1,0), sumSquares, maxAnchor = TRUE)
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1],
                         sigma  = fit$par[2],
                         gamma = fit$par[3])
    minimum <- fit$value
  }
  
  idx <- which(abs(beta) == min(abs(beta)))
  ll <- ll - ll[idx]
  fit <- nlm(sumSquares, c(0,1,0), maxAnchor = FALSE)
  if (fit$minimum < minimum) {
    result <- data.frame(mu = fit$estimate[1],
                         sigma  = fit$estimate[2],
                         gamma = fit$estimate[3])
    minimum <- fit$minimum
  }
  fit <- optim(c(0,1,0), sumSquares, maxAnchor = FALSE)
  if (fit$value < minimum) {
    result <- data.frame(mu = fit$par[1],
                         sigma  = fit$par[2],
                         gamma = fit$par[3])
    minimum <- fit$value
  }
  return(result)
}

# Remove this
plotLikelihoodFit <- function(coxProfile, pseudoCoxFit) {
  result <- fitFlexFun(coxProfile$beta, coxProfile$ll)
  pseudoCoxFit = result

  pseudoCoxLl <- flexFun(x = coxProfile$beta, 
                         mu = pseudoCoxFit$mu,
                         sigma = pseudoCoxFit$sigma,
                         gamma = pseudoCoxFit$gamma)
  # mean <- attr(coxProfile, "mode")
  # se <- attr(coxProfile, "se")
  # normalLl <- dnorm(coxProfile$beta, mean, se, log = TRUE) 
  coxLl <- coxProfile$ll
  
  # Scale to standard:
  coxLl <- coxLl - max(coxLl)
  # coxLl <- coxLl / max(coxLl)
  # normalLl <- normalLl - min(normalLl)
  # normalLl <- normalLl / max(normalLl)
  pseudoCoxLl <- pseudoCoxLl - max(pseudoCoxLl)
  # pseudoCoxLl <- pseudoCoxLl / max(pseudoCoxLl)
  
  
  plotData <- rbind(data.frame(beta = coxProfile$beta,
                               ll = coxLl,
                               label = "Cox"),
                    data.frame(beta = coxProfile$beta,
                               ll = pseudoCoxLl,
                               label = "Pseudo Cox"))
                    # data.frame(beta = coxProfile$beta,
                    #            ll = normalLl,
                    #            label = "Normal"))
  ggplot2::ggplot(plotData, ggplot2::aes(x = beta, y = ll, group = label, color = label)) +
    ggplot2::geom_line(size = 1)
}