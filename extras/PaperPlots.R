library(SmallSampleEvidenceSynthesis)


# Equivalence with normal ------------------------
beta <- -10:10
ll <- dnorm(beta, 0.3, 0.2, log = TRUE)
fit <- fitFlexFun(beta, ll)
fit


# Flex function ------------------------------------
x <- seq(log(0.1), log(10), length.out = 100)
vizData <- data.frame(x = c(x, 
                            x, 
                            x),
                      y = c(flexFun(x, 0, 1, 0),
                            flexFun(x, 0, 1, 0.2),
                            flexFun(x, 0, 1, -1)),
                      label = c(rep("mu = 0; sigma = 1; gamma = 0", length(x)),
                                    rep("mu = 0; sigma = 1; gamma = 0.2", length(x)),
                                        rep("mu = 0; sigma = 1; gamma = -1", length(x))))
ggplot2::ggplot(vizData, ggplot2::aes(x = x, y = y, group = label, color = label)) +
  ggplot2::geom_line() +
  ggplot2::scale_y_continuous(limits = c(-5,0))

