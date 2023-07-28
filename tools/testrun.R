
x <- snglmm(Sepal.Length ~ Sepal.Width + (1 | Species), iris, inits = list(b = c(1, 1)))
summary(x)
report(x)

# y <- lme4::lmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
#                 REML = F)
# summary(y)
#
# z <- glmmTMB::glmmTMB(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
#                       REML = F)
# summary(z)


# spatial -----------------------------------------------------------------

set.seed(1)
coord <- cbind(rnorm(150), rnorm(150))
D <- dist(coord)
D <- data.matrix(D)
R <- exp(-D)
u <- MASS::mvrnorm(n = 1, rep(0, ncol(D)), Sigma = R)
iris$y <- iris$Sepal.Length + u

(x <- snglmm(y ~ Sepal.Width,
             iris,
             spatial = "exp",
             D = D))

x0 <- runif(5, min(iris$Sepal.Width), max(iris$Sepal.Width))
X0 <- model.matrix(~x0)

newcoord <- cbind(rnorm(5), rnorm(5))

cD <- sapply(seq_len(nrow(newcoord)), function(x) {
  d0 <- data.matrix(dist(rbind(coord, newcoord[x, ])))
  cout <- d0[-nrow(d0), ncol(d0)]
  return(cout)
})

cW <- matrix(1, dim(cD)[1], dim(cD)[2])

kriging(x, newdata = X0, cD = cD[-1,], cW = cW)

# pos <- glmmTMB::numFactor(coord)
# iris$g <- 1
# iris$pos <- pos
# m <- glmmTMB::glmmTMB(y ~ Sepal.Width + exp(0 + pos|g), iris, dispformula = ~0)
#
