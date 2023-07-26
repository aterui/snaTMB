x <- snaTMB(Sepal.Length ~ Sepal.Width + (1 + Sepal.Width| Species), iris, inits = list(b = c(1, 1)))
summary(x)

y <- lme4::lmer(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                REML = F)
summary(y)

z <- glmmTMB::glmmTMB(Sepal.Length ~ Sepal.Width + (1 | Species), iris,
                      REML = F)
summary(z)


# spatial -----------------------------------------------------------------


set.seed(1)
coord <- cbind(rnorm(150), rnorm(150))
D <- dist(coord)
D <- data.matrix(D)
R <- exp(-D)
u <- MASS::mvrnorm(n = 1, rep(0, ncol(D)), Sigma = R)
iris$y <- iris$Sepal.Length + u

x <- snaTMB(y ~ Sepal.Width,
            iris,
            D = D)
x
summary(x)


pos <- glmmTMB::numFactor(coord)
iris$g <- 1
iris$pos <- pos
fit <- glmmTMB::glmmTMB(y ~ Sepal.Width + exp(0 + pos|g), iris, dispformula = ~0)

