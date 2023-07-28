
# build package -----------------------------------------------------------

usethis::use_mit_license(copyright_holder = "Akira Terui")
usethis::use_roxygen_md()
usethis::use_namespace()
#usethis::use_package_doc()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)


# check syntax ------------------------------------------------------------

lintr::lint_package()

#usethis::use_coverage()
# devtools::build(path='.')
# covr::package_coverage()
# file.remove("mcbrnet_1.2.3.tar.gz")


# build website -----------------------------------------------------------

# Run once to configure package to use pkgdown
#usethis::use_pkgdown()

# Run to build the website
# pkgdown::build_site()


# -------------------------------------------------------------------------

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

x <- snglmm(y ~ Sepal.Width,
            iris,
            spatial = "exp",
            D = D)
x
summary(x)


pos <- glmmTMB::numFactor(coord)
iris$g <- 1
iris$pos <- pos
fit <- glmmTMB::glmmTMB(y ~ Sepal.Width + exp(0 + pos|g), iris, dispformula = ~0)

