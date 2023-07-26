
# build package -----------------------------------------------------------

usethis::use_mit_license(copyright_holder = "Akira Terui")
usethis::use_roxygen_md()
usethis::use_namespace()
usethis::use_package_doc()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)


# check syntax ------------------------------------------------------------

# lintr::lint_package()

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

library(dplyr)
set.seed(12)
N <- 200
theta <- 10
eps <- rnorm(N, sd = 0.1)
X <- cbind(1, rnorm(N), rnorm(N))
beta <- c(0.01, 0.8, 0.2)

D <- cbind(runif(N, 0, 1000), runif(N, 0, 1000)) %>%
  dist(diag = TRUE, upper = TRUE) %>%
  data.matrix()
S <- 0.1 * exp(-D / theta)

y_hat <- X %*% beta
y <- mvtnorm::rmvnorm(1, mean = y_hat, sigma = S) %>% c()

df0 <- dplyr::tibble(x1 = X[,2], x2 = X[,3], y = y)

semtmb(y ~ x1 + x2, data = df0)
