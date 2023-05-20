# Tests for milags_functions
#
# Author: MNECKI
###############################################################################


# tests for get_n0

test_that("Getting initial biomass works", {

  # data
  biom <- c(27727155, 28619420, 31695970, 35766010, 40113020, 44429830, 47754910, 49705080)
  method_1 <- "first.observation"
  method_2 <- "minimal.observation"
  first_calc <- biom[1]
  minim_calc <- min(biom)

  expect_equal(get_n0(biom, method_1), first_calc)
  expect_equal(get_n0(biom, method_2), minim_calc)
})


# tests for fit_exp_lag_to_curve

test_that("", {}

)


# tests for compare_algorithms

test_that("Basic algorithm comparison works", {


  expect_equal(get_n0(biom, method_1), first_calc)
  expect_equal(get_n0(biom, method_1), first_calc)
  expect_equal(get_n0(biom, method_1), first_calc)
})


# tests for choose_lag_fit_algorithm_baranyi

test_that("Calculating nls of given model works", {
  ref_gr <- data.frame(LOG10N = 40:60,
                    t = 0:20)
  ref_log10n <- sample(1:5, 1)
  ref_lag <- runif(1, min = 0, max = 1)
  ref_mumax <- runif(1, min = 0, max = 1)
  ref_log10max <- runif(1, min = 0, max = 1)
  ref_maxiter <- sample(10:20, 1)


  expect_equal()
})
