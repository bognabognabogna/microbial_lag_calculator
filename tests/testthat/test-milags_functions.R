# Tests for milags_functions
#
# Author: jungwirt
###############################################################################


# tests for get_n0
context("Test the get_n0 function")


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

context("Test the calc_lag function")

test_that("Calculating if lag of given model works", {


  expect_equal()
})



context("Test the get_def_pars function")

test_that("Calculating if getting default parameters works", {

  expect_equal(get_def_pars(), list("logistic",
                                    "first.observation",
                                    "local.regression",
                                    10^2,
                                    3,
                                    NULL,
                                    NULL,
                                    "auto",
                                    100))
})




context("Test the lag_biomass_incr function")

test_that("Calculating if fitting the lag to multiple growth curves based on the biomass increase method works", {

  database <- "testing_data.csv"
  test_df <- read.csv2(database) %>% filter(grepl('biomass_incr', curve_id))
  test_threshold <- 5000
  test_n0 <- 0
  data_test <- test_df %>% filter(FALSE) %>% mutate(lag = numeric(0))
  for (this_curve_id in unique(test_df$curve_id)) {
    data_this_curve <- test_df %>%
      filter(curve_id == this_curve_id) %>%
      left_join(test_n0, by = "curve_id") %>%
      mutate(incr_from_n0 = biomass - test_n0)
    
    
    #find where second derivative is maximal
    threshold_diff <- which(data_this_curve$incr_from_n0 >= test_threshold)
    first_threshold_diff <- threshold_diff[1]
    lag_this_curve <- data_this_curve$time[first_threshold_diff]
    data_this_curve <- data_this_curve %>%
      mutate(
        lag = round(lag_this_curve,1))
    data_test <- rbind(data_test, data_this_curve)
}
  expect_equal(lag_biomass_incr(test_df, test_threshold, test_n0), data_test)
})


context("Test the plot_data function")

test_that("Plotting growth curve works", {
  
  # data
  database <- "testing_data.csv"
  test_df <- read.csv2(database) %>% filter(grepl('exponential', curve_id)) %>% select(time, biomass)
  data_new <- test_df %>%
    mutate(log10_biomass = log10(biomass))
  g_test <- ggplot(data_new)  +
    geom_line(aes(x = time, y = log10_biomass), col = "blue") +
    geom_point(aes(x = time, y = log10_biomass), col = "blue") +
    xlab("time [h]") +
    ylab("Log10(biomass)") +
    theme(axis.text.y.right = element_text(colour="black"),
          axis.text.y = element_text(colour="blue"),
          axis.title.y = element_text(colour="blue"),
          axis.title.y.right = element_text(colour="black"))
  
  expect_equal(plot_data(test_df), g_test)
 
})


