library(microbenchmark)

microbenchmark::microbenchmark(
  test_threshold_2(HLCYG_physeq_data, repeats = 10, t_max = 200, t_step = 10, group = "location"),
  times=100
)

microbenchmark(
  test_threshold_2(HLCYG_physeq_data, repeats = c(5,10), t_min = 20, t_max = 250, t_step = 1, group = "location"),
  times=1
)