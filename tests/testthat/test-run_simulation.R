MESA <- dplyr::tibble(
  id = sort(rep(1:100,3)),
  visit_number = rep(1:3,100),
  Distance = rexp(300),
  Time = rexp(300)
)

pars <- list(alpha = 23,
             shape_one = 1,
             scale_one = .5,
             shape_two = 1,
             scale_two = .3,
             delta = -2.2,
             beta = 1.2,
             sigma = 1)

test_that("run_simulation works as intended", {
  
 out <- dplyr::filter(run_simulation(data = MESA),
                      stringr::str_detect(Parameter, "shape", negate = TRUE))
 expect_true(all(out$Coverage))
 
})
