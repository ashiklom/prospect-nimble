# Test a PROSPECT simulation
source("nimble-functions.R")

run_prospect5 <- nimbleCode({
  Cab <- 40
  Car <- 10
  Cw <- 0.002
  Cm <- 0.001
  N <- 2
  reflectance[1:2101] <- prospect5(N,Cab,Car,Cw,Cm,
                             dataspec_p5[,], talf[],t12[],t21[])
})

Nwl <- 2101
model <- nimbleModel(run_prospect5, dimensions = list(
  reflectance = c(Nwl),
  dataspec_p5 = c(Nwl,5),
  talf = c(Nwl),
  t12 = c(Nwl),
  t21 = c(Nwl)
))
model$setData(list(
  talf = rrtm:::p45_talf,
  t12 = rrtm:::p45_t12,
  t21 = rrtm:::p45_t21,
  dataspec_p5 = rrtm:::dataspec_p5
))
model$simulate()

# plot(model$reflectance, type = 'l')
# lines(rrtm::prospect5(
#   N = 2, Cab = 40, Car = 10,
#   Cbrown = 0, Cw = 0.002, Cm = 0.001
# )[["reflectance"]], col = 'red', lty = "dashed")

compiled_Pmodel <- compileNimble(
  model, prospect5, e1_approx,
  showCompilerOutput = TRUE,
  dirName = "nimble_cppcode",
  projectName = "MYPROJ"
)
compiled_Pmodel$model$simulate()
plot(compiled_Pmodel$model$reflectance, type = 'l')
