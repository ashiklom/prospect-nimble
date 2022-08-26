library(nimble)

stopifnot(
  requireNamespace("rrtm", quietly = TRUE)
)

e1_approx <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(1))

    A <- log((0.56146 / x + 0.65) * (1 + x))
    B <- x^4 * exp(7.7 * x) * (2 + x)^3.7

    return((A^-7.7 + B)^-0.13)
  }
)

prospect5 <- nimbleFunction(
  run = function(N = double(0), Cab = double(0), Car = double(0), Cw = double(0), Cm = double(0),
                 dataspec_p5 = double(2), talf = double(1), t12 = double(1), t21 = double(1)) {

    cc <- matrix(NA, nrow = 5, ncol = 1)
    k <- numeric(length = 2101)

    Cbrown <- 0
    cc[,] <- t(c(Cab, Car, Cbrown, Cw, Cm) / N)

    k[] <- dataspec_p5[,] %*% cc[,]

    trans <- (1 - k)*exp(-k) + k^2 *e1_approx(k)

    ralf <- 1 - talf
    r12 <- 1 - t12
    r21 <- 1 - t21

    denom <- 1 - (r21 ^ 2) * (trans ^ 2)
    Ta <- talf * trans * t21 / denom
    Ra <- ralf + r21 * trans * Ta

    tt <- t12 * trans * t21 / denom
    rr <- r12 + r21 * trans * tt

    r <- rr
    t <- tt
    D <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    r2 <- r ^ 2
    t2 <- t ^ 2
    va <- (1 + r2 - t2 + D) / (2 * r)
    vb <- (1 - r2 + t2 + D) / (2 * t)

    vbNN <- vb ^ (N - 1)
    vbNN2 <- vbNN ^ 2
    va2 <- va ^ 2
    denomx <- va2 * vbNN2 - 1
    Rsub <- va * (vbNN2 - 1) / denomx
    Tsub <- vbNN * (va2 - 1) / denomx

    denomy <- 1 - Rsub * rr

    # transmittance <- Ta * Tsub / denomy
    reflectance <- Ra + Ta * Rsub * tt / denomy

    returnType(double(1))
    return(reflectance)
  }
)

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

plot(model$reflectance, type = 'l')
lines(rrtm::prospect5(
  N = 2, Cab = 40, Car = 10,
  Cbrown = 0, Cw = 0.002, Cm = 0.001
)[["reflectance"]], col = 'red', lty = "dashed")

compiled_Pmodel <- compileNimble(model, showCompilerOutput = TRUE)
compiled_Pmodel$setData(list(
  talf = rrtm:::p45_talf,
  t12 = rrtm:::p45_t12,
  t21 = rrtm:::p45_t21,
  dataspec_p5 = rrtm:::dataspec_p5
))
compiled_Pmodel$simulate()

plot(compiled_Pmodel$reflectance, type = 'l')

##################################################
# Begin my code...
##################################################

obs <- read.csv(
  "https://raw.githubusercontent.com/ashiklom/spectra_db/master/nasa_fft/spectra/nasa_fft%7CAK06_ABBA_M2%7C2009.csvy",
  comment.char = "#",
  col.names = c("wavelength", "observed")
)
