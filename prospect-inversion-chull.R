# Try single PROSPECT inversion based on convex hull of observation
# TODO: Try custom convex hull points based on known absorption features
source("nimble-functions.R")

obs <- read.csv(
  "https://raw.githubusercontent.com/ashiklom/spectra_db/master/nasa_fft/spectra/nasa_fft%7CAK06_ABBA_M2%7C2009.csvy",
  comment.char = "#",
  col.names = c("wavelength", "observed")
)

fit_prospect5_chull <- nimbleCode({
  N ~ T(dnorm(1.0, sd=1.0), 1.0, Inf)
  Cab ~ T(dnorm(40, sd=10), 0, Inf)
  Car ~ T(dnorm(10, sd=10), 0, Inf)
  Cw ~ T(dnorm(0.01, sd=0.01), 0, Inf)
  Cm ~ T(dnorm(0.01, sd=0.01), 0, Inf)

  tau ~ dgamma(0.01, 0.01)

  mod[1:2101] <- prospect5(N, Cab, Car, Cw, Cm,
                           dataspec_p5[,], talf[], t12[], t21[]) / hull[1:2101]
  for (i in 1:2101) {
    obs[i] ~ dnorm(mod[i], tau = tau)
  }
})

obsr <- subset(obs, wavelength >= 400)[["observed"]]
plot(400:2500, obsr, type = 'l')
Nwl <- 2101
stopifnot(length(obsr) == Nwl)

# Compute convex hull
# TODO: Compare results with/without 0-padding
chid <- sort(chull(c(0, obsr, 0)))
plot(c(0, obsr, 0), type = "l")
points(chid, c(0, obsr, 0)[chid], col = "red", pch = 19)

af <- approx(chid, c(0, obsr, 0)[chid], 2:(Nwl+1))
hull <- af$y
plot(400:2500, obsr / hull, type = "l")

fpm <- nimbleModel(fit_prospect5_chull, dimensions = list(
  hull = Nwl,
  obs = Nwl,
  dataspec_p5 = c(Nwl, 5),
  talf = Nwl,
  t12 = Nwl,
  t21 = Nwl
))
fpm$setData(list(
  obs = obsr / hull,
  hull = hull,
  talf = rrtm:::p45_talf,
  t12 = rrtm:::p45_t12,
  t21 = rrtm:::p45_t21,
  dataspec_p5 = rrtm:::dataspec_p5
))
cfp <- compileNimble(fpm)
rmcmc <- configureMCMC(fpm)
bmcmc <- buildMCMC(rmcmc)
cmcmc <- compileNimble(bmcmc, project=cfp)

# samps <- runMCMC(cmcmc, nburnin = 2000, niter = 15000, nchains = 3)
samps <- runMCMC(cmcmc, samplesAsCodaMCMC = TRUE)
samps_mcmc <- window(samps, start = 2000)
plot(samps_mcmc)

sampsum <- summary(samps_mcmc)
means <- sampsum$statistics[1:5, "Mean"]

bestfit <- do.call(rrtm::prospect5, c(as.list(means), Cbrown = 0))$reflectance
plot(observed ~ wavelength, data = obs, type = 'l')
lines(400:2500, bestfit, col = "red")
