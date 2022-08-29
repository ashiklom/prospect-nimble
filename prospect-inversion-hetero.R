# Try single PROSPECT inversion
source("nimble-functions.R")

obs <- read.csv(
  "https://raw.githubusercontent.com/ashiklom/spectra_db/master/nasa_fft/spectra/nasa_fft%7CAK06_ABBA_M2%7C2009.csvy",
  comment.char = "#",
  col.names = c("wavelength", "observed")
)

fit_prospect5_hetero <- nimbleCode({
  N ~ T(dnorm(1.0, sd=1.0), 1.0, Inf)
  Cab ~ T(dnorm(40, sd=10), 0, Inf)
  Car ~ T(dnorm(10, sd=10), 0, Inf)
  Cw ~ T(dnorm(0.01, sd=0.01), 0, Inf)
  Cm ~ T(dnorm(0.01, sd=0.01), 0, Inf)

  rslope ~ dnorm(0, sd=1)
  rint ~ dexp(10)

  mod[1:2101] <- prospect5(N, Cab, Car, Cw, Cm,
                           dataspec_p5[,], talf[], t12[], t21[])
  for (i in 1:2101) {
    # Heteroskedastic variance model
    obs[i] ~ dnorm(mod[i], sd=rint + rslope * mod[i])
  }
})

Nwl <- 2101
fpm <- nimbleModel(fit_prospect5_hetero, dimensions = list(
  reflectance = Nwl,
  obs = Nwl,
  dataspec_p5 = c(Nwl, 5),
  talf = Nwl,
  t12 = Nwl,
  t21 = Nwl
))
fpm$setData(list(
  obs = subset(obs, wavelength >= 400)[["observed"]],
  talf = rrtm:::p45_talf,
  t12 = rrtm:::p45_t12,
  t21 = rrtm:::p45_t21,
  dataspec_p5 = rrtm:::dataspec_p5
))
cfp <- compileNimble(fpm)
rmcmc <- configureMCMC(fpm)
bmcmc <- buildMCMC(rmcmc)
cmcmc <- compileNimble(bmcmc, project=cfp)

samps <- runMCMC(cmcmc, nburnin = 2000, niter = 15000, nchains = 3)

# Evaluate predictions
samps_mcmc <- coda::as.mcmc.list(lapply(samps, coda::as.mcmc))
plot(samps_mcmc)

sampsum <- summary(samps_mcmc)

bestfit <- do.call(rrtm::prospect5, c(as.list(sampsum$statistics[1:5,"Mean"]), Cbrown = 0))$reflectance
plot(observed ~ wavelength, data = obs, type = 'l')
lines(400:2500, bestfit, col = "red")
