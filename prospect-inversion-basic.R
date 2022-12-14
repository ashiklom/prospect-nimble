# Try single PROSPECT inversion
source("nimble-functions.R")

sample_id <- "nasa_fft|AK01_ACRU_B|2009.csvy"

obs <- read.csv(
  paste0("https://raw.githubusercontent.com/ashiklom/spectra_db/master/nasa_fft/spectra/", sample_id),
  comment.char = "#",
  col.names = c("wavelength", "observed")
)
obsr <- subset(obs, wavelength >= 400)[["observed"]]

fit_prospect5 <- nimbleCode({
  N ~ T(dnorm(1.0, sd=1.0), 1.0)
  Cab ~ T(dnorm(40, sd=10), 0)
  Car ~ T(dnorm(10, sd=10), 0)
  Cw ~ T(dnorm(0.01, sd=0.01), 0)
  Cm ~ T(dnorm(0.01, sd=0.01), 0)
  tau ~ dgamma(0.01, 0.01)

  mod[1:2101] <- prospect5(N, Cab, Car, Cw, Cm,
                           dataspec_p5[,], talf[], t12[], t21[])
  for (i in 1:2101) {
    obs[i] ~ dnorm(mod[i], tau=tau)
  }
})

Nwl <- 2101
dims <- list(
  reflectance = Nwl,
  obs = Nwl,
  dataspec_p5 = c(Nwl,5),
  talf = Nwl,
  t12 = Nwl,
  t21 = Nwl
)
ndata <- list(
  obs = obsr,
  talf = rrtm:::p45_talf,
  t12 = rrtm:::p45_t12,
  t21 = rrtm:::p45_t21,
  dataspec_p5 = rrtm:::dataspec_p5
)

samps <- nimbleMCMC(
  fit_prospect5,
  data = ndata,
  dimensions = dims,
  samplesAsCodaMCMC = TRUE
)

# Evaluate predictions
iburn <- 2000:nrow(samps)
par(mfrow = c(2, 2))
for (j in c("N", "Cab", "Cw", "Cm")) {
  plot(samps[iburn, j], type = 'l', ylab = j)
}

par(mfrow = c(2, 2))
for (j in c("N", "Cab", "Cw", "Cm")) {
  hist(samps[iburn, j], main = j, xlab = j)
}

bestfit <- colMeans(samps[iburn, c("N", "Cab", "Car", "Cw", "Cm")])
plot(observed ~ wavelength, data = subset(obs, wavelength > 400 & wavelength < 700), type = 'l')
lines(400:2500, do.call(rrtm::prospect5, c(as.list(bestfit), Cbrown = 0.00))$reflectance,
      col = "red")
