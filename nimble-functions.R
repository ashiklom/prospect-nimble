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

    Cbrown <- 0

    k <- dataspec_p5[,1] * Cab/N +
      dataspec_p5[,2] * Car/N +
      dataspec_p5[,3] * Cbrown/N +
      dataspec_p5[,4] * Cw/N +
      dataspec_p5[,5] * Cm/N

    trans <- e1_approx(k) * k^2 + (1 - k) * exp(-k)

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
