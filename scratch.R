# Create hull piecewise, so it captues more absorption features
plot(obsr[1530:2101], type = "l")
ch1 <- sort(chull(c(obsr[1:269])))
ch2 <- sort(chull(obsr[270:1040])) + 269
ch3 <- sort(chull(obsr[1041:1529])) + 1040
ch4 <- sort(chull(c(obsr[1530:2101], 0))) + 1529
chx <- obsr

chid <- c(ch1, ch2, ch3, ch4)
chid <- chid[chid < 1000]
plot(chx, type = "l")
points(chid, chx[chid], col = "red", pch = 19)
af <- approx(chid, chx[chid], 1:2101)
lines(af$y, col = "blue")

# Custom hull
myx <- c(1, 158, 369, 517, 677, 728, 857, 879, 898, 1235, 1259,
         1289, 1431, 1439, 1800, 1827, 1846, 2097, 2101)
myy <- obsr[myx]
myaf <- approx(myx, myy, 1:2101)
myhull <- myaf$y
plot(obsr, type = "l")
points(myx, myy, pch = 19, col = "red")
lines(y ~ x, data = myaf, col = "blue")

plot(400:2500, obsr / myhull, type = "l")

pnull <- rrtm::prospect5(1.4, 50, 10, 0, 0.005, 0.01, e1fun = gsl::expint_E1)[["reflectance"]]
plot(400:2500, pnull, type = "l", ylim = c(0, 0.45))

plot(400:2500, obsr / pnull, type = "l")

d5 <- rrtm:::dataspec_p5
d5n <- apply(d5, 2, function(x) x / max(x))
matplot(400:2500, d5n, type = "l", lty = 1, xlim = c(400, 600))
legend("top", colnames(d5n), col = 1:5, lty = 1)
