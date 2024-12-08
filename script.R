
## 3-d density and densityratio plots

library(mvtnorm)
library(rgl)

mu1 <- c(0, 0)
mu2 <- c(0.25, -0.25)
sigma1 <- diag(2)
sigma2 <- matrix(c(2,0.5,0.5, 2),2)

f1 <- \(x, y) dmvnorm(cbind(x, y), mu1, sigma1)
f2 <- \(x, y) dmvnorm(cbind(x, y), mu2, sigma2)

x <- seq(-5, 5, length.out = 100)
y <- seq(-5, 5, length.out = 100)
z1 <- outer(x, y, f1)
z2 <- outer(x, y, f2)

viewmat <- matrix(c(0.78263121843338,0.173177808523178,-0.597911298274994,0,-0.621819615364075,0.261917144060135,-0.738064885139465,0,0.0287867486476898,0.949425637722015,0.312669903039932,0,0,0,0,1),4)

open3d(userMatrix = viewmat,windowRect = c(50,50,1600,1200))
par3d(cex = 1.5)
persp3d(x, y, z1, col = "#93effa",
        front = "lines", back = "lines",
        lit = FALSE, zlim = c(0, 0.2),
        xlab = expression(x[1]),
        ylab = expression(x[2]),
        zlab = expression(P(x)))

persp3d(x, y, z2, col = "#034e57",
        front = "lines", back = "lines",
        lit = FALSE, zlim = c(0, 0.2),
        xlab = expression(x[1]),
        ylab = expression(x[2]),
        zlab = expression(P(x)),
        cex = 3)
legend3d(0.57, 0.65, legend = c(expression(italic(P[scriptstyle("nu")](x))), expression(italic(P[de](x)))),
         col = c("#93effa", "#034e57"), lty = 1, lwd = 2, bty = "n", cex = 2)

rgl.snapshot("densities.png")
close3d()

open3d(userMatrix = viewmat,windowRect = c(50,50,1600,1200))
par3d(cex = 1.5)
persp3d(x, y, z1/z2, col = "lightgrey",
        front = "lines", back = "lines",
        lit = FALSE,
        xlab = expression(x[1]),
        ylab = expression(x[2]),
        zlab = expression(r(x)),
        cex = 3)

set.seed(123)
fit <- densityratio::ulsif(
  rmvnorm(1000, mu1, sigma1),
  rmvnorm(1000, mu2, sigma2),
  parallel = TRUE, nthreads = 18
)
rhat <- outer(x, y, \(x, y) predict(fit, cbind(x, y) |> unname()))
persp3d(x, y, rhat, col = "#de0277",
        front = "lines", back = "lines",
        lit = FALSE, add = TRUE,
        xlab = expression(x[1]),
        ylab = expression(x[2]),
        zlab = expression(r(x)),
        cex = 3)
legend3d(0.57, 0.65, legend = c(expression(italic(r(x))), expression(italic(hat(r)(x)))),
         col = c("lightgrey", "#de0277"), lty = 1, lwd = 2, bty = "n", cex = 2)
rgl.snapshot("densityratioplot.png")
close3d()

## Computation times
set.seed(123)

n.iter <- 5
n <- c(500, 1000, 2000, 4000, 8000, 16000)
dim <- c(3, 6, 12, 24, 48, 96)

mu <- c(1,2,3)
sigma <- diag(3)

n_times <- sapply(n, function(n) {
  d1 <- rmvnorm(n, mu, sigma)
  d2 <- rmvnorm(n, mu, sigma)
  bench::mark({
    f = densityratio::ulsif(
      d1, d2,
      parallel = TRUE, nthreads = 18,
      progressbar = FALSE
    )
  }, min_iterations = n.iter)$time[[1]] |> mean()
})

dim_times <- sapply(dim, function(d) {
  set.seed(123)
  d1 <- rmvnorm(1000, 1:d, diag(d))
  d2 <- rmvnorm(1000, 1:d, diag(d))
  bench::mark({
    f = densityratio::ulsif(
      d1, d2,
      parallel = TRUE, nthreads = 18,
      progressbar = FALSE
    )
  }, min_iterations = n.iter)$time[[1]] |> mean()
})

saveRDS(
  list(n = n,
       dim = dim,
       n_times = n_times,
       dim_times = dim_times),
  "computation_times.rds"
)
