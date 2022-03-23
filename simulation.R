################################################################################
# Future_apply
sim <- function(x, scen) {
  set.seed(scen$seed[x])
  alpha <- scen$alpha[x]
  deg <- scen$deg[x]
  nsim <- scen$nsim[x]
  n <- scen$n[x]
  d <- scen$d[x]
  method <- scen$method[x]
  la <- c(rep("Continuous", d / 3), rep("Binary", d / 3), rep("Ordinal", d / 3))
  data <- generate_DAG(n, d, deg, label = la, lB = -0.5, uB = 0.5)
  X <- data$X
  A <- as(dag2cpdag(data$A), "matrix")
  t0 <- as.numeric(t(A))
  
  start_time <- proc.time()
  if (method == 0) {
    # Oracle PC
    sig <- data$sigmahat
    sig_diff <- 0
    pc.fit <- pc(suffStat = list(C = sig, n = n),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha = alpha, labels = colnames(X), skel.method = "stable.fast",
                 verbose = F)
    tsig <- as.numeric(as(pc.fit, "amat"))
  } else if (method == 1) {
    # Vanilla PC
    sig <- cor(X)
    sig_diff <- sum(abs(sig - data$sigmahat))
    pc.fit <- pc(suffStat = list(C = sig, n = n),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha = alpha, labels = colnames(X), skel.method = "stable.fast",
                 verbose = F)
    tsig <- as.numeric(as(pc.fit, "amat"))
  } else if (method == 2) {
    # Rank PC
    sig <- sin(pi / 2 * pcaPP::cor.fk(X))
    sig_diff <- sum(abs(sig - data$sigmahat))
    pc.fit <- pc(suffStat = list(C = sig, n = n),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha = alpha, labels = colnames(X), skel.method = "stable.fast",
                 verbose = F)
    tsig <- as.numeric(as(pc.fit, "amat"))
  } else if (method == 3) {
    # Copula PC
    cop.obj <- inferCopulaModel(X, nsamp = 1000, S0 = diag(d) / n, verb = F)
    C_samples <- cop.obj$C.psamp[, , 501:1000]
    corr.cop <- apply(C_samples, c(1, 2), mean)
    sig_diff <- sum(abs(corr.cop - data$sigmahat))
    pc.fit <- pc(suffStat = list(C = corr.cop, n = n),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha = alpha, labels = colnames(X), skel.method = "stable.fast",
                 verbose = F)
    tsig <- as.numeric(as(pc.fit, "amat"))
  } else if (method == 4) {
    # LR PC
    skel <- suppressWarnings(pc.skel(X, method = "comb.mm", alpha = alpha))
    pc.fit <- pc.or(skel)
    tsig <- pc.fit$G
    sig_diff <- NA
    tsig[tsig == 2] <- 0
    tsig[tsig == 3] <- 1
    tsig <- as.numeric(tsig)
  } else if (method == 5) {
    labels <- label_fun(X)
    # deltas <- delta_fun(X, labels)
    # sig <- latent_pc(X, labels, deltas)
    sig <- latent_pc(X, labels)
    sig_diff <- sum(abs(sig - data$sigmahat))
    sig <- as.matrix(nearPD(sig, corr = TRUE, maxit = 10000)$mat)
    pc.fit <- pc(suffStat = list(C = sig, n = n),
                 indepTest = gaussCItest, ## indep.test: partial correlations
                 alpha = alpha, labels = colnames(X), skel.method = "stable.fast",
                 verbose = F)
    tsig <- as.numeric(as(pc.fit, "amat"))
  }
  
  time <- c(proc.time() - start_time)[3]
  tp <- sum(t0[which(tsig == 1)] == 1)
  tn <- sum(t0[which(tsig == 0)] == 0)
  fp <- sum(tsig[which(t0 == 0)] == 1)
  fn <- sum(tsig[which(t0 == 1)] == 0)
  TPR <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  FPR <- ifelse(fp + tn == 0, 0, fp / (fp + tn))
  SHM <- sum(abs(tsig - t0))
  out <- matrix(c(scen$seed[x], alpha, deg, n, d, method, time, sig_diff, TPR, FPR, SHM), nrow = 1)
  colnames(out) <- c("seed", "alpha", "deg", "n", "d", "method", "time", "sig_diff", "TPR", "FPR", "SHM")
  return(out)
}

# # Parallel
# # library(pbivnorm)
# # library(rootSolve)
# library(pcalg)
# library(Matrix)
# # library(Rcpp)
# # library(RcppArmadillo)
# # library(RcppNumerical)
# library(MXM)
# library(sbgcop)
# # library(mixcausal)
# library(mvtnorm)
# library(pcaPP)


# alpha <- sort(c(0.01, 0.9^(seq(50, 1, by = -1))))
alpha <- 0.05
deg <- c(3, 5)
seed <- 1234:(1234 + 49)
method <- c(0:5)
scen <- expand.grid(seed, alpha, deg, method)
# n <- c(50, 100, 150, 200, 1000, 200, 200)
# d <- c(9,  27,  81,  243, 243, 729, 2187)
n <- 200
d <- 243
scen <- data.frame(rep(scen[, 1], length(n)),
                   rep(scen[, 2], length(n)),
                   rep(scen[, 3], length(n)),
                   rep(scen[, 4], length(n)),
                   rep(n, each = nrow(scen)), rep(d, each = nrow(scen)))
colnames(scen) <- c("seed", "alpha", "deg", "method", "n", "d")
scen <- subset(scen, !(method == 3 & n < d))
rownames(scen) <- paste0(1:nrow(scen))

library(future.apply)
source("utility.R")
source("copula_pc.R")
source("latent_pc.R")
source("data_generation.R")
plan(multisession)
start <- proc.time()
result <- future_lapply(1:nrow(scen), FUN = sim, future.seed = NULL,
                        future.packages = c("pcalg", "Matrix", "MXM", "sbgcop",
                                            "mvtnorm", "pcaPP"), scen = scen)
c(proc.time() - start)[3]
aaa <- as.data.frame(do.call(rbind, result))
write.csv(aaa, file="result_M3_200_243_0.05_0.5_1234.csv")

out <- data.frame()
for (j in 1:2) {
  for (i in 1:6) {
    b <- subset(aaa, method == (0:5)[i] & deg == c(3, 5)[j])
    out <- rbind(out, colMeans(b, na.rm = T)[3:11])
  }
}
colnames(out) <- colnames(aaa)[3:11]
# out$time <- out$time / 60
round(out, 3)
write.csv(out, file="summary_M3_200_243_0.05_0.5_1234.csv")
