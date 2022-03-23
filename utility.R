
### Expanded data set, index and delta
## Input
# X: data set with each column representing a variable
# X should not contain NA or non-number
# Binary variable should only contain 0 or 1
# Ordinal variable should only contain 0, 1, 2, ...
# label: Labels for variables
## Output
# A list of the following
# Y: Data set with ordinal variables expanded to dummies based on the binary cuts (>=1, >=2, ...)
# index: A vector connecting variables in X and Y
# delta: A list of means for each variable. For continuous variables, the list element contains the mean
# For binary variables, the list element contains the normal quantile of 1 - mean
# For ordinal variables, the list elements contain the normal quantile of 1 - mean
# of binary cuts of the ordinal variables for >=1, >=2, ...
transform_fun <- function(X, label) {
  Y <- data.frame(rep(NA, nrow(X)))
  index <- c()
  delta <- list()
  for (i in 1:ncol(X)) {
    if (label[i] == "Continuous") {
      index <- c(index, i)
      Y <- data.frame(Y, X[, i])
      delta <- append(delta, mean(X[, i]))
    } else if (label[i] == "Binary") {
      index <- c(index, i)
      Y <- data.frame(Y, X[, i])
      delta <- append(delta, qnorm(1 - mean(X[, i])))
    } else {
      level <- sort(unique(X[, i]))
      temp <- c()
      for (j in 1:max(level)) {
        a <- X[, i] >= j
        Y <- data.frame(Y, as.numeric(a))
        temp <- c(temp, qnorm(1 - mean(a)))
      }
      index <- c(index, rep(i, max(level)))
      delta <- append(delta, list(temp))
    }
  }
  return(list(Y[, -1], index, delta))
}

### Kendall's tau
## Input
# X: data set with more than two columns
## Output
# Kendall's tau correlation coefficient between two columns
# Function modified from the Kendalltau function in the latentcor R package
tau_fun <- function(X) {
  # Function from the latentcor package
  n_x <- function(x, n) {
    if (length(unique(x) != n)) {
      x.info <- rle(sort(x))
      t_x <- x.info$lengths[x.info$lengths > 1]
      n_x <- sum(t_x * (t_x - 1) / 2)
    } else {
      n_x <- 0
    }
    return(n_x)
  }
  
  n <- nrow(X)
  n0 <- n * (n - 1) / 2
  n_X <- apply(X, 2, function(x) {n_x(x = x, n)})
  n_X_sqd <- sqrt(n0 - n_X)
  K_b <- pcaPP::cor.fk(X)
  K_b.lower <- K_b[lower.tri(K_b)]
  btoa <- n_X_sqd[row(K_b)[lower.tri(K_b)]] * n_X_sqd[col(K_b)[lower.tri(K_b)]] / n0
  K_a.lower <- K_b.lower * btoa
  out <- matrix(0, ncol(X), ncol(X))
  out[lower.tri(out, diag = F)] <- K_a.lower
  out <- out + t(out)
  diag(out) <- 1
  return(out)
}

######## Converting Kendall's tau correlation to Pearson's correlation ########
### CC: Continuous and continuous
## Input
# tau: Kendall's tau correlation coefficient
## Output
# Pearson's correlation coefficient
CC_fun <- function(tau) {
  sin(pi / 2 * tau)
}

### BB: Binary and binary
## Input
# x: Pearson's correlation coefficient to be solve for
# tau_x: Kendall's tau correlation coefficient
# delta_x: List of means of the two binary variables
## Output
# To be solved for
BB_fun <- function(x, tau_x, delta_x) {
  delta_x <- unlist(delta_x)
  (2 * pmvnorm(upper = delta_x,
               mean = 0,
               sigma = matrix(c(1, x, x, 1), nrow = 2),
               algorithm = Miwa(steps = 128)) - 
      2 * pnorm(delta_x[1]) * pnorm(delta_x[2]) - tau_x)^2
}

### BC: Binary and continuous, or vice versa
## Input
# x: Pearson's correlation coefficient to be solve for
# tau_x: Kendall's tau correlation coefficient
# delta_x: List of mean of the binary variable
## Output
# To be solved for
BC_fun <- function(x, tau_x, delta_x) {
  delta_x <- unlist(delta_x)
  (4 * pmvnorm(upper = c(delta_x, 0),
               mean = 0,
               sigma = matrix(c(1, x / sqrt(2), x / sqrt(2), 1), nrow = 2),
               algorithm = Miwa(steps = 128)) -
      2 * pnorm(delta_x) - tau_x)^2
}

# CO: Continuous and ordinal, or vice versa
## Input
# tau_o: Kendall's tau correlation between continuous and binary cuts from the ordinal variable
# label_o: Labels for the two columns
# delta_o: List of means of the binary cuts (>=1, >=2, ...) from the ordinal variable
## Output
# Pearson's correlation coefficient
CO_fun <- function(tau_o, label_o, delta_o) {
  delta_o <- unlist(delta_o)
  z <- c()
  id <- which(label_o == "Ordinal")
  level <- length(delta_o)
  for (h in 1:level) {
    z <- c(z, optimize(BC_fun, c(-0.9999, 0.9999),
                       tau_x = ifelse(id == 1, tau_o[max(level) + 1, h], tau_o[1, h + 1]),
                       delta_x = delta_o[h])$minimum)
  }
  mean(z)
}

# BO: Binary and ordinal, or vice versa
## Input
# tau_o: Kendall's tau correlation between continuous and binary cuts from the ordinal variable
# label_o: Labels for the two columns
# delta_o: List of mean of the binary variable and
# the normal quantile of 1 - means of the binary cuts (>=1, >=2, ...) from the ordinal variable
## Output
# Pearson's correlation coefficient
BO_fun <- function(tau_o, label_o, delta_o) {
  z <- c()
  id <- which(label_o == "Ordinal")
  level <- length(delta_o[[id]])
  for (h in 1:level) {
    z <- c(z, optimize(BB_fun, c(-0.9999, 0.9999),
                       tau_x = ifelse(id == 1, tau_o[max(level) + 1, h], tau_o[1, h + 1]),
                       delta_x = c(delta_o[[-id]], unlist(delta_o[[id]])[h]))$minimum)
  }
  mean(z)
}

# OO: Ordinal and ordinal, or vice versa
## Input
# tau_o: Kendall's tau correlation between binary cuts from the two ordinal variable
# delta_o: List of the normal quantile of 1 - means of
# the binary cuts (>=1, >=2, ...) from two ordinal variables
## Output
# Pearson's correlation coefficient
OO_fun <- function(tau_o, delta_o) {
  z <- c()
  level1 <- length(delta_o[[1]])
  level2 <- length(delta_o[[2]])
  for (h in 1:max(level1)) { 
    for (k in 1:max(level2)) {
      z <- c(z, optimize(BB_fun, c(-0.9999, 0.9999),
                         tau_x = tau_o[h, level1 + k],
                         delta_x = c(delta_o[[1]][h], delta_o[[2]][k]))$minimum)
    }
  }
  mean(z)
}
