### Label
## Input
# X: data set with each column representing a variable
# X should not contain NA or non-number
# Binary variable should only contain 0 or 1
# Ordinal variable should only contain 0, 1, 2, ...
## Output
# A vector of labels for all variables: "Binary", "Continuous", "Ordinal"
label_fun <- function(X) {
  if (any(is.na(X))) {
    stop("Data should not contain NA")
  } else if (!all(is.numeric(as.matrix(X)))) {
    stop("Data should only contain numbers")
  } else {
    y <- c()
    for (i in 1:ncol(X)) {
      temp <- X[, i]
      if (length(table(temp)) == 0) {
        stop("Data contain a variable with no value")
      } else if (length(table(temp)) == 1) {
        stop("Data contain a variable with a constant value")
      } else if (length(table(temp)) == 2) {
        if (all(temp %in% c(0, 1))) {
          y <- c(y, "Binary")
        } else {
          y <- c(y, "Continuous")
        }
      } else {
        level <- sort(unique(temp))
        if (all(temp %in% 0:max(level))) {
          y <- c(y, "Ordinal")
        } else {
          y <- c(y, "Continuous")
        }
      }
    }
  }
  return(y)
}

### Latent PC algorithm
## Input
# X: Data set with each column representing a variable
# label: Label for each column
## Output
# Pearson correlation matrix
latent_pc <- function(X, label) {
  temp <- transform_fun(X, label)
  Y <- temp[[1]]
  index <- temp[[2]]
  delta <- temp[[3]]
  tau <- tau_fun(Y)
  sig <- matrix(1, ncol(X), ncol(X))
  for (i in 1:(ncol(X) - 1)) {
    for (j in (i + 1):ncol(X)) {
      temp_label <- label[c(i, j)]
      id1 <- which(index == i)
      id2 <- which(index == j)
      temp <- Y[, c(id1, id2)]
      if (all(temp_label == c("Continuous", "Continuous"))) {
        sig[i, j] <- CC_fun(tau[id1, id2])
      } else if (all(temp_label == c("Continuous", "Binary"))) {
        sig[i, j] <- optimize(BC_fun, c(-0.9999, 0.9999),
                              tau_x = tau[id1, id2],
                              delta_x = delta[j])$minimum
      } else if (all(temp_label == c("Binary", "Continuous"))) {
        sig[i, j] <- optimize(BC_fun, c(-0.9999, 0.9999),
                              tau_x = tau[id1, id2],
                              delta_x = delta[i])$minimum
      } else if (all(temp_label == c("Binary", "Binary"))) {
        sig[i, j] <- optimize(BB_fun, c(-0.9999, 0.9999),
                              tau_x = tau[id1, id2],
                              delta_x = delta[c(i, j)])$minimum
      } else if (all(temp_label == c("Continuous", "Ordinal"))) {
        sig[i, j] <- CO_fun(tau[c(id1, id2), c(id1, id2)], label[c(i, j)], delta[j])
      } else if (all(temp_label == c("Ordinal", "Continuous"))) {
        sig[i, j] <- CO_fun(tau[c(id1, id2), c(id1, id2)], label[c(i, j)], delta[i])
      } else if (all(temp_label == c("Binary", "Ordinal"))) {
        sig[i, j] <- BO_fun(tau[c(id1, id2), c(id1, id2)], label[c(i, j)], delta[c(i, j)])
      } else if (all(temp_label == c("Ordinal", "Binary"))) {
        sig[i, j] <- BO_fun(tau[c(id1, id2), c(id1, id2)], label[c(i, j)], delta[c(i, j)])
      } else if (all(temp_label == c("Ordinal", "Ordinal"))) {
        sig[i, j] <- OO_fun(tau[c(id1, id2), c(id1, id2)], delta[c(i, j)])
      }
      sig[j, i] <- sig[i, j]
    }
  }
  return(sig)
}

# 
