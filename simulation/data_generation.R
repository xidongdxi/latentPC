
# Function to generate a data set with mixed data from a DAG
## Input
# n: Number of observations (rows)
# d: Number of variables (columns)
# deg: Expected number of outgoing edges from each variable
# label: Vector of data types: "continuous", "binary", "ordinal"
## Output
# X: Data set
# A: Adjacency matrix with weighted edges
# Note that the adjacency matrix is transposed, i.e., (i, j) > 0 if there is a
# directed edge from j to i
# sigmahat: Correlation matrix for the latent normal variables (Oracle)
generate_DAG = function(n, d, deg, label, lB = -0.5, uB = 0.5) {
  s <- deg / (d - 1)
  if (length(label) != d) {
    stop("Number of variables and number of labels do not match.")
  }
  A <- randomDAG(d, s, lB = lB, uB = uB)
  X <- rmvDAG(n, A, errDist = "normal")
  sigmahat = cor(X)
  for (j in 1:d) {
    if (label[j] == "continuous") {
      X[, j] <- X[, j]^3
    } else if (label[j] == "binary") {
      temp <- quantile(X[, j], c(0.25, 0.75))
      cut <- runif(1, min = temp[1], max = temp[2])
      X[, j] <- as.numeric(X[, j] > cut)
    } else if (label[j] == "ordinal") {
      temp <- quantile(X[, j], c(0.2, 0.4, 0.6, 0.8))
      cut1 <- runif(1, min = temp[1], max = temp[2])
      cut2 <- runif(1, min = temp[3], max = temp[4])
      X[, j] <- (X[, j] > cut1) + (X[, j] > cut2)
    } else {
      stop('Only include these data types: "continuous", "binary", "ordinal"')
    }
  }
  X <- data.frame(X)
  colnames(X) <- paste0("V", 1:d)
  list(X = X, A = A, sigmahat = sigmahat)
}
