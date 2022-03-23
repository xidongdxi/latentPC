
# Modify the randomDAG function from the pcalg package
## Input
# n: Number of observations (rows)
# prob: Probability of an outgoing edge from each variable to another variable
# lB, uB: Lower and upper bounds from which an edge weight is drawn from
# the uniform distribution
# V: Names of variables (columns)
## Output
# A DAG with weighted edges which preserves the topological order, i.e.,
# an edge can only go from a lower order variable (e.g., H1) to
# a higher order variable (e.g., H2)
randomDAG2 <- function (n, prob, lB = 0.5, uB = 1, V = as.character(1:n)) {
  stopifnot(n >= 2, is.numeric(prob), length(prob) == 1, 0 <= 
              prob, prob <= 1, is.numeric(lB), is.numeric(uB), lB <= 
              uB)
  edL <- vector("list", n)
  nmbEdges <- 0L
  for (i in seq_len(n - 2)) {
    listSize <- rbinom(1, n - i, prob)
    nmbEdges <- nmbEdges + listSize
    edgeList <- sample(seq(i + 1, n), size = listSize)
    y <- runif(length(edgeList), min = 0, max = 2 * (uB - lB))
    weightList <- sapply(y, function(x) ifelse(x < uB - lB, x - uB, x - uB + 2 * lB))
    edL[[i]] <- list(edges = edgeList, weights = weightList)
  }
  listSize <- rbinom(1, 1, prob)
  if (listSize > 0) {
    nmbEdges <- nmbEdges + 1
    edgeList <- n
    y <- runif(1, min = 0, max = 2 * (uB - lB))
    weightList <- ifelse(y < uB - lB,
                         y - uB, y - uB + 2 * lB)
  }
  else {
    edgeList <- integer(0)
    weightList <- numeric(0)
  }
  edL[[n - 1]] <- list(edges = edgeList, weights = weightList)
  if (nmbEdges > 0) {
    edL[[n]] <- list(edges = integer(0), weights = numeric(0))
    names(edL) <- V
    new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
  }
  else new("graphNEL", nodes = V, edgemode = "directed")
}

# Function to generate a data set with mixed data from a DAG
## Input
# n: Number of observations (rows)
# d: Number of variables (columns)
# deg: Expected number of outgoing edges from each variable
# label: Vector of data types: "Continuous", "Binary", "Ordinal"
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
    if (label[j] == "Continuous") {
      X[, j] <- X[, j]^3
    } else if (label[j] == "Binary") {
      temp <- quantile(X[, j], c(0.25, 0.75))
      cut <- runif(1, min = temp[1], max = temp[2])
      X[, j] <- as.numeric(X[, j] > cut)
    } else if (label[j] == "Ordinal") {
      temp <- quantile(X[, j], c(0.2, 0.4, 0.6, 0.8))
      cut1 <- runif(1, min = temp[1], max = temp[2])
      cut2 <- runif(1, min = temp[3], max = temp[4])
      X[, j] <- (X[, j] > cut1) + (X[, j] > cut2)
    } else {
      stop('Only include these data types: "Continuous", "Binary", "Ordinal"')
    }
  }
  X <- data.frame(X)
  colnames(X) <- paste0("V", 1:d)
  list(X = X, A = A, sigmahat = sigmahat)
}
