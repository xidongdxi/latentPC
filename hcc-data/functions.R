library(pbivnorm)
library(rootSolve)
library(pcalg)
library(Matrix)

camel <- function(x){ #function for camel case
  capit <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2, nchar(x)))
  x1 = sapply(strsplit(x, "\\."), function(x) paste(capit(x), collapse=""))
  sapply(strsplit(x1, " "), function(x) paste(capit(x), collapse=""))
}
### Tunning parameter selection with the StARS method.
StARS = function(X, a1 = seq(0.01, 0.1, by = 0.01), TunModel = "LGC", label, ratio, rep.times)
{
  l1 = length(a1)
  n = dim(X)[1]
  p = dim(X)[2]
  b = round(n*ratio)
  D = rep(0, l1)
  N = rep.times
  for(i in 1:l1)
  {
    cat(a1[i], "\r")
    all.edges = matrix(0, N, p*p)
    for(j in 1:N)
    {
      cat(j, "\r")
      ind1 = sample(1:n, size = b, replace = F)
      X1 = X[ind1,]
      
      if(TunModel == "LGC") {
        sig1 <- latent_pc(X1, label)
        sig1 <- as.matrix(nearPD(sig1, corr = TRUE, maxit = 10000)$mat)
        }
      if(TunModel == "Gaussian") sig1 = cor(X1)
      if(TunModel == "Rank") sig1 = sin(cor(X1, method = "kendall")*pi/2)
      
      pc.fit1 <- pc(suffStat = list(C = sig1, n = b), indepTest = gaussCItest, alpha=a1[i], 
                    labels = colnames(X), skel.method = "stable",verbose = F)
      t1 = as.numeric(as(pc.fit1, "amat"))
      all.edges[j,] = t1
    }
    thetahat = apply(all.edges, 2, function(s){2*mean(s)*(1-mean(s))})
    D[i] = mean(thetahat, na.rm = T)
  }
  list(alpha = a1, D = D, cumD = cummax(D))
}

# est = approx(cummax(D), a1, xout=0.05)
