library(pcalg)
library(mvtnorm)
library(Matrix)
library(sbgcop)
library(MXM)
library(igraph)
library(readxl)
library(xtable)
source("/hcc-data/real_functions.R")

# Read data
p <- read.csv("/hcc-data/data/hcc-data.txt", sep = ",", header = FALSE)
p[p == "?"] <- NA
n1 <- read.table("/hcc-data/data/HCCnames.txt", sep = ":", header = FALSE)
c1 <- read.csv("/hcc-data/data/colnames.csv", header = F)
colnames(c1) <- c("Names", "Type", "Range", "Mean or mode", "Missingness")
c1$Names[1] <- "Gender"
colnames(p) <- c(c1[, 1], "class")
colnames(p) <- camel(colnames(p))
c1 <- rbind(c1, c("Class", "binary", "0/1", 1, 0))
mm <- c1[, 4]

# Label
label.p <- n1[,2]
label.p[50] <- " nominal"
label.p <- gsub(" ", "", label.p)
label.p[label.p == "nominal"] <- "binary"
label.p[label.p == "integer"] <- "continuous"

### Filter out the columns with many missing values
ind <- which(colMeans(is.na(p)) < 0.05)
p1 <- p[, ind]
la <- label.p[ind]
c1$Fullname <- n1[, 1]
descp <- c1[ind, ]
descp$Label <- la

p2 <- na.omit(p1)
p3 <- lapply(p2, function(x) {
  if(is.character(x)) as.numeric(x) else x
})
p3 <- as.data.frame(p3)
descp$Abbreviation <- abbreviate(colnames(p3), 6)
res <- xtable(descp[, c(6, 8, 7, 3, 4)])
# Description of the hepatocellular carcinoma data set, where n = 138, d = 28.
# print(res, include.rownames = FALSE)

# Re-categorize ordinal variables to 0, 1, 2, ... 
p3[, res$Abbreviation == "Encflp"] <- as.numeric(p3[, res$Abbreviation == "Encflp"]) - 1
p3[, res$Abbreviation == "Ascits"] <- as.numeric(p3[, res$Abbreviation == "Ascits"]) - 1
newhcc = list(hcc = p3, label = la)

# Tuning parameter
X = as.matrix(newhcc$hcc)
label = newhcc$label
alpha <- seq(0.01, 0.1, by = 0.01)
ratio = 0.7
rep.times = 100
scen <- expand.grid(alpha, ratio, rep.times)
colnames(scen) <- c("alpha", "ratio", "rep.times")
scen <- list(scen, X, label)

library(future.apply)
source("utility.R")
source("copula_pc.R")
source("latent_pc.R")
source("/hcc-data/real_functions.R")
plan(multisession)
start <- proc.time()
result <- future_lapply(1:nrow(scen[[1]]), FUN = StARS_parallel, future.seed = NULL,
                        future.packages = c("pcalg", "Matrix", "mvtnorm", "pcaPP"),
                        scen = scen)
c(proc.time() - start)[3]
aaa <- as.data.frame(do.call(rbind, result))
alpha[tail(which(diff(aaa$D) > 0), 1) - 1] # 0.07

# Latent PC
sig1 <- latent_pc(X, label)
sig1 <- as.matrix(nearPD(sig1, corr = TRUE, maxit = 10000)$mat)
pc.fit <- pc(suffStat = list(C = sig1, n = nrow(X)),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha = 0.07,
             labels = abbreviate(colnames(X), 6), skel.method = "stable.fast",
             verbose = F)
gg = attributes(pc.fit)$graph
dd = Rgraphviz::getDefaultAttrs()
dd$node$fontsize = "20"
pdf("/hcc-data/results/HCC_latent_pc.pdf", height = 20, width = 15)
par(mar = c(0.01, 0.01, 0.01, 0.01))
Rgraphviz::plot(gg, attrs = dd)
dev.off()

# Vanilla PC
sig <- cor(X)
pc.fit <- pc(suffStat = list(C = sig, n = nrow(X)),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha = 0.07,
             labels = abbreviate(colnames(X), 6), skel.method = "stable.fast",
             verbose = F)
gg = attributes(pc.fit)$graph
dd = Rgraphviz::getDefaultAttrs()
dd$node$fontsize = "20"
pdf("/hcc-data/results/HCC_vanilla_pc.pdf", height = 20, width = 15)
par(mar = c(0.01, 0.01, 0.01, 0.01))
Rgraphviz::plot(gg, attrs = dd)
dev.off()

# Rank PC
sig <- sin(pi / 2 * pcaPP::cor.fk(X))
pc.fit <- pc(suffStat = list(C = sig, n = nrow(X)),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha = 0.07,
             labels = abbreviate(colnames(X), 6), skel.method = "stable.fast",
             verbose = F)
gg = attributes(pc.fit)$graph
dd = Rgraphviz::getDefaultAttrs()
dd$node$fontsize = "20"
pdf("/hcc-data/results/HCC_rank_pc.pdf", height = 20, width = 15)
par(mar = c(0.01, 0.01, 0.01, 0.01))
Rgraphviz::plot(gg, attrs = dd)
dev.off()

# Copula PC
cop.obj <- inferCopulaModel(X, nsamp = 1000, S0 = diag(ncol(X)) / nrow(X), verb = F)
C_samples <- cop.obj$C.psamp[, , 501:1000]
corr.cop <- apply(C_samples, c(1, 2), mean)
pc.fit <- pc(suffStat = list(C = corr.cop, n = nrow(X)),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha = 0.07,
             labels = abbreviate(colnames(X), 6), skel.method = "stable.fast",
             verbose = F)
gg = attributes(pc.fit)$graph
dd = Rgraphviz::getDefaultAttrs()
dd$node$fontsize = "20"
pdf("/hcc-data/results/HCC_copula_pc.pdf", height = 20, width = 15)
par(mar = c(0.01, 0.01, 0.01, 0.01))
Rgraphviz::plot(gg, attrs = dd)
dev.off()

# MM-PC
skel <- suppressWarnings(pc.skel(as.data.frame(X), method = "comb.mm", alpha = 0.07))
pc.fit <- pc.or(skel)
tsig <- pc.fit$G
tsig[tsig == 2] <- 0
tsig[tsig == 3] <- 1
edL <- vector("list", length = ncol(X))
names(edL) <- abbreviate(colnames(X), 6)
for (i in 1:ncol(X)) {
  temp <- tsig[, i]
  if (any(temp == 1)) {
    edL[[i]] <- list(edges = as.numeric(which(temp == 1)), weights = rep(1, sum(temp == 1)))
  } else {
    edL[[i]] <- list()
  }
}
gg <- new("graphNEL", nodes = abbreviate(colnames(X), 6), edgeL = edL, edgemode = "directed")
dd = Rgraphviz::getDefaultAttrs()
dd$node$fontsize = "20"
pdf("/hcc-data/results/HCC_mm-pc.pdf", height = 20, width = 15)
par(mar = c(0.01, 0.01, 0.01, 0.01))
Rgraphviz::plot(gg, attrs = dd)
dev.off()

