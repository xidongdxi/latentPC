rm(list= ls())
library(igraph)
library(readxl)
setwd("C:/Users/dxi1/OneDrive - Gilead Sciences/Paper/Latent PC/Code/Update/shd/hcc")
library(xtable)
source("real_functions.R")


p <- read.csv("hcc-data.txt", sep = ",", header = FALSE)
p[p=="?"] = NA
n1 <- read.table("HCCnames.txt", sep = ":", header = FALSE)
c1 <- read.csv("colnames.csv", header = F)
colnames(c1) <- c("Names", "Type", "Range", "Mean or mode", "Missingness")
c1$Names[1] <- "Gender"
colnames(p) <- c(c1[,1], "class")
colnames(p) <- camel(colnames(p))
c1 <- rbind(c1, c("Class", "binary", "0/1", 1, 0))
c1$Range <- gsub("-", "--", c1$Range)
mm <- c1[,4]

label.p <- n1[,2]
label.p[50] <- " nominal"
label.p <- gsub(" ", "", label.p)
label.p[label.p == "nominal"] <- "binary"
label.p[label.p == "integer"] <- "continuous"

p1 <- lapply(p, function(x) {
  if(is.character(x)) as.numeric(x) else x
})
p1 <- as.data.frame(p1)
for(j in 1:49) {
  p1[is.na(p1[, j]), j] <- mm[j]
}

### Filter out the columns with many missing values
ind1 <- which(colMeans(is.na(p)) < 0.05)
p2 <- p[, ind1]
la <- label.p[ind1]
c1$Fullname <- n1[, 1]
descp <- c1[ind1, ]
descp$Label <- la

p2 <- na.omit(p2)
p3 <- lapply(p2, function(x) {
  if(is.character(x)) as.numeric(x) else x
})
p3 <- as.data.frame(p3)
descp$Abbreviation <- abbreviate(colnames(p3), 6)
res <- xtable(descp[, c(6,8, 7, 3, 4)])
print(res, include.rownames = FALSE)

newhcc = list(hcc = p3, label = la)
save(newhcc, file = "hcc.rda")

# Tuning parameter
X = as.matrix(newhcc$hcc)
label = newhcc$label
r1 = 0.7
rep.t = 100
ss = StARS(X, label = label, ratio = r1, rep.times = rep.t)
filename1 = paste("HCC_tune_", r1, ".rda", sep = "")
save(ss, file = filename1)
filename2 = paste("HCC_tune_plot_", r1, ".pdf", sep = "")
pdf(filename2)
plot(ss$alpha, ss$D)
dev.off()

# Latent PC
labels <- label_fun(p3)
sig1 <- latent_pc(p3, labels)
sig1 <- as.matrix(nearPD(sig1, corr = TRUE, maxit = 10000)$mat)
pc.fit <- pc(suffStat = list(C = sig1, n = nrow(p3)),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha = alpha, labels = colnames(X), skel.method = "stable.fast",
             verbose = F)

# Vanilla PC
sig <- cor(p3)
pc.fit <- pc(suffStat = list(C = sig, n = nrow(p3)),
             indepTest = gaussCItest, ## indep.test: partial correlations
             alpha = alpha, labels = colnames(X), skel.method = "stable.fast",
             verbose = F)



sig1 = LGC_sigma_AllType(as.matrix(p3), la)
#sig1[which(is.na(sig1), arr.ind = T)] = cor(p)[which(is.na(sig1), arr.ind = T)]
sum(is.na(sig1))
sig1 = nearPD(sig1, corr = TRUE)$mat
pc.fit1 <- pc(suffStat = list(C = sig1, n = nrow(p3)), indepTest = gaussCItest, maj.rule = T, 
              alpha=0.05, labels = abbreviate(colnames(p3), 6), skel.method = "stable", verbose = F)
#summary(pc.fit1)
adj1 = t(as(pc.fit1, "amat"))
sum(adj1)
gg = attributes(pc.fit1)$graph
dd = getDefaultAttrs()
dd$node$fontsize = "20"
pdf("/Users/zhanruicai/Library/Mobile Documents/com~apple~CloudDocs/Research/MixedCausalDoc/HCC_latent.pdf",  height = 20, width = 15)
par(mar = c(0.01, 0.01, 0.01, 0.01))
plot(gg, attrs = dd)
dev.off()



sig2 = cor(p3, method = "kendall")
sig3 = cor(p3)
pc.fit2 <- pc(suffStat = list(C = sig2, n = nrow(p3)), indepTest = gaussCItest, 
              alpha=0.07, labels = abbreviate(colnames(p3), 6), skel.method = "stable",verbose = F)
g2 = attributes(pc.fit2)$graph
adj2 = t(as(pc.fit2, "amat"))
sum(adj2)
pdf("/Users/zhanruicai/Library/Mobile Documents/com~apple~CloudDocs/Research/MixedCausalDoc/HCC_rank.pdf",  height = 20, width = 15)
par(mar = c(0.01, 0.01, 0.01, 0.01))
plot(g2, attrs = dd)
dev.off()

pc.fit3 <- pc(suffStat = list(C = sig3, n = nrow(p3)), indepTest = gaussCItest, 
              alpha=0.07, labels = abbreviate(colnames(p3), 6), skel.method = "stable",verbose = F)
adj3 = t(as(pc.fit3, "amat"))
sum(adj3)
g3 = attributes(pc.fit3)$graph
pdf("/Users/zhanruicai/Library/Mobile Documents/com~apple~CloudDocs/Research/MixedCausalDoc/HCC_vani.pdf",  height = 20, width = 15)
par(mar = c(0.01, 0.01, 0.01, 0.01))
plot(g3, attrs = dd)
dev.off()






