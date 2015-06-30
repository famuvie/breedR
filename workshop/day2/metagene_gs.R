library(breedR)
library(pedigree)
library(pedigreemm)
library(pedantics)
library(car)
library(hexbin)

# import pedigree file already including founders
pheno_ped <- read.table("pheno_ped_case2.txt", header = T)
names(pheno_ped) <- c("self", "dad", "mum", "gen", "BV_X", "phe_X")
pheno_ped4 <- subset(pheno_ped, gen==4)
pedig <- build_pedigree(c('self', 'dad', 'mum'), data = pheno_ped)
num_tot <- as.numeric(dim(pheno_ped)[1])
num_rec <- as.numeric(dim(pheno_ped4)[1])
founders <- num_tot-num_rec

# check variables in data.frame
summary(pheno_ped)

# plot distributions for traits
colores <- rainbow(5)
plot(density(pheno_ped4$phe_X), col=colores[1], main="phe_X")

# run genetic model for each trait
pedmod_phe_X <- remlf90(fixed  = phe_X ~ 1,
                      genetic = list(model    = 'add_animal',
                                     pedigree = pedig,
                                     id       = 'self'),
                      data   = pheno_ped4)
summary(pedmod_phe_X)
hh_phe_X <- pedmod_phe_X$var[1,1]/(pedmod_phe_X$var[1,1]+pedmod_phe_X$var[2,1])

# # do generic model
# # recover incidence matrix
# ped.inc <- model.matrix(pedmod_phe_X)$genetic
# # alternative way of building up the incidence matrix
# ped.inc.alt0 <- matrix(0,num_tot,1)
# ped.inc.alt1 <- diag(1,num_tot,num_tot)
# # Make it index matrix, stored as 1-based integer index vectors.
# # An index matrix is a matrix with exactly one non-zero entry per row.
# # Index matrices are useful for mapping observations 
# ped.inc.alt <- as(cbind(ped.inc.alt0,ped.inc.alt1),"indMatrix")
# # recover var/covar matrix
# ped.cov <- pedigreemm::getA(pedig)
# # run model with given incidences and structure matrices
# pedmod_generic_phe_X <- remlf90(fixed   = phe_X ~ 1,
#                               generic = list(gen = list(ped.inc,ped.cov)),
#                               data    = pheno_ped4)
# pedmod_generic_alt_phe_X <- remlf90(fixed   = phe_X ~ 1,
#                                   generic = list(gen = list(ped.inc.alt,ped.cov)),
#                                   data    = pheno_ped4)
# summary(pedmod_generic_phe_X)
# summary(pedmod_generic_alt_phe_X)

# import genotypic data
genot <- read.table("genotypes_case2.txt", sep = " ", header = F)
# removes self from file, follows same order as production file
genot_no_self <- genot[c(-1)]
# number of markers
n_SNP <- as.numeric(dim.data.frame(genot_no_self)[2])
# matrix of genotypes
X <- as.matrix(genot_no_self,nrow=num_rec,ncol=n_SNP,dimnames=NULL)
# calculates frequency of favourable allele per marker
Pi <- apply(X,2,sum)/(2*num_rec)
mat_Pi <- matrix(rep(Pi,num_rec),ncol=n_SNP,byrow=T)
W <- matrix(0,nrow=num_rec,ncol=n_SNP)
W <- X - (2*mat_Pi)
het <- 2*sum(Pi*(1-Pi))
G <- W%*%t(W) / het
# inverse is not needed
G_inv <- solve(G+diag(num_rec)*0.01)

# Make it index matrix, stored as 1-based integer index vectors.
# An index matrix is a matrix with exactly one non-zero entry per row.
# Index matrices are useful for mapping observations 
gen.inc.alt0 <- matrix(0,num_rec,1)
gen.inc.alt1 <- diag(1,num_rec,num_rec)
gen.inc.alt <- as(cbind(gen.inc.alt1),"indMatrix")
# Makes G symmetric, sparse numeric matrix in compressed, column-oriented format
gen.cov <- as(G, "dsCMatrix")

gsmod_phe_X <- remlf90(fixed   = phe_X ~ 1,
                     generic = list(gen = list(gen.inc.alt,gen.cov)),
                     data    = pheno_ped4)
summary(gsmod_phe_X)
hh_gs_phe_X <- gsmod_phe_X$var[1,1]/(gsmod_phe_X$var[1,1]+gsmod_phe_X$var[2,1])
# calculates correlations with true BV
PBV.full.gen <- ranef(gsmod_phe_X)$gen
PBV.gen <- model.matrix(gsmod_phe_X)$gen %*% PBV.full.gen

PBV.full.ped <- ranef(pedmod_phe_X)$genetic
PBV.ped <- model.matrix(pedmod_phe_X)$genetic %*% PBV.full.ped

cor_gen <- cor(pheno_ped4$BV_X,PBV.gen)
cor_ped <- cor(pheno_ped4$BV_X,PBV.ped)
max_BV <- max(pheno_ped4$BV_X)
min_BV <- min(pheno_ped4$BV_X)
max_EBV <- max(PBV.gen,PBV.ped)
min_EBV <- min(PBV.gen,PBV.ped)
# plots scatterplots showing correlated patterns
tiff(filename="scatterplot_TBV_vs_EBV_high_hh.tif",width=500,height=500,compression = c("none"))
correl <- round(cor_gen, digits = 4)
heading <- paste("corr =",correl," [TBV-EBV]")
scatterplot(pheno_ped4$BV_X,PBV.gen, xlab = "True_BV", ylab = "EBV",
            main = heading, xlim = c(min_BV,max_BV), ylim = c(min_EBV,max_EBV))
dev.off()

tiff(filename="scatterplot_TBV_vs_EBV_low_hh.tif",width=500,height=500,compression = c("none"))
correl <- round(cor_ped, digits = 4)
heading <- paste("corr =",correl," [TBV-EBV]")
scatterplot(pheno_ped4$BV_X,PBV.ped, xlab = "True_BV", ylab = "EBV",
            main = heading, xlim = c(min_BV,max_BV), ylim = c(min_EBV,max_EBV))
dev.off()

# does hitmaps of both A and G matrices
tiff(filename="heatmap_G.tif",width=500,height=500,compression = c("none"))
  heatmap(G)
dev.off()
A <- as(getA(pedig),"matrix")
tiff(filename="heatmap_A.tif",width=500,height=500,compression = c("none"))
  heatmap(A[97:672,97:672])
dev.off()


