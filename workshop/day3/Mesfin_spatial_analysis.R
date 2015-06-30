
###BreedR workshop (Jaca, Spain, 2015)
#R script: prepared by Dr. Vincent Segura
#R script: to perform spatial and genetic analysis of forest genetic trials using breedR
#the data file "DW.predicted.practice.txt" represents the predicted total dry biomass yield of black poplar experimental design in Orléans (France) in 2011

###Part1
#set working directory
setwd("F:/mngebreselas/Doc-INRA/Workshops/BreedR_Jaca_Spain_2015/Practical_session")
 
#import data into R
my_data<-read.table("DW.predicted.practice.txt", header=T, sep="\t")
 
#check the phenotype data 
summary(my_data$DW.predicted.model1) 
#there are NA's. Total number of NA's = 1704 out of 6745 observations of one variable.

#remove the NA's
data_ok <- subset(my_data, !is.na(DW.predicted.model1))
rm(my_data)

#After removing the NA's, the levels of the variable "Ident" is not changing accordingly. So correct the levels of "Ident". 
#so re-set the factor "Ident"
data_ok$Ident <- as.factor(as.character(data_ok$Ident))

#conversion of coordinates into meters (tree spacing: 1m*2m)
data_ok$X_ok <- data_ok$X * 2
data_ok$Y_ok <- data_ok$Y * 1

#2D map of phenotypic data
library(ggplot2)
p <- qplot(Y_ok, X_ok, fill = DW.predicted.model1, geom = "tile", data = data_ok, main = "DW.predicted.model1",
           xlab = "Y (m)", ylab = "X (m)") + scale_fill_gradient(low = "green", high = "red")

bloc_limits <- cbind(aggregate(data_ok$Y_ok, list(data_ok$Bloc), min)$x,
                     aggregate(data_ok$Y_ok, list(data_ok$Bloc), max)$x,
                     aggregate(data_ok$X_ok, list(data_ok$Bloc), min)$x,
                     aggregate(data_ok$X_ok, list(data_ok$Bloc), max)$x)

rownames(bloc_limits) <- paste("Bloc", 1:6, sep = "")
colnames(bloc_limits) <- c("Y_min", "Y_max", "X_min", "X_max")

p + annotate("rect", xmin = bloc_limits[, 1] - 0.5, xmax = bloc_limits[, 2] + 0.5,
             ymin = bloc_limits[, 3] - 1, ymax = bloc_limits[, 4] + 1, col = "black", alpha = .1)

#histogram:to show distribution of phenotypic data
hist(data_ok$DW.predicted.model1, col = "grey", xlab = "DW.predicted.model1", main = "")
#data are skewed towards lower values. 

###Part2
#fixed effect ANOVA 
linmod <- lm(DW.predicted.model1 ~ as.factor(Bloc) + Ident, data = data_ok)
anova(linmod)

#model QC:
#boxplot per block
boxplot(DW.predicted.model1 ~ as.factor(Bloc), data = data_ok, col = "grey", xlab = "Block", ylab = "DW.predicted.model1")
#suggests within block variance heterogeneity

#look for distribution of residuals
op <- par(mfrow = c(2, 2))
plot(linmod)
par(op)
#model QC shows that residuals are not iid
#it suggests that data should be transformed

#box-cox transformation
library(MASS)
boxcox_transf <- boxcox(linmod)
lambda <- boxcox_transf$x[which.max(boxcox_transf$y)]
lambda
#boxcox procedure suggests the need for data transformation  
#boxcox procedure suggests that the best lambda value = 0.3434343

#box-cox transformation: Y_trsf = (Y^lambda -1)/lambda
data_ok$DW.predicted.model1.trsf <- (data_ok$DW.predicted.model1^lambda-1)/lambda

#check distribution of phenotype data
hist(data_ok$DW.predicted.model1.trsf, col = "grey", xlab = "DW.predicted.model1.trsf", main = "")

#check if the anova model is improved after data transformation:
linmod_trsf <- lm(DW.predicted.model1.trsf ~ as.factor(Bloc) + Ident, data = data_ok)
anova(linmod_trsf)

#boxplot per block for phenotypic data after transformation
boxplot(DW.predicted.model1.trsf ~ as.factor(Bloc), data = data_ok, col = "grey", xlab = "Block", ylab = "DW.predicted.model1.trsf")
#now within block variances seem more homogeneous

#check distribution of residuals
op <- par(mfrow = c(2, 2))
plot(linmod_trsf)
par(op)
#it looks better!

###Part3
#Linear mixed model analysis without a spatial effect using breedR (classical model)
library(breedR)
mixmod_breedR <- remlf90(fixed = DW.predicted.model1.trsf ~ 1 + as.factor(Bloc),
                         random = ~ Ident,
                         data = data_ok,
                         method = "ai")

summary(mixmod_breedR)

#create a new data frame "resid" that includes the data frame "data_ok" and a new variable "resid_breedR"
resid <- data.frame(data_ok, "resid_breedR" = residuals(mixmod_breedR))

##spatial effect diagnosis
#2D plot of residuals from the classical model "mixmod_breedR"
r <- qplot(Y_ok, X_ok, fill = resid_breedR, geom = "tile", data = resid,
           main = "mixmod_breedR_residuals", xlab = "Y (m)", ylab = "X (m)") +
  scale_fill_gradient(low = "green", high = "red")

r + annotate("rect", xmin = bloc_limits[, 1] - 0.5, xmax = bloc_limits[, 2] + 0.5,
             ymin = bloc_limits[, 3] - 1, ymax = bloc_limits[, 4] + 1, col = "black", alpha = .1)

#variograms of residuals from the classical model "mixmod_breedR"
variogram(mixmod_breedR, coord = data_ok[, colnames(data_ok) %in% c("X_ok", "Y_ok")], R = 60)
#it confirms the existence of spatial effects. 

#AIC from the mixed model w/o spatial effects
extractAIC(mixmod_breedR)

###Part 4
#Linear mixed model with a spatial effect: autoregressive with Block effects.
mixmod_breedR_AR1_bloc_grid0 <- remlf90(fixed = DW.predicted.model1.trsf ~ 1 + as.factor(Bloc),
                                        random = ~ Ident,
                                        spatial = list(model = "AR",
                                                       coordinates = data_ok[, c("X_ok", "Y_ok")]),
                                        data = data_ok,
                                        method = "ai")

save(mixmod_breedR_AR1_bloc_grid0, file = "mixmod_breedR_AR1_bloc_grid0_DW.predicted_model1_trsf.Rdata")

#Once the 2 previous commands are executed, one can skip them and replace them by the next one to save time...
load("mixmod_breedR_AR1_bloc_grid0_DW.predicted_model1_trsf.Rdata")
qplot(rho_r, rho_c, fill = loglik, geom = "tile", data = mixmod_breedR_AR1_bloc_grid0$rho)
extractAIC(mixmod_breedR_AR1_bloc_grid0)
summary(mixmod_breedR_AR1_bloc_grid0)

#let's try to refine the grid of rho to improve model fit and thus decrease the AIC
#grid1
rho.grid <- expand.grid(rho_r = seq(0.5, 0.99, length = 4),
                        rho_c = seq(0.5, 0.99, length = 4))

mixmod_breedR_AR1_bloc_grid1 <- remlf90(fixed = DW.predicted.model1.trsf ~ 1 + as.factor(Bloc),
                                        random = ~ Ident,
                                        spatial = list(model = "AR",
                                                       coordinates = data_ok[, c("X_ok", "Y_ok")],
                                                       rho = rho.grid),
                                        data = data_ok,
                                        method = "ai")

save(mixmod_breedR_AR1_bloc_grid1, file = "mixmod_breedR_AR1_bloc_grid1_DW.predicted_model1_trsf.Rdata")

#
load("mixmod_breedR_AR1_bloc_grid1_DW.predicted_model1_trsf.Rdata")
qplot(rho_r, rho_c, fill = loglik, geom = "tile", data = mixmod_breedR_AR1_bloc_grid1$rho)
extractAIC(mixmod_breedR_AR1_bloc_grid1)
summary(mixmod_breedR_AR1_bloc_grid1)

#grid2:refine further
rho.grid <- expand.grid(rho_r = seq(0.81, 0.84, length = 4),
                        rho_c = seq(0.98, 0.99, length = 4))

mixmod_breedR_AR1_bloc_grid2 <- remlf90(fixed = DW.predicted.model1.trsf ~ 1 + as.factor(Bloc),
                                        random = ~ Ident,
                                        spatial = list(model = "AR",
                                                       coordinates = data_ok[, c("X_ok", "Y_ok")],
                                                       rho = rho.grid),
                                        data = data_ok,
                                        method = "ai")

save(mixmod_breedR_AR1_bloc_grid2, file = "mixmod_breedR_AR1_bloc_grid2_DW.predicted_model1_trsf.Rdata")

#
load("mixmod_breedR_AR1_bloc_grid2_DW.predicted_model1_trsf.Rdata")
qplot(rho_r, rho_c, fill = loglik, geom = "tile", data = mixmod_breedR_AR1_bloc_grid2$rho)
extractAIC(mixmod_breedR_AR1_bloc_grid2)
summary(mixmod_breedR_AR1_bloc_grid2)

#grid3
rho.grid <- expand.grid(rho_r = seq(0.83, 0.85, length = 4),
                        rho_c = seq(0.97, 0.99, length = 4))

mixmod_breedR_AR1_bloc_grid3 <- remlf90(fixed = DW.predicted.model1.trsf ~ 1 + as.factor(Bloc),
                                        random = ~ Ident,
                                        spatial = list(model = "AR",
                                                       coordinates = data_ok[, c("X_ok", "Y_ok")],
                                                       rho = rho.grid),
                                        data = data_ok,
                                        method = "ai")

save(mixmod_breedR_AR1_bloc_grid3, file = "mixmod_breedR_AR1_bloc_grid3_DW.predicted_model1_trsf.Rdata")

#
load("mixmod_breedR_AR1_bloc_grid3_DW.predicted_model1_trsf.Rdata")
qplot(rho_r, rho_c, fill = loglik, geom = "tile", data = mixmod_breedR_AR1_bloc_grid3$rho)
extractAIC(mixmod_breedR_AR1_bloc_grid3)
summary(mixmod_breedR_AR1_bloc_grid3)

#grid4
rho.grid <- expand.grid(rho_r = seq(0.84, 0.88, length = 4),
                        rho_c = seq(0.96, 0.99, length = 4))

mixmod_breedR_AR1_bloc_grid4 <- remlf90(fixed = DW.predicted.model1.trsf ~ 1 + as.factor(Bloc),
                                        random = ~ Ident,
                                        spatial = list(model = "AR",
                                                       coordinates = data_ok[, c("X_ok", "Y_ok")],
                                                       rho = rho.grid),
                                        data = data_ok,
                                        method = "ai")

save(mixmod_breedR_AR1_bloc_grid4, file = "mixmod_breedR_AR1_bloc_grid4_DW.predicted_model1_trsf.Rdata")

#
load("mixmod_breedR_AR1_bloc_grid4_DW.predicted_model1_trsf.Rdata")
qplot(rho_r, rho_c, fill = loglik, geom = "tile", data = mixmod_breedR_AR1_bloc_grid4$rho)
extractAIC(mixmod_breedR_AR1_bloc_grid4)
summary(mixmod_breedR_AR1_bloc_grid4)

#grid5: 
rho.grid <- expand.grid(rho_r = seq(0.83, 0.88, length = 4),
                        rho_c = seq(0.95, 0.99, length = 4))

mixmod_breedR_AR1_bloc_grid5 <- remlf90(fixed = DW.predicted.model1.trsf ~ 1 + as.factor(Bloc),
                                        random = ~ Ident,
                                        spatial = list(model = "AR",
                                                       coordinates = data_ok[, c("X_ok", "Y_ok")],
                                                       rho = rho.grid),
                                        data = data_ok,
                                        method = "ai")

save(mixmod_breedR_AR1_bloc_grid5, file = "mixmod_breedR_AR1_bloc_grid5_DW.predicted_model1_trsf.Rdata")

#
load("mixmod_breedR_AR1_bloc_grid5_DW.predicted_model1_trsf.Rdata")
qplot(rho_r, rho_c, fill = loglik, geom = "tile", data = mixmod_breedR_AR1_bloc_grid5$rho)
extractAIC(mixmod_breedR_AR1_bloc_grid5)
summary(mixmod_breedR_AR1_bloc_grid5)

#compare model AICs
extractAIC(mixmod_breedR)

extractAIC(mixmod_breedR_AR1_bloc_grid0)
extractAIC(mixmod_breedR_AR1_bloc_grid1)
extractAIC(mixmod_breedR_AR1_bloc_grid2)
extractAIC(mixmod_breedR_AR1_bloc_grid3)
extractAIC(mixmod_breedR_AR1_bloc_grid4)
extractAIC(mixmod_breedR_AR1_bloc_grid5)

#mixmod_breedR_AR1_bloc_grid5:is the best model(minimum AIC)
#best autoregressive parameters for rows and columns: (0.846, 0.976)

##mixed model with a spatial effect with the best autoregressive parameters for rows and columns: autoregressive with block effect
selmod <- remlf90(fixed = DW.predicted.model1.trsf ~ 1 + as.factor(Bloc),
                  random = ~ Ident,
                  spatial = list(model = "AR",
                                 coordinates = data_ok[, c("X_ok", "Y_ok")],
                                 rho = c(rho_r= 0.846, rho_c = 0.976)),
                  data = data_ok,
                  method = "ai")

save(selmod, file = "selmod_DW.predicted_model1_trsf.Rdata")

#
load("selmod_DW.predicted_model1_trsf.Rdata")
summary(selmod)


##Part 5: Focus on selected spatial model "selmod"
#extraction of predicted values for fixed and random effects
BLUP_genot <- data.frame("Ident" = rownames(selmod$ranef$Ident),
                         GBLUP = selmod$ranef$Ident$value)

Output <- data.frame(na.omit(data_ok, "Coord" = paste(data_ok$X_ok, data_ok$Y_ok, sep = "-")))

Output$Coord <- paste(Output$X_ok, Output$Y_ok, sep = "-")

Block_effect <- as.vector(sapply(selmod$fixed, function(x){x$value}))
Output$Block_effect <- Block_effect[1]
for (i in 2:length(Block_effect)){
  Output$Block_effect[Output$Bloc == i] <- Block_effect[i]
}

Output <- merge(Output, BLUP_genot, by = "Ident")

spatial <- data.frame(selmod$effects$spatial$effects$ar$coordinates,
                      Spatial = selmod$effects$spatial$effects$ar$incidence.matrix %*%
                        ranef(selmod)$spatial)

Output <- merge(Output,
                data.frame("Coord" = paste(spatial$X_ok, spatial$Y_ok, sep = "-"),
                           "Spatial" = spatial$Spatial),
                by = "Coord")

Output$Fitted <- Output$Block_effect + Output$GBLUP + Output$Spatial

resid$Selmod <- na.omit(residuals(selmod))

Output <- merge(Output,
                data.frame("Coord" = paste(resid$X_ok, resid$Y_ok, sep = "-"),
                           "Resid" = resid$Selmod),
                by = "Coord")

Output$Phen_Adj <- mean(Output$DW.predicted.model1.trsf) + Output$GBLUP + Output$Resid

##Save the predicted values for further analysis
save(list = c("Output", "BLUP_genot"), file = "Output_BreedR_Orl_DW.predicted.model1.trsf.Rdata")
write.table(Output, file = "Output_BreedR_Orl_DW.predicted.model1.trsf.txt", sep = "\t", row.names = FALSE)

#2D map of spatial effects
s <- qplot(Output$Y_ok, Output$X_ok, fill = Output$Spatial, geom = "tile", main ="Spatial effects",
           xlab = "Y (m)", ylab = "X (m)") + scale_fill_gradient(low = "green", high = "red")
s + annotate("rect", xmin = bloc_limits[, 1] - 0.5, xmax = bloc_limits[, 2] + 0.5,
             ymin = bloc_limits[, 3] - 1, ymax = bloc_limits[, 4] + 1, col = "black", alpha = .1)

#2D map of residuals from selected model "selmod"
t <- qplot(Output$Y_ok, Output$X_ok, fill = Output$Resid, geom = "tile", data = resid, main ="residuals",
           xlab = "Y (m)", ylab = "X (m)") + scale_fill_gradient(low = "green", high = "red")
t + annotate("rect", xmin = bloc_limits[, 1] - 0.5, xmax = bloc_limits[, 2] + 0.5,
             ymin = bloc_limits[, 3] - 1, ymax = bloc_limits[, 4] + 1, col = "black", alpha = .1)

#2D map of genotype BLUPs from selected model "selmod"
u <- qplot(Output$Y_ok, Output$X_ok, fill = Output$GBLUP, geom = "tile",
           main ="Genotype BLUPs", xlab = "Y (m)", ylab = "X (m)") +
  scale_fill_gradient(low = "green", high = "red")
u + annotate("rect", xmin = bloc_limits[, 1] - 0.5, xmax = bloc_limits[, 2] + 0.5,
             ymin = bloc_limits[, 3] - 1, ymax = bloc_limits[, 4] + 1, col = "black", alpha = .1)

#check distribution of genotype BLUPs per Bloc
boxplot(GBLUP ~ Bloc, data = Output, col = "grey", xlab = "Block", ylab = "Genotype BLUP")

#compare variograms of residuals from classical model (mixmod_breedR) & selected spatial model(selmod)
#plot the 2 variograms on the same scale using the variogram function of breedR
variogram_null <- variogram(mixmod_breedR, coord = data_ok[, colnames(data_ok) %in% c("X_ok", "Y_ok")], R = 60)
variogram_selmod <- variogram(selmod, R = 60)

# Extract both isotropic variograms to be compared
isotropic_variograms <- rbind(cbind(model = "Block", variogram_null[["isotropic"]]),
                              cbind(model = "AR1_Block", variogram_selmod[["isotropic"]]))

# Then, you can plot the 2 variograms together under the same scale
ggplot(isotropic_variograms, aes(distance, variogram)) +
  geom_point() +
  geom_line() +
  stat_smooth(se = FALSE, method = 'auto') + 
  facet_wrap(~ model)

#histogram of random effect blups from selected spatial model
op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 1))
hist(BLUP_genot$GBLUP, main = "Genotype", xlab = "BLUPs", col = "grey")
hist(Output$Spatial, main = "Spatial effect", xlab = "BLUPs", col = "grey")
hist(Output$Resid, col = "grey", "xlab" = "Resid", main = "residuals")
hist(Output$Phen_Adj, col = "grey", "xlab" = "Phenotype", main = "Adjusted phenotype")
par(op)

##Part 6: Genetic analysis based on the selected model "selmod" 
#extraction of the average information (AI) matrix from breedR output object "selmod"
aimat <- which(selmod$reml$output == " inverse of AI matrix (Sampling Variance)")

varcov_mat_breedR <- matrix(na.omit(as.numeric(unlist(apply(data.frame(selmod$reml$output[
  (aimat + 1):(aimat + 3)]), 1, function(x){strsplit(x, " ")})))), 3, 3)
colnames(varcov_mat_breedR) <- c("Ident", "spatial", "Residual")
rownames(varcov_mat_breedR) <- c("Ident", "spatial", "Residual")

#delta method for H²ind 
library(msm)
H2 <- selmod$var["Ident", "Estimated variances"] /
  (selmod$var["Ident", "Estimated variances"] + selmod$var["Residual", "Estimated variances"])
se_H2 <- deltamethod(~ x1 / (x1 + x2),
                       + c(selmod$var["Ident", "Estimated variances"],
                           + selmod$var["Residual", "Estimated variances"]),
                       + varcov_mat_breedR[c("Ident", "Residual"), c("Ident", "Residual")])
round(c(H2, 1.96*se_H2), 2)

#H²mean with se
nb_obs_genot <- length(na.omit(selmod$mf$DW.predicted.model1.trsf)) / length(levels(selmod$mf$Ident))
H2_mean <- selmod$var["Ident", "Estimated variances"] /
  (selmod$var["Ident", "Estimated variances"] +
     selmod$var["Residual", "Estimated variances"] / nb_obs_genot)
se_H2_mean <- deltamethod(~ x1 / (x1 + x3 / nb_obs_genot),
                          selmod$var[, "Estimated variances"],
                          varcov_mat_breedR)
round(c(H2_mean, 1.96*se_H2_mean), 2)

#coefficient of genetic variation (CV)
CVg <- sqrt(selmod$var["Ident", "Estimated variances"]) /
  mean(selmod$mf$DW.predicted.model1.trsf, na.rm = TRUE) * 100
round(CVg, 2)

##########        ###########            ##############
#Tips!
#stacking the AICs into a vector
mixmods<-list(mixmod_breedR, selmod)
names(mixmods)<-c("Block", "AR1_Block")
AICs<-sapply(mixmods,extractAIC)

#classical model vs selected spatial model: estimation of individual broad-sense heritability with each model
herit <- sapply(mixmods, function(x){x$var["Ident", "Estimated variances"] /
                                       (x$var["Ident", "Estimated variances"]
                                        + x$var["Residual", "Estimated variances"])})

##############Hope you will enjoy breedR!!!
 
 