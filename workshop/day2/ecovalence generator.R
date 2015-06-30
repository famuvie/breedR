library(breedR)

# family effects
num_fam <- 35
mean_eff_fam <- 0
sd_eff_fam <- 0.15
fam_eff <- rnorm(num_fam,mean_eff_fam,sd_eff_fam)      #generate some random normal deviates

# family size and residuals
fam_size <- 20
mean_residuals <- 0
sd_residuals <- 0.55

# site effects
num_sites <- 5
mean_eff_site <- 0
site_eff <- rnorm(num_sites,mean_eff_site,0.5)

# global_mean
global_mean <- 0

# family by site interaction
mean_eff_inter <- 0
sd_eff_inter <- 1-sd_eff_fam-sd_residuals

total_obs <- num_fam*num_sites*fam_size

self <- matrix(data = NA, nrow = total_obs, ncol = 1, dimnames = NULL)
fam <- matrix(data = NA, nrow = total_obs, ncol = 1, dimnames = NULL)
site <- matrix(data = NA, nrow = total_obs, ncol = 1, dimnames = NULL)
inter <- matrix(data = NA, nrow = total_obs, ncol = 1, dimnames = NULL)
rep <- matrix(data = NA, nrow = total_obs, ncol = 1, dimnames = NULL)
phenot <- matrix(data = NA, nrow = total_obs, ncol = 1, dimnames = NULL)
counter <- 0

# generate data
for (i in 1:num_fam){
  # generate interaction terms for family across sites
  inter_eff <- rnorm(num_sites,mean_eff_inter,sd_eff_inter)      #generate some random normal deviates
  for (j in 1:num_sites){
    # generates residuals for family at given site
    res_eff <- rnorm(fam_size,mean_residuals,sd_residuals)      #generate some random normal deviates
    for (k in 1:fam_size){
      counter <- counter+1
      self[counter] <- counter
      fam[counter] <- i
      site[counter] <- j
      inter[counter] <- paste(i, j, sep = "_")
      rep[counter] <- k
      phenot[counter] <- global_mean+
        site_eff[j]+
        fam_eff[i]+
        inter_eff[j]+
        res_eff[k]
    }
  }
}

data_toy <- data.frame(cbind(self,fam,site,inter,rep,as.numeric(phenot)))
names(data_toy) <- c("self","fam","site","inter","rep","phenot")
write.table(data_toy, file="data_toy.txt", quote = F, sep = " ", row.names = F)
rm(data_toy)
data_toy <- read.table(file="data_toy.txt",header = T)

# TWO-WAY anova
two_way <- remlf90(fixed = phenot ~ factor(site),
                   random = ~ factor(fam) + factor(inter),
                   data = data_toy)
summary(two_way)
# str(two_way)

# adds blups to data.frame
PBV.full <- ranef(two_way)$`factor(fam)`
PBV <- model.matrix(two_way)$`factor(fam)` %*% PBV.full
data_toy$blup_fam <- PBV
PBV.full <- ranef(two_way)$`factor(inter)`
PBV <- model.matrix(two_way)$`factor(inter)` %*% PBV.full
data_toy$blup_inter <- PBV
data_toy$blup_gi <- data_toy$blup_fam+data_toy$blup_inter

# adds blues to data.frame
fixed_eff <- data.frame(fixef(two_way))
data_toy$ef_site <- fixed_eff[data_toy$site,1]
data_toy$fitted <- fitted(two_way)

# heading="Interactions of family G+I between sites"
heading="family by site NoRs"
x_labels <- as.character(c(1:num_sites))
for (i in 1:num_sites){
  x_labels[i] <- paste("site",i,sep = "")
}
colores <- rainbow(n = num_fam)
linetype <- c(1:num_fam) 

aggdata1 <-aggregate(data_toy$site, by=list(data_toy$fam,data_toy$site),FUN=mean)
aggdata2 <-aggregate(data_toy$fam, by=list(data_toy$fam,data_toy$site),FUN=mean)
aggdata3 <-aggregate(data_toy$fitted, by=list(data_toy$fam,data_toy$site),FUN=mean)
aggdata4 <-aggregate(data_toy$blup_inter, by=list(data_toy$fam,data_toy$site),FUN=mean)
aggdata <- data.frame(cbind(aggdata1$x,aggdata2$x,aggdata3$x,aggdata4$x))
names(aggdata) <- c("site","fam","fitted","inter")

max_val <- max(aggdata$fitted)
min_val <- min(aggdata$fitted)

tiff(filename="fitted_family_effects_by_site.tif",width=500,height=500,compression = c("none"))
plot(aggdata$site, aggdata$fitted, type="n", ylim=c(min_val,max_val), xaxt="n",
  xlab="environments", ylab="fitted family", main=heading, cex=2)
  axis(1, at=(1:num_sites),labels=x_labels, col.axis="black", las=1)
# add lines
for (i in 1:num_fam) {
  fam_set <- subset(aggdata, fam==i)
  lines(fam_set$site, fam_set$fitted, type="b", lwd=1.5,
        lty=linetype[i], col=colores[i])
}
dev.off()

tiff(filename="family_interaction_effects_by_site.tif",width=500,height=500,compression = c("none"))
plot(aggdata$site, aggdata$inter, type="n", ylim=c(min_val,max_val), xaxt="n",
     xlab="environments", ylab="Interaction", main=heading, cex=2)
axis(1, at=(1:num_sites),labels=x_labels, col.axis="black", las=1)
# add lines
for (i in 1:num_fam) {
  fam_set <- subset(aggdata, fam==i)
  lines(fam_set$site, fam_set$inter, type="b", lwd=1.5,
        lty=linetype[i], col=colores[i])
}
dev.off()

# calculates ecovalences from interaction blups
aggdata$sq_inter <- aggdata$inter^2
aggdata_inter <- data.frame(aggregate(aggdata$sq_inter, by=list(aggdata$fam),FUN=sum))
names(aggdata_inter) <- c("fam","sq_inter")
total_sq_inter <- sum(aggdata_inter$sq_inter)
aggdata_inter$ecoval <- aggdata_inter$sq_inter/total_sq_inter

# rank in ascending order families by their ecovalence
rank <- aggdata_inter$fam[order(aggdata_inter$ecoval)]
upper_rank <- 4
upper_threshold <- aggdata_inter$ecoval[rank[upper_rank]]
lower_rank <- num_fam-upper_rank-1
lower_threshold <- aggdata_inter$ecoval[rank[lower_rank]]

heading <- "ecovalences per family"
tiff(filename="ecovalences_by_family.tif",width=500,height=500,compression = c("none"))
plot(aggdata_inter$ecoval, xlab="family", ylab="ecovalence", type="p", main=heading, cex=1)
# add lines
  abline(h = upper_threshold, col="green")
  abline(h = lower_threshold, col="red")
dev.off()

heading="family by site NoRs for highest/lowest ecovalences"
x_labels <- as.character(c(1:num_sites))
for (i in 1:num_sites){
  x_labels[i] <- paste("site",i,sep = "")
}

tiff(filename="extreme_fitted_family_effects_by_site.tif",width=500,height=500,compression = c("none"))
plot(aggdata$site, aggdata$fitted, type="n", ylim=c(min_val,max_val), xaxt="n",
  xlab="environments", ylab="fitted family", main=heading, cex=2)
  axis(1, at=(1:num_sites),labels=x_labels, col.axis="black", las=1)
  # add lines
  for (i in 1:upper_rank) {
    fam_set <- subset(aggdata, fam==rank[i])
    lines(fam_set$site, fam_set$fitted, type="b", lwd=1.5,
          lty=linetype[i], col="green")
  }
  for (i in lower_rank:num_fam) {
    fam_set <- subset(aggdata, fam==rank[i])
    lines(fam_set$site, fam_set$fitted, type="b", lwd=1.5,
          lty=linetype[i], col="red")
  }
dev.off()
  




