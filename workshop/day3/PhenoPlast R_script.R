### Case Study on temporal plasticity of height increment in Scots pine 
### (INIA)
###
### R Chambel, A Denardou, E Notivol, J Climent


library(breedR)
library(ggplot2)

# 1 Data

phenplas <- read.csv("F21RIA T4F.csv", sep=";")
View (phenplas)

# 1.1 Some plots of the data (total height)

plot(phenplas$H03) 
hist(phenplas$H03) 
qqnorm(phenplas$H03) 
boxplot(phenplas$H03) 

# 1.2 Some plots of the data (plasticity index) 

plot(phenplas$VPI2) 
hist(phenplas$VPI2) 
qqnorm(phenplas$VPI2) 
boxplot(phenplas$VPI2) 

# 2 Model without blocking

PPres <- remlf90(fixed = VPI2 ~ 1,
                 random = ~ FAM,
                 data = phenplas, method = 'em')
summary (PPres)

# 2.1 Plotting residuals

qplot(fitted(PPres), VPI2,
      data = phenplas) +
  geom_abline(int = 0,
              slope = 1,
              col = 'gray')

# 2.2 Heritability of plasticity (I)

with(PPres, 4*var["FAM",1] / sum(var))

# 3 Model with blocking and interaction

phenplasB <- transform (phenplas,
                        FAM_REP = factor(FAM:factor(REP)))
PPBres <- remlf90(fixed = VPI2 ~ 1+ factor(REP),
                  random = ~ FAM + FAM_REP,
                  data = phenplasB, method = 'em')
summary (PPBres)

# 3.1 Plotting residuals

qplot(fitted(PPBres), VPI2,
      data = phenplasB) +
  geom_abline(int = 0,
              slope = 1,
              col = 'gray')

# 3.2 Family BLUPS for plasticity

ranef(PPBres)$FAM

# 3.2 Heritability of plasticity (II)

with(PPBres,4*var["FAM",1] /sum(var))

# 4 Total height. Model with blocking and interaction 

HBres <- remlf90(fixed = H03 ~ 1+ factor(REP),
                 random = ~ FAM + FAM_REP,
                 data = phenplasB, method = 'em')
summary (HBres)

# 4.1 Plotting residuals

qplot(fitted(HBres), H03,
      data = phenplasB) +
  geom_abline(int = 0,
              slope = 1,
              col = 'gray')

# 4.2 Heritability for total height

with(HBres,4*var["FAM",1] /sum(var))

# 5 Correlation among family BLUPS for total height and plasticity

PIF <- ranef(PPBres)$FAM

HF <- ranef(HBres)$FAM

qplot (HF, PIF) 

cor (HF, PIF) 

