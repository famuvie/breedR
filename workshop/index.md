---
title: "Joint training workshop on breedR and phenotypic plasticity"
output: pdf_document
permalink: /workshop/
layout: page
---

[![T4F](../images/Logo-t4f.png)](http://www.trees4future.eu/) ![EU](../images/logo_eu.png)
*This workshop is organised with funds from the European Union’s Seventh Framework Programme for research, technological development and demonstration under grant agreement n° 284181 (“Trees4Future”)*

June 30 – July 3 2015

Jaca, Spain


#### Preliminary tasks

Participants are expected to bring their laptops and have the following
tasks done, if necessary, before the workshop

- Install and get familiar with [R](../getR)
- Install the latest version of [breedR](http://famuvie.github.io/breedR/)
- Join our discussion [group](http://groups.google.com/group/breedr)
- Bring your **datasets** and prepare your **questions**

### Provisional program
last updated March 18, 2015


#### Day 1: Tuesday June 30

**Introduction to Linear Mixed Models (LMMs) and breedR**

Participants are expected to arrive in the course of the morning.

**09.00** Workshop opening ([slides](day1/intro_WP6.pdf))

**09.30** Morning Workshop (1.5 hs.) Facundo Muñoz (INRA Orléans)

- [slides](day1/Overview.html) · [code](day1/Overview.R)
- breedR overview: 
  - Case studies in common field trials
  - The additive genetic model
  - breedR models for spatial autocorrelation and competition

**11.00** Theoretical session (1.5 hs.) Luis Varona (U. de Zaragoza)

- Statistical inference (frequentist/Bayesian; (RE)ML/MCMC)
- Linear Mixed Models

**Genetic models, spatial autocorrelation and competition**

**14.00** Theoretical session (1.5 hs.) Eduardo Cappa (INTA Buenos Aires)
- [slides](day1/ECappa_theory.pdf)
- Introduction to environment heterogeneity and competition effects
- Diagnosis of spatial and competition effects
- Autoregressive and splines spatial models
- Competition model
- Examples of application


**16.00** Afternoon Workshop (3 hs.) Facundo Muñoz (INRA Orléans); Eduardo Cappa (INTA Buenos Aires)
- scripts: [Douglas-fir](day1/Douglas-fir_script.pdf); [E. globulus](day1/E.globulus_script.pdf)
- Diagnosis of spatial and competition effects
- Analysis of spatial autocorrelation
- Analysis of Competition


#### Day 2: Wednesday July 1st


**09.00 Field excursion** visit natural populations of mixed Beech-Fir and Romanesque
monastery *San Juan de la Peña*

**Genotype by Environment and multi-trait models**

**14.00** Theoretical session (1.5 hs.) Leopoldo Sánchez (INRA Orleans)
- Genetic correlations ([slides](day2/WED genetic_correlations concepts LSanchez.pdf))
- G×E analysis ([slides](day2/WED gbye concepts LSanchez.pdf))

**16.00** GxE Workshop (2 hs.)
- Toy examples of interaction analysis ([slides](day2/ecovalence_generator.html); [script](day2/ecovalence generator.R))
- Case study in G×E interaction in a multi-site context ([slides](day2/GEI.html); [script](day2/GEI.R))

**18.00** Afternoon workshop (1.5 hs.)
- datasets: (genotypes [case-1](day2/genotypes.txt); [case-2](day2/genotypes_case2.txt); phenotypes [case-1](day2/pheno_ped.txt); [case-2](day2/pheno_ped_case2.txt))
- Simulated examples with Genomic data ([slides-intro](day2/WED Brief_intro_GS LSanchez.pdf); [excercise](day2/gs_case.html); [script](day2/metagene_gs.R))


#### Day 3: Thursday July 2

**WP3 thematic network on Phenotypic Plasticity**

**09.00** Theoretical session (1.5 hs.) José Climent (CIFOR-INIA)
- Phenotypic plasticity ([slides](day3/BreedR_Jaca_Plasticity_Climent.pdf))
- Analysis of phenotipic plasticity ([dataset](day3/F21RIA T4F.csv); [script](day3/PhenoPlast R_script.R))

**11.00** Theoretical session (1.5 hs.)
- Longitudinal models ([slides](day3/RR model Jaca.pdf))
- Larix longitudinal example ([slides](day3/Longitudinal.html); [script](day3/Longitudinal.R))

**14.00** Afternoon workshop (5 hs.)
- Discussion and analysis of participant's use cases
  - Mesfin case study ([slides](day3/Mesfin_spatial_analysis.pdf); [script](day3/Mesfin_spatial_analysis.R); [dataset](day3/DW.predicted.practice.txt))


#### Day 4: Friday July 3

**09.00** Participant's session
- Open time for excercises

**12.30** Wrap up and conclusons


#### Further possible contents

- Genfored as a data base (E. Notivol)
- Possible extensions of breedR (A. Legarra)


