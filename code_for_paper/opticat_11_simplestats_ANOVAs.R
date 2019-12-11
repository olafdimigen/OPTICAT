# olaf.dimigen@hu-berlin.de
# 10.06.2018

# This scripts provides very simple (mixed ANOVA) statistical analyses of OPTICAT results
# using .txt/ASCII output exported from MATLAB

# set work directory
setwd("Z:/OPTICA/statistic")
getwd()

library(ez)
library(lme4) # LME

rm(list=ls())

# import data
all_results   = read.csv("Statistics_20181004_asInPaper.txt", header=TRUE)

all_results$Exp     <- factor(all_results$Exp, levels = c(1:2), labels = c("Scenes", "Reading") ) 
all_results$Subject <- factor(all_results$Subject, levels = c(1:24)) 
all_results$HC      <- factor(all_results$HC, levels = c(1:2), labels = c("40 Hz", "100 Hz") ) 
all_results$OW      <- factor(all_results$OW, levels = c(1:2), labels = c("basic", "overweighted") )
all_results$LC_fac  <- factor(all_results$LC, levels = c(1:20), labels = c("0.02 Hz", "0.1 Hz", "0.25 Hz", "0.5 Hz", "0.75 Hz", "1 Hz", "1.5 Hz", "2 Hz", "2.5 Hz", "3 Hz", "3.5 Hz", "4 Hz", "5 Hz" ,"7.5 Hz", "10 Hz", "12.5 Hz", "15 Hz", "20 Hz", "25 Hz", "30 Hz") ) 


############### CR ################
stats_CR <- ezANOVA(all_results, 
                    dv       = DV_CR, 
                    wid      = Subject, 
                    within   = .(HC,LC_fac,OW), 
                    between  = .(Exp),
                    observed = .(DV_CR),
                    detailed = FALSE)
stats_CR


############### SP ################
stats_SP <- ezANOVA(all_results, 
                      dv       = DV_SP, 
                      wid      = Subject, 
                      within   = .(HC,LC_fac,OW), 
                      between  = .(Exp),
                      observed = .(DV_SP),
                      detailed = FALSE)
stats_SP

############### STIM ################
stats_STIM <- ezANOVA(all_results, 
                    dv       = DV_STIM, 
                    wid      = Subject, 
                    within   = .(HC,LC_fac,OW), 
                    between  = .(Exp),
                    observed = .(DV_STIM),
                    detailed = FALSE)
stats_STIM


################################################################
# post-hoc comparisons via t-tests (Bonferoni-corrected)
################################################################

# split dataset into Scenes & Reading
sceneLines <- which(all_results$Exp=="Scenes")
all_results_scenes <- all_results[sceneLines,]
all_results_read   <- all_results[-sceneLines,]


# CR
contrasts_scenes_CR = pairwise.t.test(all_results_scenes$DV_CR,all_results_scenes$LC, p.adj="bonferroni", paired=TRUE)
contrasts_read_CR   = pairwise.t.test(all_results_read$DV_CR,all_results_read$LC, p.adj="bonferroni", paired=TRUE)

# SP
contrasts_scenes_SP = pairwise.t.test(all_results_scenes$DV_SP,all_results_scenes$LC, p.adj="bonferroni", paired=TRUE)
contrasts_read_SP   = pairwise.t.test(all_results_read$DV_SP,all_results_read$LC, p.adj="bonferroni", paired=TRUE)

# STIM
contrasts_scenes_STIM = pairwise.t.test(all_results_scenes$DV_STIM,all_results_scenes$LC, p.adj="bonferroni", paired=TRUE)
contrasts_read_STIM   = pairwise.t.test(all_results_read$DV_STIM,all_results_read$LC, p.adj="bonferroni", paired=TRUE)

contrasts_scenes_CR
contrasts_read_CR
contrasts_scenes_SP
contrasts_read_SP
contrasts_scenes_STIM
contrasts_read_STIM