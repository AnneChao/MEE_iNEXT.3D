# ========================================================================================================== #

library(ape)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(gg.gap)


# ========================================================================================================== #

library(devtools)
install_github("AnneChao/iNEXT.3D")  # Press 'Enter' to skip number selection.
library(iNEXT.3D)

source("Source R code.txt")

# ========================================================================================================== #
# Part 1 : Abundance-based yearly data analysis (plotting Figure 1 and 2 in the MEE paper).
# See "Brief guide" for details. 

Abun <- read.csv("Fish abundance data.csv", row.names = 1, header= TRUE)
tree <- read.tree("Fish phyloTree.txt")
traits <- read.csv("Fish traits.csv", row.names = 1, header= TRUE)


Cmax <- apply(Abun, 2, function(x) iNEXT.3D:::Coverage(x, 'abundance', 2*sum(x))) %>% min %>% round(., 4)
Cmin <- apply(Abun, 2, function(x) iNEXT.3D:::Coverage(x, 'abundance', sum(x))) %>% min %>% round(., 4)


# ========================================================================================================== #
# Figure 1 - Taxonomic diversity

TD_est <- estimate3D(data = Abun, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance', base = 'coverage', 
                     level = c(Cmin, Cmax), nboot = 0)
TD_obs <- AO3D(data = Abun, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance', nboot = 0, method = 'Observed')
TD_asy <- AO3D(data = Abun, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance', nboot = 0, method = 'Asymptotic')


out_TD <- rbind(TD_est %>% select(Assemblage, Order.q, qTD, SC), 
                TD_obs %>% select(Assemblage, Order.q, qTD) %>% mutate(SC = 'Observed') , 
                TD_asy %>% select(Assemblage, Order.q, qTD) %>% mutate(SC = 'Asymptotic'))
fig_1_or_3(out_TD, y_label = 'Taxonomic diversity')


# ========================================================================================================== #
# Figure 1 - Phylogenetic diversity

PD_est <- estimate3D(data = Abun, diversity = 'PD', q = c(0, 1, 2), datatype = 'abundance', base = 'coverage',
                     level = c(Cmin, Cmax), nboot = 0, PDtree = tree, PDreftime = 1)
PD_obs <- AO3D(data = Abun, diversity = 'PD', q = c(0, 1, 2), datatype = 'abundance', 
               nboot = 0, PDtree = tree, PDreftime = 1, method = 'Observed') 
PD_asy <- AO3D(data = Abun, diversity = 'PD', q = c(0, 1, 2), datatype = 'abundance', 
               nboot = 0, PDtree = tree, PDreftime = 1, method = 'Asymptotic')


out_PD <- rbind(PD_est %>% select(Assemblage, Order.q, qPD, SC), 
                PD_obs %>% select(Assemblage, Order.q, qPD) %>% mutate(SC = 'Observed'), 
                PD_asy %>% select(Assemblage, Order.q, qPD) %>% mutate(SC = 'Asymptotic'))
fig_1_or_3(out_PD, y_label = 'Phylogenetic diversity')


# ========================================================================================================== #
# Figure 1 - Functional diversity

for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[,i] <- factor(traits[,i], levels = unique(traits[,i]))
}
distM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()


FD_est <- estimate3D(data = Abun, diversity = 'FD', q = c(0, 1, 2), datatype = 'abundance', base = 'coverage',
                     level = c(Cmin, Cmax), nboot = 0, FDdistM = distM)
FD_obs <- AO3D(data = Abun, diversity = 'FD', q = c(0, 1, 2), datatype = 'abundance', nboot = 0, FDdistM = distM, method = 'Observed')
FD_asy <- AO3D(data = Abun, diversity = 'FD', q = c(0, 1, 2), datatype = 'abundance', nboot = 0, FDdistM = distM, method = 'Asymptotic')


out_FD <- rbind(FD_est %>% select(Assemblage, Order.q, qFD, SC), 
                FD_obs %>% select(Assemblage, Order.q, qFD) %>% mutate(SC = 'Observed') , 
                FD_asy %>% select(Assemblage, Order.q, qFD) %>% mutate(SC = 'Asymptotic'))
fig_1_or_3(out_FD, y_label = 'Functional diversity')


# ========================================================================================================== #
# Figure 2

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 0)

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 1)

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 2)


# ========================================================================================================== #
# Part 2 : Three-year incidence data analysis for Figures 3 and 4.

Inci_raw <- read.csv("Fish incidence raw data.csv", row.names = 1, header= TRUE)
nT <- read.csv('nT for incidence data.csv', row.names = 1)
tree <- read.tree("Fish phyloTree.txt")
traits <- read.csv("Fish traits.csv", row.names = 1, header= TRUE)


Cmax <- sapply(1:length(nT), function(i) rowSums( Inci_raw[, (sum(nT[1:i]) - sum(nT[i]) + 1) : sum(nT[1:i])] )) %>% rbind(as.integer(nT),.) %>% 
  apply(., 2, function(x) iNEXT.3D:::Coverage(x, 'incidence_freq', 2*x[1])) %>% min %>% round(., 4)
Cmin <- sapply(1:length(nT), function(i) rowSums( Inci_raw[, (sum(nT[1:i]) - sum(nT[i]) + 1) : sum(nT[1:i])] )) %>% rbind(as.integer(nT),.) %>% 
  apply(., 2, function(x) iNEXT.3D:::Coverage(x, 'incidence_freq', x[1])) %>% min %>% round(., 4)


# ========================================================================================================== #
# Figure 3 - Taxonomic diversity

TD_est <- estimate3D(data = Inci_raw, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_raw', base = 'coverage',
                     level = c(Cmin, Cmax), nboot = 0, nT = nT)
TD_obs <- AO3D(data = Inci_raw, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_raw', nboot = 0, nT = nT, method = 'Observed')
TD_asy <- AO3D(data = Inci_raw, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_raw', nboot = 0, nT = nT, method = 'Asymptotic')


out_TD <- rbind(TD_est %>% select(Assemblage, Order.q, qTD, SC), 
                TD_obs %>% select(Assemblage, Order.q, qTD) %>% mutate(SC = 'Observed'), 
                TD_asy %>% select(Assemblage, Order.q, qTD) %>% mutate(SC = 'Asymptotic'))
fig_1_or_3(out_TD, y_label = 'Taxonomic diversity')


# ========================================================================================================== #
# Figure 3 - Phylogenetic diversity

PD_est <- estimate3D(data = Inci_raw, diversity = 'PD', q = c(0, 1, 2), datatype = 'incidence_raw', base = 'coverage',
                     level = c(Cmin, Cmax), nboot = 0, nT = nT, PDtree = tree, PDreftime = 1)
PD_obs <- AO3D(data = Inci_raw, diversity = 'PD', q = c(0, 1, 2), datatype = 'incidence_raw',
               nboot = 0, nT = nT, PDtree = tree, PDreftime = 1, method = 'Observed')
PD_asy <- AO3D(data = Inci_raw, diversity = 'PD', q = c(0, 1, 2), datatype = 'incidence_raw',
               nboot = 0, nT = nT, PDtree = tree, PDreftime = 1, method = 'Asymptotic')


out_PD <- rbind(PD_est %>% select(Assemblage, Order.q, qPD, SC), 
                PD_obs %>% select(Assemblage, Order.q, qPD ) %>% mutate(SC = 'Observed'), 
                PD_asy %>% select(Assemblage, Order.q, qPD) %>% mutate(SC = 'Asymptotic'))
fig_1_or_3(out_PD, y_label = 'Phylogenetic diversity')


# ========================================================================================================== #
# Figure 3 - Functional diversity

for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[, i] <- factor(traits[,i], levels = unique(traits[, i]))
}
distM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()


FD_est <- estimate3D(data = Inci_raw, diversity = 'FD', q = c(0, 1, 2), datatype = 'incidence_raw', base = 'coverage',
                     level = c(Cmin, Cmax), nboot = 0, nT = nT, FDdistM = distM)
FD_obs <- AO3D(data = Inci_raw, diversity = 'FD', q = c(0, 1, 2), datatype = 'incidence_raw', 
               nboot = 0, nT = nT, FDdistM = distM, method = 'Observed')
FD_asy <- AO3D(data = Inci_raw, diversity = 'FD', q = c(0, 1, 2), datatype = 'incidence_raw',
               nboot = 0, nT = nT, FDdistM = distM, method = 'Asymptotic')

out_FD <- rbind(FD_est %>% select(Assemblage, Order.q, qFD, SC), 
                FD_obs %>% select(Assemblage, Order.q, qFD) %>% mutate(SC = 'Observed'), 
                FD_asy %>% select(Assemblage, Order.q, qFD) %>% mutate(SC = 'Asymptotic'))
fig_1_or_3(out_FD, y_label = 'Functional diversity')


# ========================================================================================================== #
# Figure 4

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 0)

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 1)

fig_2_or_4(TD.output = out_TD, PD.output = out_PD, FD.output = out_FD, q = 2)


# ========================================================================================================== #
# Alternatively, for Figure 3a (TD) and Fifure 3c (FD), one can also use incidence-frequency data, but not 
# for PD because the incidence frequency of any interior node of a phylogenetic tree depends on raw data. 

Inci_freq <- read.csv("Fish incidence frequency data.csv", row.names = 1, header= TRUE)
TD_est <- estimate3D(data = Inci_freq, diversity = 'TD', q = c(0, 1, 2), datatype = 'incidence_freq', 
                     base = 'coverage', level = c(Cmin, Cmax), nboot = 0)
