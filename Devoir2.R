############################
##### PHYLOGÉNIE
##### BIF-4002 H26
############################

# ------ Installation des libraries nécessaires ------ #

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DECIPHER")

install.packages('ape')

install.packages("phangorn", dependencies = TRUE)


library(DECIPHER)
library(ape)
library(phangorn)

# ----- a) Importation et alignement (DECIPHER) ------ #

