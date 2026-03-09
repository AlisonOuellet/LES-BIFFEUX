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
# Importation des données
data <- readDNAStringSet("Cetacea COI.fasta")

# Alignement
aligned_data <- AlignSeqs(data)

# Retire les gaps de l'alignement pour retrouver les séquences brutes
sequences_brutes <- RemoveGaps(aligned_data)

# Alignement avec les pénalités à zéro
aligned_no_gap <- AlignSeqs(sequences_brutes, gapOpening = 0, gapExtension = 0)


# Visualiser l'alignement dans un browser (optionnel)
BrowseSeqs(aligned_no_gap, highlight = 0)

# ----- b) Traduction et cadre de lecture ------ #

