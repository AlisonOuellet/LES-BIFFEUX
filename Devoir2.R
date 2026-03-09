############################
##### PHYLOGÉNIE
##### BIF-4002 H26
############################

# ------ Installation des libraries nécessaires ------ #

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DECIPHER")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")

install.packages('ape')

install.packages("phangorn", dependencies = TRUE)


library(DECIPHER)
library(Biostrings)
library(ape)
library(phangorn)


########################################################
# ----- a) Importation et alignement (DECIPHER) ------ #
########################################################

# Importation des séquences
seqs <- readDNAStringSet("Cetacea COI.fasta")

# Alignement avec paramètres par défaut
alignment_default <- AlignSeqs(seqs)

# Alignement avec pénalités d’indels = 0
alignment_no_gap_penalty <- AlignSeqs(seqs, gapOpening = 0, gapExtension = 0)

# Export pour vérification
writeXStringSet(alignment_default, "alignment_default.fasta")
writeXStringSet(alignment_no_gap_penalty, "alignment_no_gap_penalty.fasta")

##################################################
# ----- b) Traduction et cadre de lecture ------ #
##################################################

# Traduction des séquences nucléotidiques
aa_seqs <- translate(seqs, genetic.code = getGeneticCode("1"))  # code génétique mitochondrial de vertébrés

# Alignement des protéines
aa_alignment <- AlignSeqs(aa_seqs)
writeXStringSet(aa_alignment, "aa_alignment.fasta")

# Décalage d’un cadre de lecture +1
aa_seqs_shift <- translate(subseq(seqs, start=2), genetic.code = getGeneticCode("1"))
aa_alignment_shift <- AlignSeqs(aa_seqs_shift)


#############################################
# ------ c) Phylogénie nucléotidique ------ #
#############################################

# Alignement en format "DNAbin" pour ape
dna_bin <- as.DNAbin(alignment_default)

# Matrices de distance avec les bons codes
dist_JC   <- dist.dna(dna_bin, model = "JC69")
dist_K2P  <- dist.dna(dna_bin, model = "K80")
dist_TN93 <- dist.dna(dna_bin, model = "TN93")
dist_GG   <- dist.dna(dna_bin, model = "GG95") # Correction ici

# Force la conversion en matrice.
dna_mat <- as.matrix(dna_bin)

# Vérifie que c'est bien une matrice
is.matrix(dna_mat) # Doit retourner TRUE

# Analyse Bootstrap
run_bootstrap <- function(dna_input, mod) {
    # Arbre de référence
    d <- dist.dna(dna_input, model = mod)
    tree <- nj(d)

    # Bootstrap (1000 itérations)
    boot_scores <- boot.phylo(tree, dna_input, function(x) nj(dist.dna(x, model = mod)), B = 1000)

    return(boot_scores)
}

# Lance les calculs sur dna_mat
boot_JC   <- run_bootstrap(dna_mat, "JC69")
boot_K2P  <- run_bootstrap(dna_mat, "K80")
boot_TN93 <- run_bootstrap(dna_mat, "TN93")
boot_GG   <- run_bootstrap(dna_mat, "GG95")

# Calcul des moyennes en ignorant les NA
scores_moyens <- c(
    JC   = mean(boot_JC, na.rm = TRUE),
    K2P  = mean(boot_K2P, na.rm = TRUE),
    TN93 = mean(boot_TN93, na.rm = TRUE),
    GG   = mean(boot_GG, na.rm = TRUE)
)

print(scores_moyens)

# Identifie le modèle le plus robuste
meilleur_ape <- names(which.max(scores_moyens))
cat("Le modèle le plus robuste selon le bootstrap est :", meilleur_ape, "\n")


# Conversion au format phyDat
dna_pd <- as.phyDat(dna_mat)

# Test de modèles
mt <- modelTest(dna_pd)

# Trouve le meilleur modèle selon le critère BIC
best_model <- mt[which.min(mt$BIC), ]
print(best_model)

# Transforme le meilleur résultat du test en objet PML
fit_best <- as.pml(mt, model = "BIC")

alpha_val <- fit_best$shape
inv_val   <- fit_best$inv

cat("Paramètre Gamma (alpha) :", alpha_val, "\n")
cat("Sites invariants (I) :", inv_val, "\n")

# Optimisation de l'arbre
fit_optim <- optim.pml(fit_best, optGamma = TRUE, optInv = TRUE, optEdge = TRUE)

# Affichage de l'arbre final
par(mar = c(1, 1, 3, 1))
plot(fit_optim$tree, main = paste("Arbre Final (ML) - Modèle :", best_model$Model), cex = 0.8)
add.scale.bar()
