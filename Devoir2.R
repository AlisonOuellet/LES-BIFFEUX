############################
##### PHYLOGÉNIE
##### BIF-4002 H26
############################

# ------ Installation des libraries nécessaires ------ #

if (!requireNamespace("DECIPHER", quietly = TRUE)) BiocManager::install("DECIPHER")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn")

library(DECIPHER)
library(Biostrings)
library(ape)
library(phangorn)


########################################################
# ----- a) Importation et alignement (DECIPHER) ------ #
########################################################

clean_labels_standard <- function(labels) {

    # Extraire "Genre espèce"
    species <- sub(".*?([A-Z][a-z]+\\s[a-z]+).*", "\\1", labels)

    # Extraire UAM:Mamm même sans parenthèses
    uam <- ifelse(grepl("UAM:Mamm:[0-9]+", labels),
                  sub(".*(UAM:Mamm:[0-9]+).*", "\\1", labels),
                  "")

    # Construire label final
    final_labels <- ifelse(uam != "",
                           paste0(species, " (", uam, ")"),
                           species)

    return(final_labels)
}

# Importation des séquences
seqs <- readDNAStringSet("Cetacea COI.fasta", format = "fasta")
names(seqs) <- clean_labels_standard(names(seqs))

# Alignement avec paramètres par défaut
alignment_default <- AlignSeqs(seqs)

# Alignement avec pénalités d’indels = 0
alignment_no_gap <- AlignSeqs(seqs, gapOpening = 0, gapExtension = 0)

# Comparer les alignements
width(alignment_default)
width(alignment_no_gap)

# Nombre total de gaps
sum(letterFrequency(alignment_default, "-"))
sum(letterFrequency(alignment_no_gap, "-"))

# Export pour vérification
writeXStringSet(alignment_default, "alignment_default.fasta")
writeXStringSet(alignment_no_gap, "alignment_no_gap.fasta")

##################################################
# ----- b) Traduction et cadre de lecture ------ #
##################################################

# Alignement d'acides aminées
aa_alignment_NA <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"), type = "AAStringSet", readingFrame = NA)
writeXStringSet(aa_alignment_NA, "aa_alignment_NA.fasta")

aa_alignment_1 <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"), type = "AAStringSet", readingFrame = 1)
writeXStringSet(aa_alignment_1, "aa_alignment_1.fasta")

# Vérification
# Nombre de codons stop par séquence
stops_NA <- letterFrequency(aa_alignment_NA, "*")
stops_frame1 <- letterFrequency(aa_alignment_1, "*")

# Résumé statistique (min, médiane, max)
summary(stops_NA)
summary(stops_frame1)

# Distribution : combien de séquences ont 0, 1, 2 stops, etc.
table(stops_NA)
table(stops_frame1)

# Identifier les séquences contenant au moins un codon stop
names(aa_alignment_NA)[stops_NA > 0]
names(aa_alignment_1)[stops_frame1 > 0]


#############################################
# ------ c) Phylogénie nucléotidique ------ #
#############################################

# Conversion de l'alignement en format DNAbin (ape)
dna <- as.DNAbin(alignment_default)

# Matrices de distances avec les différents modèles
dist_JC  <- dist.dna(dna, pairwise.deletion = FALSE, model = "JC69")    # Jukes-Cantor
dist_K80 <- dist.dna(dna, pairwise.deletion = FALSE, model = "K80")     # Kimura 2 paramètres
dist_TN  <- dist.dna(dna, pairwise.deletion = FALSE, model = "TN93")    # Tamura-Nei
dist_GG95 <- dist.dna(dna, pairwise.deletion = FALSE, model = "GG95")  # Galtier-Gouy (LogDet)

# Construction des arbres
arbre_JC  <- nj(dist_JC)
arbre_K80 <- nj(dist_K80)
arbre_TN  <- nj(dist_TN)
arbre_GG95 <- nj(dist_GG95)

# Visualisation
par(mar = c(1, 1, 1, 1))
plot(arbre_JC, cex = 0.9, main="Jukes et Cantor", edge.width = 1.5, x.lim = c(0, max(node.depth.edgelength(arbre_JC)) * 2.5))
plot(arbre_K80, cex = 0.9, main="Kimura 2 parametres 1980", edge.width = 1.5, x.lim = c(0, max(node.depth.edgelength(arbre_JC)) * 2.5))
plot(arbre_TN, cex = 0.9, main="Tamura et Nei 1993", edge.width = 1.5, x.lim = c(0, max(node.depth.edgelength(arbre_JC)) * 2.5))
plot(arbre_GG95, cex = 0.9, main="Galtier and Gouy 1995", edge.width = 1.5)

#Comparaison des modèles entre eux par corrélation (Pearson)
round(cor(cbind(dist_JC, dist_TN,dist_K80, dist_GG95)),3)

# boot.phylo nécessite une matrice
dna_matrix <- as.matrix(dna)

# Bootstrap (1000 itérations)
boot_JC <- boot.phylo(phy = arbre_JC, x = as.matrix(dna), FUN = function(xx) nj(dist.dna(xx,pairwise.deletion = FALSE, model = "JC69")), B = 1000)
boot_K80 <- boot.phylo(phy = arbre_K80, x = as.matrix(dna), FUN = function(xx) nj(dist.dna(xx,pairwise.deletion = FALSE, model = "K80")), B = 1000)
boot_TN <- boot.phylo(phy = arbre_TN, x = as.matrix(dna), FUN = function(xx) nj(dist.dna(xx,pairwise.deletion = FALSE, model = "TN93")), B = 1000)
boot_GG95 <- boot.phylo(phy = arbre_GG95, x = as.matrix(dna), FUN = function(xx) nj(dist.dna(xx,pairwise.deletion = FALSE, model = "GG95")), B = 1000)

summary(boot_JC)
summary(boot_K80)
summary(boot_TN)
summary(boot_GG95)

plot(arbre_JC, main = "Jukes et Cantor", cex = 0.9, edge.width = 1.5, label.offset = 0.002)
nodelabels(boot_JC/10, frame = "circle", bg = "#FFFFFFCC", cex = 0.5, adj = c(0.5), col = "red")

plot(arbre_K80, main="Kimura 2 paramètres", cex = 0.9, edge.width = 1.5, label.offset = 0.002)
nodelabels(boot_K80/10, frame="circle", bg = "#FFFFFFCC", cex = 0.5, adj = c(0.5), col = "red")

plot(arbre_TN, main="Tamura-Nei", cex = 0.9, edge.width = 1.5, label.offset = 0.002)
nodelabels(boot_TN/10, frame="circle", bg = "#FFFFFFCC", cex = 0.5, adj = c(0.5), col = "red")

plot(arbre_GG95, main="Galtier-Gouy", cex = 0.9, edge.width = 1.5, label.offset = 0.015, x.lim = 1.5)
nodelabels(boot_GG95/10, frame="circle", bg = "#FFFFFFCC", cex = 0.5, adj = c(0.5), col = "red")

round(cor(cbind(boot_JC, boot_TN, boot_K80, boot_GG95),
          use = "pairwise.complete.obs"),3)

# Conversion vers format phangorn
dna_phy <- as.phyDat(dna)

# Test de modèles
model_test <- phangorn::modelTest(dna_phy)
model_test

# tester les modèles évolutifs
env <- attr(model_test, "env")
ls(env = env)

best_model <- model_test[which.min(model_test$AIC), ]
best_model

dm <- dist.ml(dna_phy)
treeNJ <- NJ(dm)

plot(treeNJ, main="Neighbor-Joining")

fit <- pml(treeNJ, data = dna_phy)
fitGTR <- update(fit, k = 4, inv = 0.2)

fitGTR <- optim.pml(
    fitGTR,
    model = "GTR",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    control = pml.control(trace = 0)
)

bs <- bootstrap.pml(
    fitGTR,
    bs = 1000,
    optNni = TRUE,
    control = pml.control(trace = 0)
)

par(mar = c(1, 1, 1, 1))
plotBS(
    midpoint(fitGTR$tree),
    bs,
    p = 0,
    type = "p",
    frame = "circle",
    cex = 0.7,
    bs.adj = c(0.7, 0.7),
    bg = "#FFFFFFCC",
    bs.col = "red"
)



#################################################
# ------ d) Alignement des acides aminés ------ #
#################################################

# Conversion des alignements AA en format phangorn
aa_auto_ape <- as.AAbin(aa_alignment_NA)
aa_auto_phy <- as.phyDat.AAbin(aa_auto_ape, type = "AA")

aa_frame1_ape <- as.AAbin(aa_alignment_1)
aa_frame1_phy <- as.phyDat.AAbin(aa_frame1_ape, type = "AA")

model_test_aa <- phangorn::modelTest(aa_auto_phy)
model_test_aa

# Meilleur modèle
model_test_aa[which.min(model_test_aa$AIC),]

dist_LG <- dist.ml(aa_auto_phy, model="LG")
dist_JTT <- dist.ml(aa_auto_phy, model="JTT")
dist_BLOSUM <- dist.ml(aa_auto_phy, model="Blosum62")

tree_LG <- NJ(dist_LG)
tree_JTT <- NJ(dist_JTT)
tree_BLOSUM <- NJ(dist_BLOSUM)

par(mar = c(1, 1, 1, 1))
plot(tree_LG, main="NJ - LG")
plot(tree_JTT, main="NJ - JTT")
plot(tree_BLOSUM, main="NJ - BLOSUM62")

bs_LG <- bootstrap.phyDat(
    aa_auto_phy,
    FUN=function(x) NJ(dist.ml(x, model="LG")),
    bs=1000
)

bs_JTT <- bootstrap.phyDat(
    aa_auto_phy,
    FUN=function(x) NJ(dist.ml(x, model="JTT")),
    bs=1000
)

bs_BLOSUM <- bootstrap.phyDat(
    aa_auto_phy,
    FUN=function(x) NJ(dist.ml(x, model="Blosum62")),
    bs=1000
)

par(mar = c(1, 1, 1, 1))

plotBS(tree_LG, bs_LG, main="NJ LG + bootstrap", frame="circle",     bg = "#FFFFFFCC",
       bs.col = "red", cex=0.6, p=0)

plotBS(tree_JTT, bs_JTT, main="NJ JTT + bootstrap",frame="circle",     bg = "#FFFFFFCC",
       bs.col = "red", cex=0.6, p=0)

plotBS(tree_BLOSUM, bs_BLOSUM, main="NJ BLOSUM62 + bootstrap",frame="circle",     bg = "#FFFFFFCC",
       bs.col = "red", cex=0.6, p=0)

dist_LG_f1 <- dist.ml(aa_frame1_phy, model="LG")
tree_LG_f1 <- NJ(dist_LG_f1)

bs_LG_f1 <- bootstrap.phyDat(
    aa_frame1_phy,
    FUN=function(x) NJ(dist.ml(x, model="LG")),
    bs=1000
)

par(mar = c(1, 1, 1, 1))
plotBS(tree_LG_f1, bs_LG_f1, main="Frame 1 corrigé",frame="circle", bg = "#FFFFFFCC",
       bs.col = "red", cex=0.6, p=0)
###########################################################
# ------ e) Analyse modelTest par cadre de lecture ------ #
###########################################################

# Alignements AA pour chaque cadre
aa_frame1 <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"),
                              type = "AAStringSet", readingFrame = 1)
aa_frame2 <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"),
                              type = "AAStringSet", readingFrame = 2)

aa_frame3 <- AlignTranslation(seqs, geneticCode = getGeneticCode("2"),
                              type = "AAStringSet", readingFrame = 3)


# Conversion en phyDat
aa1_phy <- as.phyDat(as.AAbin(aa_frame1))
aa2_phy <- as.phyDat(as.AAbin(aa_frame2))
aa3_phy <- as.phyDat(as.AAbin(aa_frame3))
#################################
# ModelTest pour chaque cadre
#################################

mt_aa1 <- modelTest(aa1_phy)
mt_aa2 <- modelTest(aa2_phy)
mt_aa3 <- modelTest(aa3_phy)

# Meilleurs modèles
best_aa1 <- mt_aa1[which.min(mt_aa1$AIC),]
best_aa2 <- mt_aa2[which.min(mt_aa2$AIC),]
best_aa3 <- mt_aa3[which.min(mt_aa3$AIC),]

best_aa1
best_aa2
best_aa3

#################################
# Reconstruction arbres + bootstrap
#################################

# Créer les arbres de base (Neighbor-Joining) pour initialiser pml()
tree1 <- NJ(dist.ml(aa1_phy))
tree2 <- NJ(dist.ml(aa2_phy))
tree3 <- NJ(dist.ml(aa3_phy))

model1 <- gsub("\\+.*", "", best_aa1$Model)
model2 <- gsub("\\+.*", "", best_aa2$Model)
model3 <- gsub("\\+.*", "", best_aa3$Model)

fit1 <- pml(tree1, data = aa1_phy)
fit2 <- pml(tree2, data = aa2_phy)
fit3 <- pml(tree3, data = aa3_phy)

fit1_opt <- optim.pml(fit1, model = model1, optGamma = TRUE)
fit2_opt <- optim.pml(fit2, model = model2, optGamma = TRUE)
fit3_opt <- optim.pml(fit3, model = model3, optGamma = TRUE)

bs1 <- bootstrap.pml(fit1_opt, bs = 1000)
bs2 <- bootstrap.pml(fit2_opt, bs = 1000)
bs3 <- bootstrap.pml(fit3_opt, bs = 1000)

par(mar = c(1,1,1,1))

plotBS(fit1_opt$tree, bs1,
       main="Frame 1",
       frame="circle",
       bg = "#FFFFFFCC",
       bs.col = "red",
       cex=0.6,
       p=0,
       use.edge.length = FALSE)

plotBS(fit1_opt$tree, bs1,
       main="Frame 1",
       frame="circle",
       bg = "#FFFFFFCC",
       bs.col = "red",
       cex=0.6,
       p=0)

plotBS(fit2_opt$tree, bs2,
       main="Frame 2",
       frame="circle",
       bg = "#FFFFFFCC",
       bs.col = "red",
       cex=0.6,
       p=0)
plotBS(fit2_opt$tree, bs2,
       main="Frame 2",
       frame="circle",
       bg = "#FFFFFFCC",
       bs.col = "red",
       cex=0.6,
       p=0,
       use.edge.length = FALSE)

plotBS(fit3_opt$tree, bs3,
       main="Frame 3",
       frame="circle",
       bg = "#FFFFFFCC",
       bs.col = "red",
       cex=0.7,
       p=0)
plotBS(fit3_opt$tree, bs3,
       main="Frame 3",
       frame="circle",
       bg = "#FFFFFFCC",
       bs.col = "red",
       cex=0.7,
       p=0,
       use.edge.length = FALSE)
