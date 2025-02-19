library(phangorn)
library(ggtree)
library(ggplot2)
library(ape)

# Analysis for reference tree
setwd("/Users/corinna/Documents/Work/Schwartz_Lab/Plant_paralog_evolution/Bryophyta/CESA/CESA_trees_refs")

#Tree with de novo assembled contigs for reference species
data <- read.table("tree_metadata.txt", header=T)
col <- c("Anthoceros" = "red2", "Ceratodon" = "cornflowerblue", "Entodon" = "mediumturquoise", "Fontinalis" = "darkcyan", "Hypnum" = "lightseagreen", "Marchantia" = "mediumpurple", "Penium" = "royalblue3", "Physcomitrium" = "plum", "Sphagnum" = "violetred4", "Takakia" = "darkgoldenrod1")

tree <- read.tree("RAxML_bipartitions.CESA_ref_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node)) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=109, resolve.root = TRUE, edgelabel=TRUE)
 
pdf("RAxML_bipartitions.CESA_ref_exon.pdf", width=10, height=10)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + geom_tiplab(align=FALSE, hjust=-.15, fontface=3) + xlim(0, 2.5) + geom_treescale(x=0, y=53) + geom_nodelab(geom="label", color="black", fill="whitesmoke", size=3, label.size=NA, nudge_x=-0.05)
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Genus)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Genus)) + scale_color_manual(values = col, name="Genus", na.value="black")
t2
dev.off()

tree <- read.tree("RAxML_bestTree.CESA_ref_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node)) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=110, resolve.root = TRUE, edgelabel=TRUE)
 
pdf("RAxML_bestTree.CESA_ref_exon.pdf", width=10, height=10)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + geom_tiplab(align=FALSE, hjust=-.15, fontface=3) + xlim(0, 2.5) + geom_treescale(x=0, y=53)
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Genus)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Genus)) + scale_color_manual(values = col, name="Genus", na.value="black")
t2
dev.off()

#Tree with all reference genes (de novo assembled and BLASTx searches)
data <- read.table("tree_metadata2.txt", header=T)
col <- c("Anthoceros" = "red2", "Ceratodon" = "cornflowerblue", "Entodon" = "mediumturquoise", "Fontinalis" = "darkcyan", "Hypnum" = "lightseagreen", "Marchantia" = "mediumpurple", "Penium" = "royalblue3", "Physcomitrium" = "plum", "Sphagnum" = "violetred4", "Takakia" = "darkgoldenrod1")

tree <- read.tree("RAxML_bipartitions.CESA_all_refs_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node)) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=200, resolve.root = TRUE, edgelabel=TRUE)

pdf("RAxML_bipartitions.CESA_all_refs_exon.pdf", width=10, height=15)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + geom_tiplab(align=FALSE, hjust=-.15, fontface=3) + xlim(0, 2.5) + geom_treescale(x=0, y=100) + geom_nodepoint(aes(subset = label > 75), size=3, shape=21, fill="darkgray", color="black")
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Genus)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Genus)) + scale_color_manual(values = col, name="Genus", na.value="black")
t2
dev.off()


# Analysis for full tree
setwd("/Users/corinna/Documents/Work/Schwartz_Lab/Plant_paralog_evolution/Bryophyta/CESA/CESA_trees_all_samples")

data <- read.table("tree_metadata.txt", header=T)
col <- c("Andreaeales" = "royalblue4", "Bartramiales" = "mediumturquoise", "Bryales" = "darkcyan", "Buxbaumiales" = "orange", "Dicranales" = "pink1", "Diphysciales" = "slateblue", "Encalyptales" = "plum4", "Grimmiales" = "palevioletred", "Orthotrichales" = "paleturquoise4", "Polytrichales" = "tomato", "Pottiales" = "hotpink", "Rhizogoniales" = "paleturquoise3", "Scouleriales" = "deeppink3", "Tetraphidales" = "gold", "Timmiales" = "salmon", "Anthocerotales" = "red2", "Pseudoditrichales" = "cornflowerblue", "Hypnales" = "lightseagreen", "Marchantiales" = "mediumpurple", "Desmidiales" = "royalblue3", "Funariales" = "plum", "Sphagnales" = "violetred4", "Takakiales" = "darkgoldenrod1")

tree <- read.tree("RAxML_bipartitions.CESA_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node), size=3) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=662, resolve.root = TRUE, edgelabel=TRUE)
 
pdf("RAxML_bipartitions.CESA_exon.pdf", width=15, height=60)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + geom_tiplab(align=FALSE, hjust=-.15, fontface=3) + xlim(0, 2.5) + geom_treescale(x=0, y=480) + geom_nodelab(geom="label", color="black", fill="whitesmoke", size=3, label.size=NA, nudge_x=-0.05)
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black")
t2
dev.off()

t <- ggtree(rooted_tree, layout="rectangular", size=1) + geom_tiplab(align=FALSE, hjust=-.15, fontface=3) + xlim(0, 2.5) + geom_nodepoint(aes(subset = label > 75), size=3, shape=21, fill="darkgray", color="black")
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black")

MRCA(t2, "Sphagnum_compactum_4", "Sphagnum_fuscum_5")
t3 <- t2 %>% collapse(node=835) + geom_point2(aes(subset=(node==835)), shape=18, size=5, color="violetred4")
MRCA(t3, "Leucodon_julaceus_4", "Brachythecium_rivulare_4")
t4 <- t3 %>% collapse(node=683) + geom_point2(aes(subset=(node==683)), shape=18, size=5, color="lightseagreen")
MRCA(t4, "Climacium_americanum_2", "Platyhypnidium_sp._2")
t5 <- t4 %>% collapse(node=748) + geom_point2(aes(subset=(node==748)), shape=18, size=5, color="lightseagreen")
MRCA(t5, "Eosphagnum_inretortum_6", "Sphagnum_divinum_3")
t6 <- t5 %>% collapse(node=910) + geom_point2(aes(subset=(node==910)), shape=18, size=5, color="violetred4")
MRCA(t6, "Hygrohypnum_luridum_3", "Entodon_seductrix_4")
t7 <- t6 %>% collapse(node=611) + geom_point2(aes(subset=(node==611)), shape=18, size=5, color="lightseagreen")
MRCA(t7, "Fontinalis_antipyretica_2", "Campyliadelphus_chrysophyllus_2")
t8 <- t7 %>% collapse(node=520) + geom_point2(aes(subset=(node==520)), shape=18, size=5, color="lightseagreen")

pdf("RAxML_bipartitions.CESA_exon.collapsed.pdf", width=15, height=25)
t8 + geom_treescale(x=0, y=160)
dev.off()

tree <- read.tree("RAxML_bestTree.CESA_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node), size=3) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=501, resolve.root = TRUE, edgelabel=TRUE)
 
pdf("RAxML_bestTree.CESA_exon.pdf", width=15, height=60)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + geom_tiplab(align=FALSE, hjust=-.15, fontface=3) + xlim(0, 2.5) + geom_treescale(x=0, y=480)
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black")
t2
dev.off()

t <- ggtree(rooted_tree, layout="rectangular", size=1) + geom_tiplab(align=FALSE, hjust=-.15, fontface=3) + xlim(0, 2.5)
t2 <- t %<+% data + geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black")

MRCA(t2, "Sphagnum_compactum_4", "Sphagnum_fuscum_5")
t3 <- t2 %>% collapse(node=674) + geom_point2(aes(subset=(node==674)), shape=18, size=5, color="violetred4")
MRCA(t3, "Leucodon_julaceus_4", "Brachythecium_rivulare_4")
t4 <- t3 %>% collapse(node=522) + geom_point2(aes(subset=(node==522)), shape=18, size=5, color="lightseagreen")
MRCA(t4, "Climacium_americanum_2", "Platyhypnidium_sp._2")
t5 <- t4 %>% collapse(node=587) + geom_point2(aes(subset=(node==587)), shape=18, size=5, color="lightseagreen")
MRCA(t5, "Eosphagnum_inretortum_6", "Sphagnum_divinum_3")
t6 <- t5 %>% collapse(node=749) + geom_point2(aes(subset=(node==749)), shape=18, size=5, color="violetred4")
MRCA(t6, "Hygrohypnum_luridum_3", "Entodon_seductrix_4")
t7 <- t6 %>% collapse(node=941) + geom_point2(aes(subset=(node==941)), shape=18, size=5, color="lightseagreen")
MRCA(t7, "Fontinalis_antipyretica_2", "Campyliadelphus_chrysophyllus_2")
t8 <- t7 %>% collapse(node=853) + geom_point2(aes(subset=(node==853)), shape=18, size=5, color="lightseagreen")

pdf("RAxML_bestTree.CESA_exon.collapsed.pdf", width=15, height=25)
t8 + geom_treescale(x=0, y=160)
dev.off()


