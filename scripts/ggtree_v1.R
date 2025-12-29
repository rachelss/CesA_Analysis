library(phangorn)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(glue)

# Analysis for reference trees
setwd("/Users/corinna/Documents/Work/Schwartz_Lab/Plant_paralog_evolution/Bryophyta/CESA/CESA_trees_refs")

#Tree with all reference genes (de novo assembled and BLASTx searches)
data <- read.table("tree_metadata.txt", header=T)
data2 <- data %>% mutate(NewLab = ifelse(Fontface==3 & Info== "n.a.", glue("italic({Genus}~{Species})"), ifelse(Fontface==3 & Info!= "n.a.", glue("italic({Genus}~{Species})~{Info}"), ifelse(Fontface==4 & Info== "n.a.", glue("bolditalic({Genus}~{Species})"), glue("bolditalic({Genus}~{Species})~bold({Info})")))))

col <- c("Anthoceros" = "red2", "Ceratodon" = "cornflowerblue", "Entodon" = "mediumturquoise", "Fontinalis" = "darkcyan", "Hypnum" = "lightseagreen", "Marchantia" = "mediumpurple", "Penium" = "royalblue3", "Physcomitrium" = "plum", "Sphagnum" = "violetred4", "Takakia" = "darkgoldenrod1")

tree <- read.tree("RAxML_bipartitions.CESA_all_refs_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_text(aes(label=node)) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=109, resolve.root = TRUE, edgelabel=TRUE)

pdf("RAxML_bipartitions.CESA_all_refs_exon.pdf", width=15, height=15)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_treescale(x=0, y=100) + annotate("point", x=0, y=97, shape=21, fill="darkgray", color="black", size=3) + annotate("text", x=0.04, y=97, label = "> 75% node support", hjust = "left") + geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 75), size=3, shape=21, fill="darkgray", color="black")
t2 <- t %<+% data2 + geom_tippoint(aes(color=factor(Genus)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Genus)) + scale_color_manual(values = col, name="Genus", na.value="black") + geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.15, parse=T, family="Helvetica") + geom_cladelab(node=116, label="CESA A", family="Helvetica", fontface="bold", offset=0.5) + geom_cladelab(node=144, label="CESA B", family="Helvetica", fontface="bold", offset=0.6)
t2
dev.off()

# Analysis for full tree
setwd("/Users/corinna/Documents/Work/Schwartz_Lab/Plant_paralog_evolution/Bryophyta/CESA/CESA_trees_all_samples")

data <- read.table("tree_metadata.txt", header=T)
data2 <- data %>% mutate(NewLab = ifelse(Info== "n.a.", glue("italic({Genus}~{Species})"), glue("italic({Genus}~{Species})~{Info}")))

col <- c("Andreaeales" = "royalblue4", "Bartramiales" = "mediumturquoise", "Bryales" = "darkcyan", "Buxbaumiales" = "orange", "Dicranales" = "pink1", "Diphysciales" = "slateblue", "Encalyptales" = "plum4", "Grimmiales" = "palevioletred", "Orthotrichales" = "paleturquoise4", "Polytrichales" = "tomato", "Pottiales" = "hotpink", "Rhizogoniales" = "paleturquoise3", "Scouleriales" = "deeppink3", "Tetraphidales" = "gold", "Timmiales" = "salmon", "Anthocerotales" = "red2", "Pseudoditrichales" = "cornflowerblue", "Hypnales" = "lightseagreen", "Marchantiales" = "mediumpurple", "Desmidiales" = "royalblue3", "Funariales" = "plum", "Sphagnales" = "violetred4", "Takakiales" = "darkgoldenrod1")

tree <- read.tree("RAxML_bipartitions.CESA_exon.tre")

t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_text(aes(label=node), size=2) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=880, resolve.root = TRUE, edgelabel=TRUE)
 
pdf("RAxML_bipartitions.CESA_exon.pdf", width=25, height=60)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_treescale(x=0, y=530) + geom_nodelab(aes(subset = !is.na(as.numeric(label))), geom="label", color="black", fill="whitesmoke", size=3, label.size=NA, nudge_x=-0.05)
t2 <- t %<+% data2 + geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black") + geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.15, parse=T, family="Helvetica") + geom_cladelab(node=614, label="CESA A", family="Helvetica", fontface="bold", offset=0.2) + geom_cladelab(node=616, label="CESA B", family="Helvetica", fontface="bold", offset=0.25)
t2
dev.off()

t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + annotate("point", x=0, y=175, shape=21, fill="darkgray", color="black", size=3) + annotate("text", x=0.04, y=175, label = "> 75% node support", hjust = "left") + geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 75), size=3, shape=21, fill="darkgray", color="black")
t2 <- t %<+% data2 + geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black") + geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.15, parse=T, family="Helvetica")

MRCA(t2, "Eosphagnum_inretortum_6", "Sphagnum_magellanicum_2")
t3 <- t2 %>% collapse(node=588) + geom_point2(aes(subset=(node==588)), shape=18, size=5, color="violetred4") + geom_cladelabel(588, "Sphagnales", color="violetred4", family="Helvetica")
MRCA(t3, "Palustriella_falcata_4", "Entodon_seductrix_4")
t4 <- t3 %>% collapse(node=1027) + geom_point2(aes(subset=(node==1027)), shape=18, size=5, color="lightseagreen") + geom_cladelabel(1027, "Hypnales", color="lightseagreen", family="Helvetica")
MRCA(t4, "Pilotrichella_sp._3", "Campylium_stellatum_4")
t5 <- t4 %>% collapse(node=917) + geom_point2(aes(subset=(node==917)), shape=18, size=5, color="lightseagreen") + geom_cladelabel(917, "Hypnales + Leucobryum albidum", color="lightseagreen", family="Helvetica") # Fix Leucobryum label in Illustrator
MRCA(t5, "Sphagnum_compactum_4", "Sphagnum_capillifolium_5")
t6 <- t5 %>% collapse(node=619) + geom_point2(aes(subset=(node==619)), shape=18, size=5, color="violetred4") + geom_cladelabel(619, "Sphagnales", color="violetred4", family="Helvetica")
MRCA(t6, "Flatbergium_sericeum_5", "Sphagnum_subsecundum_3")
t7 <- t6 %>% collapse(node=817) + geom_point2(aes(subset=(node==817)), shape=18, size=5, color="violetred4") + geom_cladelabel(817, "Sphagnales", color="violetred4", family="Helvetica")
MRCA(t7, "Leucodon_julaceus_4", "Hylocomium_splendens_4")
t8 <- t7 %>% collapse(node=655) + geom_point2(aes(subset=(node==655)), shape=18, size=5, color="lightseagreen") + geom_cladelabel(655, "Hypnales", color="lightseagreen", family="Helvetica")
MRCA(t8, "Loeskeobryum_brevirostre_12", "Platyhypnidium_sp._3")
t9 <- t8 %>% collapse(node=758) + geom_point2(aes(subset=(node==758)), shape=18, size=5, color="lightseagreen") + geom_cladelabel(758, "Hypnales + Leucobryum albidum", color="lightseagreen", family="Helvetica") # Fix Leucobryum label in Illustrator

pdf("RAxML_bipartitions.CESA_exon.collapsed.pdf", width=22, height=27)
t9 + geom_treescale(x=0, y=180) + geom_cladelab(node=614, label="CESA A", family="Helvetica", fontface="bold", offset=0.23) + geom_cladelab(node=616, label="CESA B", family="Helvetica", fontface="bold", offset=0.28)
dev.off()
