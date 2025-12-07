library(phangorn)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(glue)

# Analysis for reference trees
setwd("/Users/corinna/Documents/Work/Schwartz_Lab/Plant_paralog_evolution/Bryophyta/CESA/CESA_trees_refs")

########################Tree with all reference genes (de novo assembled and BLASTx searches)
data <- read.table("tree_metadata.txt", header=T)
data2 <- data %>% 
  mutate(NewLab = ifelse(Fontface==3 & Info== "n.a.", glue("italic({Genus}~{Species})"), ifelse(Fontface==3 & Info!= "n.a.", glue("italic({Genus}~{Species})~{Info}"), ifelse(Fontface==4 & Info== "n.a.", glue("bolditalic({Genus}~{Species})"), glue("bolditalic({Genus}~{Species})~bold({Info})")))))

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

col <- c("Anthoceros" = safe_colorblind_palette[1], "Ceratodon" = safe_colorblind_palette[2], "Entodon" = safe_colorblind_palette[3], 
         "Fontinalis" = safe_colorblind_palette[4], "Hypnum" = safe_colorblind_palette[5], "Marchantia" = safe_colorblind_palette[6], 
         "Penium" = safe_colorblind_palette[7], "Physcomitrium" = safe_colorblind_palette[8], "Sphagnum" = safe_colorblind_palette[10], 
         "Takakia" = safe_colorblind_palette[11])

tree <- read.tree("RAxML_bipartitions.CESA_all_refs_exon.tre")

rooted_tree <- root(tree, node=109, resolve.root = TRUE, edgelabel=TRUE)

pdf("RAxML_bipartitions.CESA_all_refs_exon_v2.pdf", width=15, height=15)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_treescale(x=0, y=100) + 
  annotate("point", x=0, y=97, shape=21, fill="darkgray", color="black", size=3) + 
  annotate("text", x=0.04, y=97, label = "> 75% node support", hjust = "left") + 
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 75), size=2, shape=21, fill="darkgray", color="black")
t2 <- t %<+% data2 +
  theme(legend.position = "none") + 
  aes(color=factor(Genus)) + scale_color_manual(values = col, name="Genus", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.05, size=2, parse=T, family="Helvetica") + 
  geom_cladelab(node=116, label="CESA A", family="Helvetica", fontface="bold", offset=0.5) + 
  geom_cladelab(node=144, label="CESA B", family="Helvetica", fontface="bold", offset=0.6)
t2
dev.off()

ggsave("RAxML_CESA_all_refs_exon_v2.png", plot = t2, width = 6.5, height = 9, units = "in")

################ Analysis for full tree
setwd("/Users/corinna/Documents/Work/Schwartz_Lab/Plant_paralog_evolution/Bryophyta/CESA/CESA_trees_all_samples")

data <- read.table("tree_metadata.txt", header=T)
data2 <- data %>% mutate(NewLab = ifelse(Info== "n.a.", glue("italic({Genus}~{Species})"), glue("italic({Genus}~{Species})~{Info}")))

col <- c("Andreaeales" = safe_colorblind_palette[2], "Bartramiales" = safe_colorblind_palette[2], "Bryales" = safe_colorblind_palette[10], 
         "Buxbaumiales" = safe_colorblind_palette[4], "Dicranales" = safe_colorblind_palette[5], "Diphysciales" = safe_colorblind_palette[6], 
         "Encalyptales" = safe_colorblind_palette[7], "Grimmiales" = safe_colorblind_palette[1], "Orthotrichales" = safe_colorblind_palette[10], 
         "Polytrichales" = safe_colorblind_palette[1], "Pottiales" = safe_colorblind_palette[6], "Rhizogoniales" = safe_colorblind_palette[11], 
         "Scouleriales" = safe_colorblind_palette[4], "Tetraphidales" = safe_colorblind_palette[5], "Timmiales" = safe_colorblind_palette[6], 
         "Anthocerotales" = safe_colorblind_palette[1], "Pseudoditrichales" = safe_colorblind_palette[2], "Hypnales" = safe_colorblind_palette[3], 
         "Marchantiales" = safe_colorblind_palette[6], "Desmidiales" = safe_colorblind_palette[7], "Funariales" = safe_colorblind_palette[8], "Sphagnales" = safe_colorblind_palette[10], 
         "Takakiales" = safe_colorblind_palette[11])

tree <- read.tree("RAxML_bipartitions.CESA_exon.tre")

rooted_tree <- root(tree, node=880, resolve.root = TRUE, edgelabel=TRUE)
rooted_tree$edge.length[1] <- 0.1 #manually change root edge for visualization purposes

pdf("RAxML_bipartitions.CESA_exon_v2.pdf", width=25, height=60)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 1.3) + geom_treescale(x=0, y=530) + 
  geom_nodelab(aes(subset = !is.na(as.numeric(label))), geom="label", color="black", fill=NA, size=2, label.size=NA, nudge_x=-0.01, nudge_y=-0.55)

t2 <- t %<+% data2 +
  theme(legend.position = "none") + 
  aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.15, parse=T, family="Helvetica", size=2) + 
  geom_cladelab(node=614, label="CESA A", family="Helvetica", fontface="bold", offset=0.06, angle=90) + 
  geom_cladelab(node=616, label="CESA B", family="Helvetica", fontface="bold", offset=0.07, angle=90)
t2
dev.off()
ggsave("RAxML_allsamples.CESA_exon.png", plot = t2, width = 30, height = 30, units = "in")

###########################collapsed tree
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 1.5) + 
  annotate("point", x=0, y=180, shape=21, fill="darkgray", color="black", size=3) + 
  annotate("text", x=0.02, y=180, label = "> 75% node support", hjust = "left") + 
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 75), size=1.8, shape=21, fill="darkgray", color="black")
t2 <- t %<+% data2 +
  theme(legend.position = "none") + 
  aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.05, parse=T, family="Helvetica", size=1.8)

M1 <- MRCA(t2, "Eosphagnum_inretortum_6", "Sphagnum_magellanicum_2")
t3 <- t2 %>% collapse(node=M1) + geom_point2(aes(subset=(node==M1)), shape=21, size=2, fill=safe_colorblind_palette[10]) +
  geom_cladelabel(M1, "Sphagnales", color=safe_colorblind_palette[10], family="Helvetica", fontsize = 1.8)
M2 <- MRCA(t3, "Palustriella_falcata_4", "Entodon_seductrix_4")
t4 <- t3 %>% collapse(node=M2) + geom_point2(aes(subset=(node==M2)), shape=21, size=2, fill=safe_colorblind_palette[3]) + 
  geom_cladelabel(M2, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica", fontsize=1.8)
M3 <- MRCA(t4, "Pilotrichella_sp._3", "Campylium_stellatum_4")
t5 <- t4 %>% collapse(node=M3) + geom_point2(aes(subset=(node==M3)), shape=21, size=2, fill=safe_colorblind_palette[3]) + 
  geom_cladelabel(M3, "Hypnales + Leucobryum albidum", color=safe_colorblind_palette[3], family="Helvetica", fontsize=1.8) #Fix label in Illustrator
M4 <- MRCA(t5, "Sphagnum_compactum_4", "Sphagnum_capillifolium_5")
t6 <- t5 %>% collapse(node=M4) + geom_point2(aes(subset=(node==M4)), shape=21, size=2, fill=safe_colorblind_palette[10]) + 
  geom_cladelabel(M4, "Sphagnales", color=safe_colorblind_palette[10], family="Helvetica", fontsize=1.8)
M5 <- MRCA(t6, "Flatbergium_sericeum_5", "Sphagnum_subsecundum_3")
t7 <- t6 %>% collapse(node=M5) + geom_point2(aes(subset=(node==M5)), shape=21, size=2, fill=safe_colorblind_palette[10]) + 
  geom_cladelabel(M5, "Sphagnales", color=safe_colorblind_palette[10], family="Helvetica", fontsize=1.8)
M6 <- MRCA(t7, "Leucodon_julaceus_4", "Hylocomium_splendens_4")
t8 <- t7 %>% collapse(node=M6) + geom_point2(aes(subset=(node==M6)), shape=21, size=2, fill=safe_colorblind_palette[3]) + 
  geom_cladelabel(M6, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica", fontsize=1.8)
M7 <- MRCA(t8, "Loeskeobryum_brevirostre_12", "Platyhypnidium_sp._3")
t9 <- t8 %>% collapse(node=M7) + geom_point2(aes(subset=(node==M7)), shape=21, size=2, fill=safe_colorblind_palette[3]) + 
  geom_cladelabel(M7, "Hypnales + Leucobryum albidum", color=safe_colorblind_palette[3], family="Helvetica", fontsize=1.8) #Fix label in Illustrator

# find node for outgroup
M_out <- MRCA(t9, "Penium_margaritaceum_1","Penium_margaritaceum_2")
pdf("RAxML_bipartitions.CESA_exon.collapsed_v2.pdf", width=15, height=15)
t10 <- groupClade(t9, M_out) + aes(linetype=group) + geom_treescale(x=0, y=185) + 
  geom_cladelab(node=614, label="CESA A", family="Helvetica", fontface="bold", offset=0.1, angle = 90) + 
  geom_cladelab(node=616, label="CESA B", family="Helvetica", fontface="bold", offset=0.125, angle = 90)
t10
dev.off()

t10 + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"))
ggsave("RAxML_bipartitions.CESA_exon.collapsed.png", width = 6.5, height = 9, units = "in")


