library(phangorn)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(glue)

# Analysis for reference trees
#setwd("/Users/corinna/Documents/Work/Schwartz_Lab/Plant_paralog_evolution/Bryophyta/CESA/CESA_trees_refs")

################Tree with de novo assembled contigs for reference species - not using
data <- read.table("data/tree_metadata.txt", header=T)
data2 <- data %>% mutate(NewLab = ifelse(Info== "n.a.", glue("italic({Genus}~{Species})"), glue("italic({Genus}~{Species})~{Info}")))

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

col <- c("Anthoceros" = safe_colorblind_palette[1], "Ceratodon" = safe_colorblind_palette[2], "Entodon" = safe_colorblind_palette[3], 
         "Fontinalis" = safe_colorblind_palette[4], "Hypnum" = safe_colorblind_palette[5], "Marchantia" = safe_colorblind_palette[6], 
         "Penium" = safe_colorblind_palette[7], "Physcomitrium" = safe_colorblind_palette[8], "Sphagnum" = safe_colorblind_palette[10], 
         "Takakia" = safe_colorblind_palette[11])

tree <- read.tree("CESA_trees_refs/RAxML_bipartitions.CESA_ref_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node)) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=109, resolve.root = TRUE, edgelabel=TRUE)
 
#pdf("RAxML_bipartitions.CESA_ref_exon.pdf", width=12, height=10)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_treescale(x=0, y=53) + 
  geom_nodelab(aes(subset = !is.na(as.numeric(label))), geom="label", color="black", 
               fill=NA, size=2, label.size=NA, nudge_x=-0.05, nudge_y=-0.5)
t2 <- t %<+% data2 + #geom_tippoint(aes(color=factor(Genus)), shape=19, size=3) + 
  theme(legend.position = "none") + aes(color=factor(Genus)) + 
  scale_color_manual(values = col, name="Genus", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, parse=T, family="Helvetica", hjust=-.05, size=2) + 
  geom_cladelab(node=94, label="CESA A", family="Helvetica", fontface="bold", offset=0.65) + 
  geom_cladelab(node=74, label="CESA B", family="Helvetica", fontface="bold", offset=0.6)
t2
#dev.off()
ggsave("RAxML_CESA_ref_exon.png", plot = t2, width = 6.5, height = 6, units = "in")

########################ignoring this
tree <- read.tree("RAxML_bestTree.CESA_ref_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node)) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=110, resolve.root = TRUE, edgelabel=TRUE)
 
pdf("RAxML_bestTree.CESA_ref_exon.pdf", width=12, height=10)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_treescale(x=0, y=53)
t2 <- t %<+% data2 + geom_tippoint(aes(color=factor(Genus)), shape=19, size=3) + theme(legend.position = "none") + aes(color=factor(Genus)) + scale_color_manual(values = col, name="Genus", na.value="black") + geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.15, parse=T, family="Helvetica") + geom_cladelab(node=95, label="CESA A", family="Helvetica", fontface="bold", offset=0.65) + geom_cladelab(node=75, label="CESA B", family="Helvetica", fontface="bold", offset=0.75)
t2
dev.off()

########################Tree with all reference genes (de novo assembled and BLASTx searches)
data <- read.table("data/tree_metadata2.txt", header=T)
data2 <- data %>% 
  mutate(NewLab = ifelse(Fontface==3 & Info== "n.a.", glue("italic({Genus}~{Species})"), ifelse(Fontface==3 & Info!= "n.a.", glue("italic({Genus}~{Species})~{Info}"), ifelse(Fontface==4 & Info== "n.a.", glue("bolditalic({Genus}~{Species})"), glue("bolditalic({Genus}~{Species})~bold({Info})")))))

#col <- c("Anthoceros" = "red2", "Ceratodon" = "cornflowerblue", "Entodon" = "mediumturquoise", "Fontinalis" = "darkcyan", 
#"Hypnum" = "lightseagreen", "Marchantia" = "mediumpurple", "Penium" = "royalblue3", "Physcomitrium" = "plum", "Sphagnum" = "violetred4", 
#"Takakia" = "darkgoldenrod1")

tree <- read.tree("CESA_trees_refs/RAxML_bipartitions.CESA_all_refs_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node)) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=200, resolve.root = TRUE, edgelabel=TRUE)

#pdf("RAxML_bipartitions.CESA_all_refs_exon.pdf", width=15, height=15)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_treescale(x=0, y=100) + 
  annotate("point", x=0, y=97, shape=21, fill="darkgray", color="black", size=3) + 
  annotate("text", x=0.04, y=97, label = "> 75% node support", hjust = "left") + 
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 75), size=2, shape=21, fill="darkgray", color="black")
t2 <- t %<+% data2 + #geom_tippoint(aes(color=factor(Genus)), shape=19, size=3) + 
  theme(legend.position = "none") + 
  aes(color=factor(Genus)) + scale_color_manual(values = col, name="Genus", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.05, size=2, parse=T, family="Helvetica") + 
  geom_cladelab(node=168, label="CESA A", family="Helvetica", fontface="bold", offset=0.5) + 
  geom_cladelab(node=114, label="CESA B", family="Helvetica", fontface="bold", offset=0.6)
t2
#dev.off()
ggsave("RAxML_CESA_all_refs_exon.png", plot = t2, width = 6.5, height = 9, units = "in")

################ Analysis for full tree
#setwd("/Users/corinna/Documents/Work/Schwartz_Lab/Plant_paralog_evolution/Bryophyta/CESA/CESA_trees_all_samples")

data <- read.table("data/tree_metadata3.txt", header=T)
data2 <- data %>% mutate(NewLab = ifelse(Info== "n.a.", glue("italic({Genus}~{Species})"), glue("italic({Genus}~{Species})~{Info}")))

#col <- c("Andreaeales" = "royalblue4", "Bartramiales" = "mediumturquoise", "Bryales" = "darkcyan", "Buxbaumiales" = "orange", "Dicranales" = "pink1", 
#         "Diphysciales" = "slateblue", "Encalyptales" = "plum4", "Grimmiales" = "palevioletred", "Orthotrichales" = "paleturquoise4", 
#         "Polytrichales" = "tomato", "Pottiales" = "hotpink", "Rhizogoniales" = "paleturquoise3", "Scouleriales" = "deeppink3", "Tetraphidales" = "gold", 
#         "Timmiales" = "salmon", "Anthocerotales" = "red2", "Pseudoditrichales" = "cornflowerblue", "Hypnales" = "lightseagreen", 
#         "Marchantiales" = "mediumpurple", "Desmidiales" = "royalblue3", "Funariales" = "plum", "Sphagnales" = "violetred4", "Takakiales" = "darkgoldenrod1")
#col <- c("Anthoceros" = safe_colorblind_palette[1], "Ceratodon" = safe_colorblind_palette[2], "Entodon" = safe_colorblind_palette[3], 
#"Fontinalis" = safe_colorblind_palette[4], "Hypnum" = safe_colorblind_palette[5], "Marchantia" = safe_colorblind_palette[6], 
#"Penium" = safe_colorblind_palette[7], "Physcomitrium" = safe_colorblind_palette[8], "Sphagnum" = safe_colorblind_palette[10], 
#"Takakia" = safe_colorblind_palette[11])
col <- c("Andreaeales" = safe_colorblind_palette[2], "Bartramiales" = safe_colorblind_palette[2], "Bryales" = safe_colorblind_palette[10], 
         "Buxbaumiales" = safe_colorblind_palette[4], "Dicranales" = safe_colorblind_palette[5], "Diphysciales" = safe_colorblind_palette[6], 
         "Encalyptales" = safe_colorblind_palette[7], "Grimmiales" = safe_colorblind_palette[1], "Orthotrichales" = safe_colorblind_palette[10], 
         "Polytrichales" = safe_colorblind_palette[1], "Pottiales" = safe_colorblind_palette[6], "Rhizogoniales" = safe_colorblind_palette[11], 
         "Scouleriales" = safe_colorblind_palette[4], "Tetraphidales" = safe_colorblind_palette[5], "Timmiales" = safe_colorblind_palette[6], 
         "Anthocerotales" = safe_colorblind_palette[1], "Pseudoditrichales" = safe_colorblind_palette[2], "Hypnales" = safe_colorblind_palette[3], 
         "Marchantiales" = safe_colorblind_palette[6], "Desmidiales" = safe_colorblind_palette[7], "Funariales" = safe_colorblind_palette[8], "Sphagnales" = safe_colorblind_palette[10], 
         "Takakiales" = safe_colorblind_palette[11])

tree <- read.tree("CESA_trees_all_samples/RAxML_bipartitions.CESA_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node), size=3) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=662, resolve.root = TRUE, edgelabel=TRUE) %>% drop.tip(c("Flowersia_sinensis_1", "Flowersia_sinensis_3")) #these are other plants
rooted_tree$edge.length[1] <- 0.1 #manually change root edge for visualization purposes

#pdf("RAxML_bipartitions.CESA_exon.pdf", width=25, height=60)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 1.3) + geom_treescale(x=0, y=480) + 
  geom_nodelab(aes(subset = !is.na(as.numeric(label))), geom="label", color="black", fill=NA, size=2, label.size=NA, nudge_x=-0.01, nudge_y=-0.55)

MRCA_CESAA <- MRCA(t, "Takakia_lepidozioides_2", "Andreaea_rupestris" )
MRCA_CESAB <- MRCA(t, "Physcomitrium_patens_1", "Sphagnum_portoricense_1"  )

t2 <- t %<+% data2 + #geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + 
  theme(legend.position = "none") + 
  aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.15, parse=T, family="Helvetica", size=2) + 
  geom_cladelab(node=MRCA_CESAA, label="CESA A", family="Helvetica", fontface="bold", offset=0.05, angle=90) + 
  geom_cladelab(node=MRCA_CESAB, label="CESA B", family="Helvetica", fontface="bold", offset=0.05, angle=90)
t2
#dev.off()
ggsave("RAxML_allsamples.CESA_exon.png", plot = t2, width = 30, height = 30, units = "in")

###########################collapsed tree
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 1.5) + 
  annotate("point", x=0, y=155, shape=21, fill="darkgray", color="black", size=3) + 
  annotate("text", x=0.04, y=155, label = "> 75% node support", hjust = "left") + 
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) > 75), size=1.8, shape=21, fill="darkgray", color="black")
t2 <- t %<+% data2 + #geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + 
  theme(legend.position = "none") + 
  aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.05, parse=T, family="Helvetica", size=1.8)

M1 <- MRCA(t2, "Sphagnum_compactum_4", "Sphagnum_fuscum_5")
t3 <- t2 %>% collapse(node=M1) + geom_point2(aes(subset=(node==M1)), shape=18, size=2, color=safe_colorblind_palette[10]) + #835
  geom_cladelabel(M1, "Sphagnales", color=safe_colorblind_palette[10], family="Helvetica", fontsize = 1.8)
M2 <- MRCA(t3, "Leucodon_julaceus_4", "Brachythecium_rivulare_4")
t4 <- t3 %>% collapse(node=M2) + geom_point2(aes(subset=(node==M2)), shape=18, size=2, color=safe_colorblind_palette[3]) + 
  geom_cladelabel(M2, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica", fontsize=1.8)
M3 <- MRCA(t4, "Climacium_americanum_2", "Platyhypnidium_sp._2")
t5 <- t4 %>% collapse(node=M3) + geom_point2(aes(subset=(node==M3)), shape=18, size=2, color=safe_colorblind_palette[3]) + 
  geom_cladelabel(M3, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica", fontsize=1.8)
M4 <- MRCA(t5, "Eosphagnum_inretortum_6", "Sphagnum_divinum_3")
t6 <- t5 %>% collapse(node=M4) + geom_point2(aes(subset=(node==910)), shape=18, size=2, color=safe_colorblind_palette[10]) + 
  geom_cladelabel(M4, "Sphagnales", color=safe_colorblind_palette[10], family="Helvetica", fontsize=1.8)
M5 <- MRCA(t6, "Hygrohypnum_luridum_3", "Entodon_seductrix_4")
t7 <- t6 %>% collapse(node=M5) + geom_point2(aes(subset=(node==611)), shape=18, size=2, color=safe_colorblind_palette[3]) + 
  geom_cladelabel(M5, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica", fontsize=1.8)
M6 <- MRCA(t7, "Fontinalis_antipyretica_2", "Campyliadelphus_chrysophyllus_2")
t8 <- t7 %>% collapse(node=M6) + geom_point2(aes(subset=(node==520)), shape=18, size=2, color=safe_colorblind_palette[3]) + 
  geom_cladelabel(M6, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica", fontsize=1.8)
t8

# find node for outgroup
M_out <- MRCA(t8, "Penium_margaritaceum_1","Penium_margaritaceum_2")
t9 <- groupClade(t8, M_out) + aes(linetype=group) + geom_treescale(x=0, y=160) + 
  geom_cladelab(node=MRCA_CESAA, label="CESA A", family="Helvetica", fontface="bold", offset=0.205, angle = 90) + 
  geom_cladelab(node=MRCA_CESAB, label="CESA B", family="Helvetica", fontface="bold", offset=0.26, angle = 90)

#the following should help but shortening the one branch makes the other longer and the ingroup brnaches don't change because x pos don't change
#outgroup_xval <- filter(t9$data, !is.na(x) & isTip==TRUE) %>% pull(x) %>% mean()
#t9$data[t9$data$node == M_out, "x"] <- 0.2 #change branch length
#t9$data[t9$data$node == 662, "x"] <- 0.2 #change branch length
#t9$data[t9$data$node == 163, "x"] <- 0.2 #change Penium position
#t9$data[t9$data$node == 164, "x"] <- 0.2 #change Penium position
t9 + theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"))

#pdf("RAxML_bipartitions.CESA_exon.collapsed.pdf", width=22, height=25)
ggsave("RAxML_bipartitions.CESA_exon.collapsed.png", width = 6.5, height = 9, units = "in")
#dev.off()

############################
tree <- read.tree("RAxML_bestTree.CESA_exon.tre")

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + geom_text(aes(label=node), size=3) + geom_tiplab(align=TRUE, hjust=-.15)
t

rooted_tree <- root(tree, node=501, resolve.root = TRUE, edgelabel=TRUE)
 
pdf("RAxML_bestTree.CESA_exon.pdf", width=25, height=60)
t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5) + geom_treescale(x=0, y=480)
t2 <- t %<+% data2 + geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + theme(legend.position = "none") + 
  aes(color=factor(Order)) + scale_color_manual(values = col, name="Order", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.15, parse=T, family="Helvetica") + 
  geom_cladelab(node=497, label="CESA A", family="Helvetica", fontface="bold", offset=0.2) + 
  geom_cladelab(node=506, label="CESA B", family="Helvetica", fontface="bold", offset=0.25)
t2
dev.off()

t <- ggtree(rooted_tree, layout="rectangular", size=1) + xlim(0, 2.5)
t2 <- t %<+% data2 + geom_tippoint(aes(color=factor(Order)), shape=19, size=3) + 
  theme(legend.position = "none") + aes(color=factor(Order)) + 
  scale_color_manual(values = col, name="Order", na.value="black") + 
  geom_tiplab(aes(label=NewLab), align=FALSE, hjust=-.15, parse=T, family="Helvetica")

MRCA(t2, "Sphagnum_compactum_4", "Sphagnum_fuscum_5")
t3 <- t2 %>% collapse(node=674) + geom_point2(aes(subset=(node==674)), shape=18, size=5, color=safe_colorblind_palette[10]) + 
  geom_cladelabel(674, "Sphagnales", color=safe_colorblind_palette[10], family="Helvetica")
MRCA(t3, "Leucodon_julaceus_4", "Brachythecium_rivulare_4")
t4 <- t3 %>% collapse(node=522) + geom_point2(aes(subset=(node==522)), shape=18, size=5, color=safe_colorblind_palette[3]) + 
  geom_cladelabel(522, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica")
MRCA(t4, "Climacium_americanum_2", "Platyhypnidium_sp._2")
t5 <- t4 %>% collapse(node=587) + geom_point2(aes(subset=(node==587)), shape=18, size=5, color=safe_colorblind_palette[3]) + 
  geom_cladelabel(587, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica")
MRCA(t5, "Eosphagnum_inretortum_6", "Sphagnum_divinum_3")
t6 <- t5 %>% collapse(node=749) + geom_point2(aes(subset=(node==749)), shape=18, size=5, color=safe_colorblind_palette[10]) + 
  geom_cladelabel(749, "Sphagnales", color=safe_colorblind_palette[10], family="Helvetica")
MRCA(t6, "Hygrohypnum_luridum_3", "Entodon_seductrix_4")
t7 <- t6 %>% collapse(node=941) + geom_point2(aes(subset=(node==941)), shape=18, size=5, color=safe_colorblind_palette[3]) + 
  geom_cladelabel(941, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica")
MRCA(t7, "Fontinalis_antipyretica_2", "Campyliadelphus_chrysophyllus_2")
t8 <- t7 %>% collapse(node=853) + geom_point2(aes(subset=(node==853)), shape=18, size=5, color=safe_colorblind_palette[3]) + 
  geom_cladelabel(853, "Hypnales", color=safe_colorblind_palette[3], family="Helvetica")

pdf("RAxML_bestTree.CESA_exon.collapsed.pdf", width=22, height=25)
t8 + geom_treescale(x=0, y=160) + geom_cladelab(node=497, label="CESA A", family="Helvetica", fontface="bold", offset=0.23) + 
  geom_cladelab(node=506, label="CESA B", family="Helvetica", fontface="bold", offset=0.28)
dev.off()
