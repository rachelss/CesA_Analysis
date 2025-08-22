library(phangorn)
library(ggtree)
library(ggplot2)
library(ape)
library(dplyr)
library(glue)
library(treeio)

tree <- treeio::read.beast.newick("dnds_analysis/Annotated_CesA_Tree")
#print(as.phylo(tree)$tip.label %>% sort())

data <- read.table("data/tree_metadata4.txt", header=T)
data2 <- data %>% mutate(NewLab = ifelse(Info== "n.a.", glue("italic({Genus}~{Species})"), glue("italic({Genus}~{Species})~{Info}")))

t <- ggtree(tree, layout="rectangular", size=1, branch.length="none") + 
  hexpand(.5) + 
  geom_nodelab(aes(label=`!name`), hjust=1.1, vjust=-.23, size=2.5)

t2 <- t %<+% data2 +
  geom_tiplab(aes(label=NewLab), align=FALSE, parse=T, family="Helvetica", size=2)

t2
ggsave("dnds_phylo.png", plot = t2, width = 6.5, height = 4.5, units = "in")
