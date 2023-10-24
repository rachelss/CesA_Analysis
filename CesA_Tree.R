#load necessary libraries to create the tree 

#these libraries are used for data manipulation
library (tidyverse)
library(dplyr)
library(glue)
#this library allows us to read in tree objects
library (treeio)
#these libraries are used in plotting
library (ggtree)
library(ggplot2)
library (ggrepel)




#Read in your RAXML output file: should be the bipartitionsBranchLabels file
file <- "RAxML_bipartitionsBranchLabels.all_genome_cesa_plusoutgroup_translatorx_aligned"
tree <-read.raxml(file)
tree

#Read in a metadata file containing information on species names
metadatafile <- "Species_labels_edited.csv"
metadata <-read.csv(metadatafile, header=TRUE, sep=",")
metadata



#transform the tree object into a datatype called a tibble that will allow you to manipulate 
#variables in the tree object
x <- as_tibble(tree)
x

#transform metadata object into a tibble
y <-as_tibble(metadata)
y

#join your two tibbles by the column "node"
z <-full_join(x,y, by ="node")
z

#save data as a csv
write.csv(z, "cesatreedata.csv")

#create variables from your tibble z that can be read by ggtree
tip.label <- z$label
bootstrap <-z$bootstrap
species <-z$species
cesa <-z$CesA 


#create a dataframe of your tip label, species, cesa variables 
d <- data.frame(label=tip.label, species=species,cesa=cesa)
d

#use dplyr to create d2, a new variable, as a function of the variables species and cesa
d2 <-dplyr::mutate(d,
                   lab = glue("{species} {cesa}"),)


#transform tibble z into a treedata object called tree1
tree1<- as.treedata(z)
tree1

#Use ggtree to create a figure from tree1 and d2, with branches colored by species
ggtree(tree1, aes(color=species), size=2) %<+% d2 + 
  geom_treescale(x=2, y=0) + #scale the figure
  geom_tiplab(aes(label=lab), hjust=-0.2) +  #add tip labels and adjust position
  geom_label(aes(label=bootstrap), size=3.5) + #add bootstrap labels to nodes
  xlim(0,3) + #define plot limits
  theme(legend.position='none') #remove the legend

#save the figure as a .png file for publication
ggsave("cesatree.png", width=10, height=10)


                           