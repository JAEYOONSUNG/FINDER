#' Visualizing the results of \code{bouth} with GraPhlAn
#' 
#' Create a tree figure for visualizing the results of \code{bouth}
#' 
#' @param bouth.out an output object from the \code{bouth} function.
#' @param show.leaf a logical value indicating whether or not the leaf nodes are
#'   shown in the output figure.
#' @param output.dir the directory for storing the output figure. The default is
#'   the current working directory.
#' @param graphlan.dir the directory where the GraPhlAn package is located.
#'   
#' @note This function is dependent on the python package GraPhlAn. Parameters
#'   \code{clade.separation}, \code{branch.thickness},
#'   \code{branch.bracket.depth}, \code{branch.bracket.width}, and
#'   \code{clade.marker.size} can be set the same as in the python package. Details can be found
#'   at
#'   (\href{https://bitbucket.org/nsegata/graphlan/overview}{https://bitbucket.org/nsegata/graphlan/overview}).
#'   
#'   
#' @references Asnicar, F, et al. "Compact graphical representation of
#'   phylogenetic data and metadata with GraPhlAn." PeerJ 3 (2015): e1029.
#'   
#' @examples
#' data(IBD)
#' 
#' ## performing the one-stage, weighted bottom-up test on the IBD data
#' test.1 = bouth(anno.table = IBD$tax.table, pvalue.leaves = IBD$pvalue.otus,
#' na.symbol = "unknown", far = 0.1, is.weighted = TRUE)
#' 
#' 
#' ## suppose the GraPhlAn package is located at graphlan_directory/
#' graphlan(bouth.out = test.1, graph.dir = graphlan_directory)
#'
#' @export
graphlan<-function(bouth.out, show.leaf = FALSE, output.dir = getwd(), graphlan.dir = NULL,
                   clade.separation = NULL, branch.thickness = NULL, branch.bracket.depth = NULL, branch.bracket.width= NULL,
                   clade.marker.size = NULL){
  
  if(is.null(graphlan.dir)){stop("Unspecified directory for GraPhlAn.")}
  if(!"graphlan"%in%dir(graphlan.dir)){stop("The GraphlAn hasn't been installed at this directory.")}
  if(is.null(clade.separation)) clade.separation = 0.5
  if(is.null(branch.thickness)) branch.thickness = 0.5
  if(is.null(branch.bracket.depth)) branch.bracket.depth = 0.5
  if(is.null(branch.bracket.width)) branch.bracket.width = 0.25
  if(is.null(clade.marker.size)) clade.marker.size = 15
  basics = c(paste("clade_separation",clade.separation,sep='\t'),
             paste("branch_thickness",branch.thickness,sep='\t'),
             paste("branch_bracket_depth",branch.bracket.depth, sep='\t'),
             paste("branch_bracket_width", branch.bracket.width, sep = '\t'),
             paste("clade_marker_size", clade.marker.size, sep='\t'),
             paste("clade_marker_edge_color", "#555555", sep='\t'),
             paste("clade_marker_edge_width", 0.5, sep='\t'),
             paste("branch_color_from_ancestor", 0, sep='\t'))
  gamma_prop1 = 0.02 # A subtree will be colored if it contains more than 2% of all leaf nodes
  gamma_prop2 = 0.1  # A subtree will have a label if it contains more than 10% of all leaf nodes
  tree = bouth.out$tree
  nn = nrow(tree@node)
  nl = sum(tree@node$level==1)
  nr = max(tree@node$level)
  if(!show.leaf){
    nr = nr - 1
  }
  ## create styles and obtain the mapping table
  sub.nodes = tree@node$num[tree@node$level == max(tree@node$level)-1]
  sub.count = array(0, dim = length(sub.nodes))
  out0 = NULL
  
  # the mapping table
  table = matrix(NA, nrow = nl, ncol = nr)
  for(i in 1:nl){
    stack = NULL; len.stack = 0
    curr = i
    while(curr <= nn){
      stack = c(curr, stack)
      len.stack = len.stack+1
      curr = tree@node$parent[curr]
    }
    if(!show.leaf & stack[len.stack] <= nl){
      stack = stack[-len.stack]
      len.stack = len.stack-1
    }
    table[i,1:len.stack] = tree@node$id[stack]
    sub.count[sub.nodes%in%stack] = sub.count[sub.nodes%in%stack] + 1
  }
  # colors for subtrees
  sub.nodes = sub.nodes[sub.count >= nl*gamma_prop1]
  sub.count = sub.count[sub.count >= nl*gamma_prop1]
  sub.nodes = tree@node$id[sub.nodes]
  
  mypalette = RColorBrewer::brewer.pal(length(sub.nodes),"Accent")
  for(i in 1: length(sub.nodes)){
    out0 = c(out0, paste(sub.nodes[i], "annotation_background_color", mypalette[i], sep = '\t'))
    if(sub.count[i] > gamma_prop2){
      out0 = c(out0, paste(sub.nodes[i], "annotation", sub.nodes[i], sep = '\t'))
    }
  }
  
  table = data.frame(table, stringsAsFactors = FALSE)
  colnames(table) = paste("attri",c(1:nr),sep='_')
  # unique lineages
  table = dplyr::distinct_(table)
  # prepare for the annotation I
  out1 = array(NA, nrow(table))
  for(i in 1:nrow(table)){
    out1[i] = paste(table[i,!is.na(table[i,])], collapse ='.')
  }
  
  ## label detected nodes
  sig.id = which(tree@test$reject)
  if(!show.leaf){
    sig.id = sig.id[sig.id > nl]
  }
  sig.id = tree@node$id[sig.id]
  out2 = array(NA, length(sig.id))
  for(i in 1:length(sig.id)){
    out2[i] = paste(sig.id[i], "clade_marker_color", "r", sep='\t')
  }
  
  ## output all files for GraPhlAn
  anno1 = paste(output.dir,"anno1.txt", sep = '/')
  write.table(c(basics,out0, out2), anno1, col.names = F, row.names = F, quote = F)
  anno2 = paste(output.dir,"anno2.txt", sep = '/')
  write.table(out1, anno2, col.names = F, row.names = F, quote = F)
  
  command = c("#!/bin/sh", paste("export PATH=", graphlan.dir ,"/graphlan/:$PATH",sep=''),
              "graphlan_annotate.py anno2.txt plot.xml --annot anno1.txt",
              "graphlan.py plot.xml plot.png --dpi 150 --size 4 --pad 0.2")
  
  script = paste(output.dir,"run.sh",sep='/')
  write.table(command, script, col.names = F, row.names = F, quote = F)
  system2("bash", args = "run.sh")
  file.remove(anno1)
  file.remove(anno2)
  file.remove(script)
}

graphlan(clade.marker.size = 3,clade.separation = 0.35, branch.bracket.depth = 0.5, branch.thickness = 1.2)


##GraPhlAn##
pacman::p_load(dplyr, plyr, tidyr, ggplot2, ggnewscale, ggpubr, GGally, cowplot, stringr, stringi, ggnewscale, gghighlight, scales, reshape, reshape2, gtools, magrittr, rlist, xlsx, openxlsx, taxa, taxize, myTAI, pipeR, BBmisc, data.table) 

ezcloud_QIIME<-read.csv("ezbiocloud_id_taxonomy.txt", sep="\t", header = FALSE)
ezcloud_DB<-separate(ezcloud_QIIME, col=V2, into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";")
ezcloud_DB<-ezcloud_DB[,-1]

ezcloud_Phylum_DB<-unique(setNames(as.data.frame(paste(ezcloud_DB$Phylum,sep=".")), "Phylum"))
ezcloud_Class_DB<-unique(setNames(as.data.frame(paste(ezcloud_DB$Phylum,ezcloud_DB$Class,sep=".")), "Class"))
ezcloud_Order_DB<-unique(setNames(as.data.frame(paste(ezcloud_DB$Phylum,ezcloud_DB$Class,ezcloud_DB$Order,sep=".")), "Order"))
ezcloud_Family_DB<-unique(setNames(as.data.frame(paste(ezcloud_DB$Phylum,ezcloud_DB$Class,ezcloud_DB$Order,ezcloud_DB$Family,sep=".")), "Family"))
ezcloud_Genus_DB<-unique(setNames(as.data.frame(paste(ezcloud_DB$Phylum,ezcloud_DB$Class,ezcloud_DB$Order,ezcloud_DB$Family,ezcloud_DB$Genus,sep=".")), "Genus"))

names(ezcloud_Phylum_DB)[1]<-"taxa"
names(ezcloud_Class_DB)[1]<-"taxa"
names(ezcloud_Order_DB)[1]<-"taxa"
names(ezcloud_Family_DB)[1]<-"taxa"
names(ezcloud_Genus_DB)[1]<-"taxa"

ezcloud_taxa_DB <- rbind(ezcloud_Phylum_DB,ezcloud_Class_DB,ezcloud_Order_DB,ezcloud_Family_DB,ezcloud_Genus_DB)
row.names(ezcloud_taxa_DB)<-NULL
ezcloud_taxa_DB<-separate(ezcloud_taxa_DB, col="taxa", into=c("supertaxa","taxa"),sep="\\.(?=[^\\.]+$)",fill="left", extra="merge")



####
require(xlsx)    

temp = list.files(pattern="*.tax")
list2env(
  lapply(setNames(temp, make.names(gsub("*.tax$", "", temp))), 
         read.csv, sep="\t", header=FALSE, 
         stringsAsFactors=FALSE, na.strings=c(""," ","NA")), envir = .GlobalEnv)

#list<-ls(all.names=T,envir=globalenv())
list<-gsub("*.tax$", "", temp)
list<-gsub("-", ".", list)
#list<-list[-c(24,48)]
taxa<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


for(i in 1:length(list)){
  assign(list[i],rename(get(list[i]) ,c("V1" = "Kingdom","V2" = "Phylum","V3" = "Class","V4" = "Order","V5" = "Family","V6" = "Genus","V7" = "Species")))
}
#for(i in 1: length(list)){
#  assign(list[i],lapply(get((list)[i]), setNames, taxa))
#}
#-lapply(get(list), setNames, nm = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
{
for(i in 1:length(list)){
  assign(list[i],separate(get(list[i]) ,
                          col= c("Kingdom"), 
                          into = c("Kingdom", "No. of Kingdom"), 
                          sep="\\: ", extra = "merge"))
}

for(i in 1:length(list)){
  assign(list[i],separate(get(list[i]) ,
                          col= c("Phylum"), 
                          into = c("Phylum", "No. of Phylum"), 
                          sep="\\: ", extra = "merge"))
}

for(i in 1:length(list)){
  assign(list[i],separate(get(list[i]) ,
                          col= c("Class"), 
                          into = c("Class", "No. of Class"), 
                          sep="\\: ", extra = "merge"))
}

for(i in 1:length(list)){
  assign(list[i],separate(get(list[i]) ,
                          col= c("Order"), 
                          into = c("Order", "No. of Order"), 
                          sep="\\: ", extra = "merge"))
}

for(i in 1:length(list)){
  assign(list[i],separate(get(list[i]) ,
                          col= c("Family"), 
                          into = c("Family", "No. of Family"), 
                          sep="\\: ", extra = "merge"))
}

for(i in 1:length(list)){
  assign(list[i],separate(get(list[i]) ,
                          col= c("Genus"), 
                          into = c("Genus", "No. of Genus"), 
                          sep="\\: ", extra = "merge"))
}


for(i in 1:length(list)){
  assign(list[i],separate(get(list[i]) ,
                          col= c("Species"), 
                          into = c("Species", "No. of Species"), 
                          sep="\\: ", extra = "merge"))
}
}

for(i in 1:length(list)){
  assign(list[i],transform(get(list[i]), 
                               "No. of Kingdom" = as.numeric(`No. of Kingdom`),
                               "No. of Phylum" = as.numeric(`No. of Phylum`), 
                               "No. of Class" = as.numeric(`No. of Class`), 
                               "No. of Order" = as.numeric(`No. of Order`), 
                               "No. of Family" = as.numeric(`No. of Family`), 
                               "No. of Genus" = as.numeric(`No. of Genus`), 
                               "No. of Species" = as.numeric(`No. of Species`)))}


for(i in 1:length(list)){
  assign(list[i],rename(get(list[i]) ,c("No..of.Kingdom" = "No. of Kingdom",
                                        "No..of.Phylum" = "No. of Phylum",
                                        "No..of.Class" = "No. of Class",
                                        "No..of.Order" = "No. of Order",
                                        "No..of.Family" = "No. of Family",
                                        "No..of.Genus" = "No. of Genus",
                                        "No..of.Species" = "No. of Species")))
}

for(i in 1:length(list)){
  assign(paste(list[i],"taxa",sep="_"),
rbind(
  melt(get(list[i])) %>% subset(!is.na(value)) %>% filter(variable == "No. of Phylum") %>% select("Phylum", "variable", "value") %>% set_names(c("taxa","variable","value")),
  melt(get(list[i])) %>% subset(!is.na(value)) %>% filter(variable == "No. of Class") %>% select("Class", "variable", "value") %>% set_names(c("taxa","variable","value")),
  melt(get(list[i])) %>% subset(!is.na(value)) %>% filter(variable == "No. of Order") %>% select("Order", "variable", "value") %>% set_names(c("taxa","variable","value")),
  melt(get(list[i])) %>% subset(!is.na(value)) %>% filter(variable == "No. of Family") %>% select("Family", "variable", "value") %>% set_names(c("taxa","variable","value")),
  melt(get(list[i])) %>% subset(!is.na(value)) %>% filter(variable == "No. of Genus") %>% select("Genus", "variable", "value") %>% set_names(c("taxa","variable","value"))
))}


##hierarchical taxa assignment_supertax##
#for(j in 1:length(list)){
#for(j in 28:40){
#  assign(paste(list[j],"taxa_hier",sep="_"), taxize::classification(get(paste(list[j],"taxa",sep="_"))$taxa, db = 'ncbi', return_id=FALSE, rows=1))
#  }
#for(i in 1:1){
##for(i in 1:length(get(paste(list[j],"taxa_hier",sep="_")))){
#    paste(
#      lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "phylum")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#      lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "class")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#      lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "order")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#      lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "family")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#      lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "genus")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#      sep = ".")
#}

#for(i in 1:1){
#  #for(i in 1:length(get(paste(list[j],"taxa_hier",sep="_")))){
#  paste(
#    lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "phylum")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "class")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "order")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "family")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    lapply(get(paste(list[j],"taxa_hier",sep="_"))[i], function(x) filter(x, rank == "genus")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    sep = ".")
#}


#for(j in 1){
#  taxize::classification(get(paste(list[j],"taxa",sep="_"))$taxa, db = 'ncbi', return_id=FALSE) %>% mutate(taxa_hier =paste(
#    lapply(z[i], function(x) filter(x, rank == "phylum")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    lapply(z[i], function(x) filter(x, rank == "class")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    lapply(z[i], function(x) filter(x, rank == "order")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    lapply(z[i], function(x) filter(x, rank == "family")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    lapply(z[i], function(x) filter(x, rank == "genus")) %>% list.select(name) %>% as.data.frame() %>%t() %>% as.character(),
#    sep = "."))}


################################################################################################################################################################################
################################################################################################################################################################################


for(i in 1:length(list)){
  assign(paste(list[i],"taxa",sep="_"), merge(get(paste(list[i],"taxa",sep="_")), ezcloud_taxa_DB, all.x = TRUE, by= "taxa"))
}


##############check
#for(i in 1:length(list)){
#  assign(paste(list[i],"taxa",sep="_"), get(paste(list[i],"taxa",sep="_")) %>% unite(supertaxa=supertaxa,get(paste(list[i],"taxa",sep="_")) %>% filter(`variable` == "No. of Phylum") %>% select(taxa)))
#}




for(i in 1:length(list)){
  assign(paste(list[i],"taxa",sep="_"), get(paste(list[i],"taxa",sep="_")) %>% left_join( get(paste(list[i],"taxa",sep="_")) %>% filter(is.na(supertaxa) & variable=="No. of Phylum") %>% mutate("supertaxa"=paste0(taxa)),by = `taxa`))
}
get(paste(list[i],"taxa",sep="_")) %>% filter(is.na(supertaxa) & variable=="No. of Phylum") %>% mutate("supertaxa"=paste0(taxa))

for(i in 1:length(list)){
  assign(paste(list[i],"taxa",sep="_"), merge(get(paste(list[i],"taxa",sep="_")), get(paste(list[i],"taxa",sep="_")) %>% filter(is.na(supertaxa) & variable=="No. of Phylum") %>% mutate("supertaxa"=paste0(taxa)), all.x = TRUE, by= "taxa"))
}




for(i in 1:length(list)){
  assign(paste(list[i],"taxa",sep="_"), get(paste(list[i],"taxa",sep="_")) %>% left_join(get(paste(list[i],"taxa",sep="_")) %>% filter(is.na(supertaxa) & variable=="No. of Phylum") %>% mutate("supertaxa"=paste0(taxa)),by=c("taxa","variable","value")))
}











for(i in 1:length(list)){
  assign(paste(list[i]), get(paste(list[i])) %>% mutate(`clade marker size` = value/`Total No. of taxa`))
}

for(i in 1:length(list)){
  assign(paste(list[i]), get(paste(list[i])) %>% mutate(`normalized size` = normalize(`clade marker size`, method = "range", range = c(10, 60), margin = 1, on.constant = "quiet")))
}

for(i in 1:length(list)){
  assign(paste(list[i]), get(paste(list[i])) %>% mutate(`for_plot` = paste(get(paste(list[i]))$taxa, "clade_marker_size",get(paste(list[i]))$`normalized size`, sep="\t")))
}

for(i in 1:length(list)){
  assign(list[i],separate(get(list[i]) ,
                          col= "variable", 
                          into = c("variable2","Rank"),
                          sep="\\ (?=[^\\ ]+$)", extra = "drop", fill="right", remove = FALSE))
}

for(i in 1:length(list)){
  assign(list[i],get(list[i])[,-2])
}

##delta species
for(i in 1:length(list)){
  assign(paste(list[i],"for_plot",sep="_"),
         rbind(
           get(paste(list[i])) %>% filter(rank=="Phylum"),
           get(paste(list[i])) %>% filter(rank=="Class"),
           get(paste(list[i])) %>% filter(rank=="Order"),
           get(paste(list[i])) %>% filter(rank=="Family"),
           get(paste(list[i])) %>% filter(rank=="Genus"),
           get(paste(list[i])) %>% filter(rank=="Species")
         ))}



## annotation assignment ##
for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_annotation",sep="_"),
         rbind(
           get(paste(list[i],"for_plot",sep="_")) %>% filter(rank=="Phylum") %>% mutate(annotation="*"),
           get(paste(list[i],"for_plot",sep="_")) %>% filter(rank=="Class") %>% mutate(annotation="*:*")
         ))}

for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_annotation",sep="_"), get(paste(list[i],"for_plot_annotation",sep="_")) %>% mutate(`for_plot_annotation` = paste(get(paste(list[i],"for_plot_annotation",sep="_"))$taxa, "annotation",get(paste(list[i],"for_plot_annotation",sep="_"))$`annotation`, sep="\t")))
}


for(i in 1:length(list)){
  assign(paste(list[i],"for_plot",sep="_"),separate(get(paste(list[i],"for_plot",sep="_")) ,
                                                    col= c("taxa"), 
                                                    into = c("Phylum", "shorten_taxa"), 
                                                    sep="\\.", extra = "merge",remove = FALSE))
}

for(i in 1:length(list)){
  assign(paste(list[i],"for_plot",sep="_"),get(paste(list[i],"for_plot",sep="_"))[,which(names(get(paste(list[i],"for_plot",sep="_"))) != "shorten_taxa")])
}




## backgroud color by phylum ##
for(i in 1){
  assign(("background_color_list"),rbind(get(paste(list[i],"for_plot",sep="_")) %>% filter(rank=="Phylum") %>% select(taxa)))
}
background_color_list<-unique(background_color_list)
#yarrr_mix<-usecol(c(piratepal("nemo"),piratepal("bugs")))
#color assignment
hcl.pals()
hcl<-hcl.colors(n=nrow(background_color_list), palette ="Dynamic")
#brew<-brewer.pal(n=nrow(background_color_list), name="Spectral")
background_color_list<-cbind(background_color_list, hcl)
colnames(background_color_list)<-c("Phylum","color")
background_color_list<-data.table(Phylum=c("Firmicutes","Bacteroidetes","Proteobacteria","Verrucomicrobia","Actinobacteria"),color=c("#FFFFB3","#8DD3C7","#BEBADA","#80B1D3","#FB8072"))
#background_color_list<-merge(background_color_list, legend_color_list, by="Phylum")
#background_color_list<-background_color_list[,-2]


for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_background_color",sep="_"), merge(get(paste(list[i],"for_plot",sep="_")), background_color_list, all.x = TRUE, by= "Phylum"))
}

for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_background_color",sep="_"), get(paste(list[i],"for_plot_background_color",sep="_")) %>% mutate(`for_plot_background_color` = paste(get(paste(list[i],"for_plot_background_color",sep="_"))$taxa, "annotation_background_color",get(paste(list[i],"for_plot_background_color",sep="_"))$`color`, sep="\t")))
}


## marker color by phylum ##
for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_marker_color",sep="_"), merge(get(paste(list[i],"for_plot",sep="_")), background_color_list, all.x = TRUE, by= "Phylum"))
}
for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_marker_color",sep="_"), get(paste(list[i],"for_plot_marker_color",sep="_")) %>% mutate(`for_plot_marker_color` = paste(get(paste(list[i],"for_plot_marker_color",sep="_"))$taxa, "clade_marker_color",get(paste(list[i],"for_plot_marker_color",sep="_"))$`color`, sep="\t")))
}

##ring_label
for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_ring_label",sep="_"),
         rbind(
           get(paste(list[i],"for_plot",sep="_")) %>% filter(rank=="Species")
         ))}

for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_ring_label",sep="_"), merge(get(paste(list[i],"for_plot_ring_label",sep="_")), ring_label, all.x = TRUE, by= "taxa"))
}


for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_ring_label",sep="_"), get(paste(list[i],"for_plot_ring_label",sep="_")) %>% mutate(`for_plot_ring_label` = paste(get(paste(list[i],"for_plot_ring_label",sep="_"))$taxa, "ring_label",get(paste(list[i],"for_plot_ring_label",sep="_"))$`ring_label`, sep="\t")))
}


## basic option assignment ##
for(i in 1:length(list)){
  assign(paste(list[i],"annot_basic",sep="_"), data.table("for_plot"=c(paste("title",list[i],sep="\t"),
                                                                       "title_font_size	13",
                                                                       "annotation_background_alpha	0.15",
                                                                       "clade_separation	0.35",
                                                                       "annotation_font_stretch	0",
                                                                       "branch_bracket_depth	0.5",
                                                                       "branch_thickness	1.2",
                                                                       "internal_label\t1\tPh.",
                                                                       "internal_label\t2\tClasses",
                                                                       "internal_label\t3\tOrders",
                                                                       "internal_label\t4\tFamilies",
                                                                       "internal_label\t5\tGenera",
                                                                       "internal_label\t6\tSpecies",
                                                                       "internal_labels_rotation\t270",
                                                                       "total_plotted_degrees\t330",
                                                                       "class_legend_font_size\t11",
                                                                       "start_rotation\t270")))}
for(i in 1:length(list)){  
  assign(paste(list[i],"taxa","annot",sep="_"),rbind((setNames(as.data.frame(get(paste(list[i],"for_plot",sep="_"))$`for_plot`),"for_plot")), 
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_"))$`for_plot_annotation`),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_background_color",sep="_"))$`for_plot_background_color`),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_marker_color",sep="_"))$`for_plot_marker_color`),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_ring_label",sep="_"))$`for_plot_ring_label`),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"annot_basic",sep="_"))$`for_plot`),"for_plot"))))
}


for(i in 1:length(list)){  
  assign(paste(list[i],"taxa","annot",sep="_"), get(paste(list[i],"taxa","annot",sep="_"))%>% filter(!grepl("\tNA$",for_plot)))
}


for(i in 1:length(list)){
  write.table(get(paste(list[i],"taxa","annot",sep="_"))$`for_plot`,paste(list[i],"annot.txt", sep="_"), row.names = FALSE, col.names = FALSE, sep = "\n",quote = FALSE)
}

rm(list=ls(pattern=".{1,}_annot_basic$"))
## combining annot db


## data export taxa for plotting ##
for(i in 1:length(list)){
  write.table(get(paste(list[i],sep="_"))$`taxa`,paste(list[i],"taxa.txt", sep="_"), row.names = FALSE, col.names = FALSE, sep = "\n",quote = FALSE)
}































































for(i in 1:length(list)){
  assign(paste(list[i],"taxa","annot",sep="_"),merge(get(paste(list[i],"taxa",sep="_")),ezcloud_taxa_DB,by="taxa"))
}

for(i in 1:length(list)){
  assign(paste(list[i],"taxa","annot",sep="_"),get(paste(list[i],"taxa","annot",sep="_"))[,c(4,1,2,3)] %>% mutate("taxa"=paste(supertaxa,taxa,sep=".")))
}

for(i in 1:length(list)){
  assign(paste(list[i],"taxa","annot",sep="_"),get(paste(list[i],"taxa","annot",sep="_"))[,-1])
}


for(i in 1:length(list)){
  tmp<-get(paste(list[i],"taxa","annot",sep="_"))
  tmp$taxa <- gsub("^NA\\.", "", tmp[,"taxa"])
  assign(paste(list[i],"taxa","annot",sep="_"),tmp)
  rm(tmp)
}

## clade marker size ##

for(i in 1:length(list)){
  assign(paste(list[i],"taxa","sum",sep="_"), get(paste(list[i],"taxa","annot",sep="_")) %>% group_by(variable) %>% dplyr::summarise(`Total No. of taxa`=sum(value)))
}

for(i in 1:length(list)){
  assign(paste(list[i],"taxa","annot",sep="_"), merge(get(paste(list[i],"taxa","annot",sep="_")), get(paste(list[i],"taxa","sum",sep="_")), all.x = TRUE, by= "variable"))
}

for(i in 1:length(list)){
  assign(paste(list[i],"taxa","annot",sep="_"), get(paste(list[i],"taxa","annot",sep="_")) %>% mutate(`clade marker size` = value/`Total No. of taxa`))
}

for(i in 1:length(list)){
  assign(paste(list[i],"taxa","annot",sep="_"), get(paste(list[i],"taxa","annot",sep="_")) %>% mutate(`normalized size` = normalize(`clade marker size`, method = "range", range = c(10, 60), margin = 1, on.constant = "quiet")))
}

for(i in 1:length(list)){
  assign(paste(list[i],"taxa","annot",sep="_"), get(paste(list[i],"taxa","annot",sep="_")) %>% mutate(`for_plot` = paste(get(paste(list[i],"taxa","annot",sep="_"))$taxa, "clade_marker_size",get(paste(list[i],"taxa","annot",sep="_"))$`normalized size`, sep="\t")))
}

## annotation assignment ##
for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_annotation",sep="_"),
         rbind(
           get(paste(list[i],"taxa","annot",sep="_")) %>% filter(variable=="No. of Phylum") %>% mutate(annotation="*"),
           get(paste(list[i],"taxa","annot",sep="_")) %>% filter(variable=="No. of Class") %>% mutate(annotation="*"),
           get(paste(list[i],"taxa","annot",sep="_")) %>% filter(variable=="No. of Order") %>% mutate(annotation="*"),
           get(paste(list[i],"taxa","annot",sep="_")) %>% filter(variable=="No. of Family") %>% mutate(annotation="*"),
           get(paste(list[i],"taxa","annot",sep="_")) %>% filter(variable=="No. of Genus") %>% mutate(annotation="*")
         ))}




for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_annotation",sep="_"), get(paste(list[i],"for_plot_annotation",sep="_")) %>% mutate(`for_plot` = paste(get(paste(list[i],"taxa","annot",sep="_"))$taxa, "annotation",get(paste(list[i],"for_plot_annotation",sep="_"))$`annotation`, sep="\t")))
}

for(i in 1:length(list)){
  assign(paste(list[i],"annot_basic",sep="_"), data.table("for_plot"=c(paste("title",list[i],sep="\t"),
                                                                       "title_font_size	13",
                                                                       "start_rotation	70",
                                                                       "annotation_background_alpha	0.15",
                                                                       "clade_separation	0.35",
                                                                       "annotation_font_stretch	0",
                                                                       "branch_bracket_depth	0.5",
                                                                       "branch_thickness	1.2")))}
for(i in 1:length(list)){  
  assign(paste(list[i],"taxa","annot",sep="_"),rbind((setNames(as.data.frame(get(paste(list[i],"taxa","annot",sep="_"))$`for_plot`),"for_plot")), 
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_"))$`for_plot`),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"annot_basic",sep="_"))$`for_plot`),"for_plot"))))
}

for(i in 1:length(list)){
  write.table(get(paste(list[i],"taxa","annot",sep="_"))$`for_plot`,paste(list[i],"annot.txt", sep="_"), row.names = FALSE, col.names = FALSE, sep = "\n",quote = FALSE)
}
rm(list=ls(pattern=".{1,}+annot_basic$"))
## combining annot db


## data export taxa for plotting ##
for(i in 1:length(list)){
  write.table(get(paste(list[i],"taxa","annot",sep="_"))$`taxa`,paste(list[i],"taxa.txt", sep="_"), row.names = FALSE, col.names = FALSE, sep = "\n",quote = FALSE)
}

## data export annot for plotting ##
for(i in 1:length(list)){
  write.table(get(paste(list[i],"taxa","annot",sep="_"))$`for_plot`,paste(list[i],"annot.txt", sep="_"), row.names = FALSE, col.names = FALSE, sep = "\n",quote = FALSE),
  rm(annot_basic)
}




