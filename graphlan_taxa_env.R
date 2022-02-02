#data from blastx tree
pacman::p_load(ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr,ggtree, magrittr, sjmisc, ggseqlogo, gridExtra, ggtree, ggtreeExtra, ggstar, seqinr, ape, rentrez, reutils, colorspace, scales, taxonomizr, XML, xml2, reshape2, gtools, BBmisc, unikn, data.table)

#scheme
#1. melting id by sseqid -> con.aggregate=length (clustering reads annotated same protein)
#2. filtering sseqid by KH's filtered reads

katE_taxa_count<-character[,1:28]


### attatched each environmental info.###
katE_Halo_count<-cbind(katE_Halo, environment ="Halo")

##melt & cast by envs.
katE_Halo_count<-melt(katE_Halo_count, id=c("sseqid","environment"))
katE_Halo_count<-dcast(katE_Halo_count, formula = sseqid ~ environment, fun.aggregate = length)
katE_Halo_count$Halo<-katE_Halo_count$Halo/(ncol(katE_Halo)-2)
katE_Halo_count_sum<-sum(katE_Halo_count$Halo)
katE_Halo_count<-merge(katE_Halo_count, katE_taxa_count, by.x="sseqid")

##melt & cast by envs. with reads count
katE_Halo_count2<-melt(katE_Halo_count, measure.vars="Halo")
#phylum
katE_Halo_count_phylum<-dcast(katE_Halo_count2, formula = phylum ~ variable, fun.aggregate = sum)
katE_Halo_count_phylum<-cbind(katE_Halo_count_phylum, Rank="Phylum")
katE_Halo_count_phylum<-cbind(katE_Halo_count_phylum,taxa=paste(katE_Halo_count_phylum$phylum,sep="."))
names(katE_Halo_count_phylum)[names(katE_Halo_count_phylum) == 'Halo'] <- 'value'
#class
katE_Halo_count_class<-dcast(katE_Halo_count2, formula = phylum + class ~ variable, fun.aggregate = sum)
katE_Halo_count_class<-cbind(katE_Halo_count_class, Rank="Class")
katE_Halo_count_class<-cbind(katE_Halo_count_class,taxa=paste(katE_Halo_count_class$phylum,katE_Halo_count_class$class,sep="."))
names(katE_Halo_count_class)[names(katE_Halo_count_class) == 'Halo'] <- 'value'
#order
katE_Halo_count_order<-dcast(katE_Halo_count2, formula = phylum + class + order ~ variable, fun.aggregate = sum)
katE_Halo_count_order<-cbind(katE_Halo_count_order, Rank="Order")
katE_Halo_count_order<-cbind(katE_Halo_count_order,taxa=paste(katE_Halo_count_order$phylum,katE_Halo_count_order$class,katE_Halo_count_order$order,sep="."))
names(katE_Halo_count_order)[names(katE_Halo_count_order) == 'Halo'] <- 'value'
#family
katE_Halo_count_family<-dcast(katE_Halo_count2, formula = phylum + class + order + family ~ variable, fun.aggregate = sum)
katE_Halo_count_family<-cbind(katE_Halo_count_family, Rank="Family")
katE_Halo_count_family<-cbind(katE_Halo_count_family,taxa=paste(katE_Halo_count_family$phylum,katE_Halo_count_family$class,katE_Halo_count_family$order,katE_Halo_count_family$family,sep="."))
names(katE_Halo_count_family)[names(katE_Halo_count_family) == 'Halo'] <- 'value'
#genus
katE_Halo_count_genus<-dcast(katE_Halo_count2, formula = phylum + class + order + family + genus ~ variable, fun.aggregate = sum)
katE_Halo_count_genus<-cbind(katE_Halo_count_genus, Rank="Genus")
katE_Halo_count_genus<-cbind(katE_Halo_count_genus,taxa=paste(katE_Halo_count_genus$phylum,katE_Halo_count_genus$class,katE_Halo_count_genus$order,katE_Halo_count_genus$family,katE_Halo_count_genus$genus,sep="."))
names(katE_Halo_count_genus)[names(katE_Halo_count_genus) == 'Halo'] <- 'value'


taxa<-c("Phylum", "Class", "Order", "Family", "Genus")
list<-c("phylum","class","order","family","genus")

## clade marker size ##
for(i in 1:length(list)){
  assign(paste("katE","Halo","count",list[i],"sum",sep="_"), get(paste("katE","Halo","count",list[i],sep="_")) %>% group_by(Rank) %>% dplyr::summarise(`Total No. of taxa`=sum(value)))
}

for(i in 1:length(list)){
  assign(paste("katE","Halo","count",list[i],sep="_"), merge(get(paste("katE","Halo","count",list[i],sep="_")), get(paste("katE","Halo","count",list[i],"sum",sep="_")), all.x = TRUE, by= "Rank"))
}

rm(list=ls(pattern=".{1,}_sum$"))


for(i in 1:length(list)){
  assign(paste("katE","Halo","count",list[i],sep="_"), get(paste("katE","Halo","count",list[i],sep="_")) %>% mutate(`clade marker size` = value/`Total No. of taxa`))
}

for(i in 1:length(list)){
  assign(paste("katE","Halo","count",list[i],sep="_"), get(paste("katE","Halo","count",list[i],sep="_")) %>% mutate(`normalized size` = normalize(`clade marker size`, method = "range", range = c(5, 60), margin = 1, on.constant = "quiet")))
}

for(i in 1:length(list)){
  assign(paste("katE","Halo","count",list[i],sep="_"), get(paste("katE","Halo","count",list[i],sep="_")) %>% mutate(`for_plot` = paste(get(paste("katE","Halo","count",list[i],sep="_"))$taxa, "clade_marker_size",get(paste("katE","Halo","count",list[i],sep="_"))$`normalized size`, sep="\t")))
}

#combine data group by taxonomic rank
katE_Halo_for_plot_annotation<-rbind.fill(katE_Halo_count_phylum,katE_Halo_count_class,katE_Halo_count_order,katE_Halo_count_family,katE_Halo_count_genus)

list<-"katE_Halo"
## backgroud color by phylum ##
for(i in 1){
  assign(("background_color_list"),rbind(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Phylum") %>% select(taxa)))
}
#for(i in 2:length(list)){
#  assign(("background_color_list"),rbind(`background_color_list`,(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Phylum") %>% select(taxa))))
#}

background_color_list<-unique(background_color_list)
#yarrr_mix<-usecol(c(piratepal("nemo"),piratepal("bugs")))
#color assignment
hcl.pals()
hcl<-hcl.colors(n=nrow(background_color_list), palette ="Dynamic")
#brew<-brewer.pal(n=nrow(background_color_list), name="Spectral")
background_color_list<-data.table(phylum=c("Acidobacteria","Actinobacteria","Aquificae","Armatimonadetes","Bacteroidetes","BRC1",
                                           "Caldiserica","Caldithrix_p","Chlamydiae","Chlorobi","Chloroflexi","Chrysiogenetes","Crenarchaeota","Cyanobacteria",
                                           "Deferribacteres","Deinococcus-Thermus","Dictyoglomi","Elusimicrobia","Euryarchaeota",
                                           "Fibrobacteres","Firmicutes","Fusobacteria","Gemmatimonadetes","Lentisphaerae","Nanoarchaeota","Nitrospirae",
                                           "Planctomycetes","Poribacteria_p","Proteobacteria","Spirochaetes","Synergistetes",
                                           "Tenericutes","Thaumarchaeota","Thermodesulfobacteria","Thermotogae","Verrucomicrobia","NA"),
                                  background_color=c("#A71B4B", "#B42B49", "#C13944", "#CE473B", "#DA542D", "#E5610A",
                                                     "#E97202", "#ED820A", "#F1911B", "#F49F2D", "#F6AD3E", "#F8BA50", "#FAC662", "#FBD274", 
                                                     "#FCDE85", "#FDE896", "#FEF2A7", "#FEFAB7", "#F6FCBB", 
                                                     "#E4F8B5", "#D0F4B1", "#BAEEAE", "#A2E7AD", "#8AE0AD", "#6FD8AE", "#52CFB0", 
                                                     "#2EC6B2", "#00BCB4", "#00B1B5", "#00A5B6", "#0099B5", 
                                                     "#008BB4", "#057DB1", "#326EAD", "#485DA7", "#584B9F","#d3d3d3"))
background_color_list<-cbind(background_color_list,abbreviation=1:nrow(background_color_list))




#merge coloring
for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_annotation",sep="_"), merge(get(paste(list[i],"for_plot_annotation",sep="_")), background_color_list, all.x = TRUE, by= "phylum"))
}

##delete sum
#rm(list=ls(pattern=".{1,}_class$|.{1,}_order$|.{1,}_family$|.{1,}_genus$"))


## data export taxa for plotting ##
for(i in 1:length(list)){
  write.table(get(paste(list[i],"for_plot_annotation",sep="_"))$`taxa`,paste(list[i],"taxa.txt", sep="_"), row.names = FALSE, col.names = FALSE, sep = "\n",quote = FALSE)
}


#------------------------------------------OK--------------------------




#for(i in 1:length(list)){
#  assign(paste(list[i],"for_plot_annotation",sep="_"), get(paste(list[i],"for_plot_annotation",sep="_")) %>% mutate_if(is.factor, as.character))
#}


for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_annotation",sep="_"), get(paste(list[i],"for_plot_annotation",sep="_")) %>% mutate(`for_plot_background_color` = paste(get(paste(list[i],"for_plot_annotation",sep="_"))$taxa, "annotation_background_color",get(paste(list[i],"for_plot_annotation",sep="_"))$`background_color`, sep="\t")))
}

for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_annotation",sep="_"), get(paste(list[i],"for_plot_annotation",sep="_")) %>% mutate(`for_plot_annotation` = paste(get(paste(list[i],"for_plot_annotation",sep="_"))$taxa, "annotation",get(paste(list[i],"for_plot_annotation",sep="_"))$`annotation`, sep="\t")))
}

for(i in 1:length(list)){
  assign(paste(list[i],"for_plot_annotation",sep="_"), get(paste(list[i],"for_plot_annotation",sep="_")) %>% mutate(`for_plot_marker_color` = paste(get(paste(list[i],"for_plot_annotation",sep="_"))$taxa, "clade_marker_color",get(paste(list[i],"for_plot_annotation",sep="_"))$`background_color`, sep="\t")))
}


## basic option assignment ##
for(i in 1:length(list)){
  assign(paste(list[i],"annot_basic",sep="_"), data.table("for_plot"=c(paste("title",list[i],sep="\t"),
                                                                       "title_font_size	13",
                                                                       "annotation_background_alpha	0.15",
                                                                       "clade_separation	0.35",
                                                                       "annotation_font_stretch	0",
                                                                       "branch_bracket_depth	0.5",
                                                                       "branch_thickness	1.0",
                                                                       "internal_label\t1\tPh.",
                                                                       "internal_label\t2\tClasses",
                                                                       "internal_label\t3\tOrders",
                                                                       "internal_label\t4\tFamilies",
                                                                       "internal_label\t5\tGenera",
                                                                       "internal_labels_rotation\t270",
                                                                       "total_plotted_degrees\t330",
                                                                       "class_legend_font_size\t11",
                                                                       "start_rotation\t270")))}
for(i in 1:length(list)){
  assign(paste(list[i],"taxa","annot",sep="_"),rbind((setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_"))$`for_plot`),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Phylum") %>% select(c(`for_plot_annotation`,`abbreviation`)) %>% unite("for_plot_annotation",c(`for_plot_annotation`,`abbreviation`),sep="") %>% select(`for_plot_annotation`)),"for_plot")), 
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Phylum") %>% select(`for_plot_background_color`)),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Class") %>% select(`for_plot_background_color`)),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Phylum") %>% select(`for_plot_marker_color`)),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Class") %>% select(`for_plot_marker_color`)),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Order") %>% select(`for_plot_marker_color`)),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Family") %>% select(`for_plot_marker_color`)),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Genus") %>% select(`for_plot_marker_color`)),"for_plot")),
                                                     (setNames(as.data.frame(get(paste(list[i],"annot_basic",sep="_"))$`for_plot`),"for_plot"))))
}
#(setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Order") %>% select(`for_plot_background_color`)),"for_plot")),
#(setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Family") %>% select(`for_plot_background_color`)),"for_plot")),
#(setNames(as.data.frame(get(paste(list[i],"for_plot_annotation",sep="_")) %>% filter(Rank=="Genus") %>% select(`for_plot_background_color`)),"for_plot")),


for(i in 1:length(list)){  
  assign(paste(list[i],"taxa","annot",sep="_"), get(paste(list[i],"taxa","annot",sep="_")) %>% filter(!grepl("\tNA$|NA\tannotation\t*",for_plot)))
}


for(i in 1:length(list)){
  write.table(get(paste(list[i],"taxa","annot",sep="_"))$`for_plot`,paste(list[i],"annot.txt", sep="_"), row.names = FALSE, col.names = FALSE, sep = "\n",quote = FALSE)
}

rm(list=ls(pattern=".{1,}_annot_basic$"))
## combining annot db



#write shell script
for(i in 1:length(list)){
  write.table(
    paste((paste("graphlan_annotate.py", paste(list[i],"taxa.txt", sep="_"), paste(list[i], "annot.xml", sep="_"), "--annot", paste(list[i], "annot.txt", sep="_"),sep=" ")),
          (paste("graphlan.py", paste(list[i], "annot.xml", sep="_"), paste(list[i], "annot.svg", sep="_"), "--dpi 300 --size 4.5 --pad 0", sep=" ")),sep="\n"), 
    paste0(list[i],"_run.sh"), quote = FALSE, row.names = FALSE, col.names = FALSE
  )}






















#filling missing taxa
#Taxonomy data
{
  for(i in 1){
    protein_search<- entrez_search(db="protein", term=katE_Halo_count$tip.label[i],retmode = "xml", retmax = 999999, use_history = TRUE)
    taxonomy_data <- entrez_link(dbfrom = "protein", id = protein_search$ids, db ="taxonomy")
    taxonomy_id <- xmlToDataFrame(entrez_fetch(db="taxonomy", id=taxonomy_data$links$protein_taxonomy[1], rettype="xml")) %>% select("TaxId")
    raw_recs <- taxize::classification(taxonomy_id[[1]], db = 'ncbi')
    #raw_recs <- as.data.frame.list(raw_recs)
    raw_recs <- do.call(rbind.data.frame, raw_recs)
    rownames(raw_recs)<-NULL
    raw_recs<-melt(raw_recs, id=c("rank","id"),variable.name="name")
    raw_recs<-dcast(raw_recs, formula = name ~ rank, fun.aggregate = function(x){paste(x,collapse = "_")})
    raw_recs<-cbind(tip.label=katE_Halo_count$tip.label[i], protein_id=protein_search$ids[1], tax_id=taxonomy_data$links$protein_taxonomy[1], raw_recs)
    taxa_assignment<-raw_recs
  }
  
  for(i in (1:nrow(katE_Halo_count))[-1]){
    protein_search<- entrez_search(db="protein", term=katE_Halo_count$tip.label[i], retmode = "xml", retmax = 999999, use_history = TRUE)
    taxonomy_data <- entrez_link(dbfrom = "protein", id = protein_search$ids, db ="taxonomy")
    if(!is.null(taxonomy_data$links$protein_taxonomy[1])){
      taxonomy_id <- xmlToDataFrame(entrez_fetch(db="taxonomy", id=taxonomy_data$links$protein_taxonomy[1], rettype="xml")) %>% select("TaxId")
      raw_recs <- taxize::classification(taxonomy_id[[1]], db = 'ncbi')
      #raw_recs <- as.data.frame.list(raw_recs)
      raw_recs <- do.call(rbind.data.frame, raw_recs)
      rownames(raw_recs)<-NULL
      raw_recs<-melt(raw_recs, id=c("rank","id"),variable.name="name")
      raw_recs<-dcast(raw_recs, formula = name ~ rank, fun.aggregate = function(x){paste(x,collapse = "_")})
      raw_recs<-cbind(tip.label=katE_Halo_count$tip.label[i], protein_id=protein_search$ids[1], tax_id=taxonomy_data$links$protein_taxonomy[1], raw_recs)
      #raw_recs <- data.frame(matrix(unlist(raw_recs), nrow=length(raw_recs), byrow=T))
      #filter(rank=c("phylum","class","order","family","genus","species","strain"))
    }
    else{raw_recs<-NULL}
    taxa_assignment<-rbind.fill(taxa_assignment,raw_recs)} 
}




taxa_assignment<-taxa_assignment[c("tip.label","protein_id","tax_id","superkingdom","phylum","class","order","family","genus","species","strain")]
names(taxa_assignment)<-paste("ncbi",c("protein_accession","protein_id","tax_id","superkingdom","phylum","class","order","family","genus","species","strain"),sep="_")
taxa_assignment$tip.label<-gsub(">","",taxa_assignment$tip.label)
names(taxa_assignment)[1]<-"tip.label"

katE_Halo_count_with_ncbi<-merge(katE_Halo_count, taxa_assignment, by ="tip.label")
##########################################################################################################################################################################




##melt & cast by envs.
tree_info<-merge(`KH_katE_2109`, katE_filtered, by ="sseqid",sort=FALSE, all.x = FALSE, all.y=FALSE, no.dups = TRUE)
tree_info = tree_info %>% distinct(sseqid,.keep_all = TRUE)
tree_info<-merge(tree_info, tree_info_env, by.x="sseqid")
rm(tree_info_env)







