pacman::p_load(ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr,ggtree, magrittr,sjmisc, ggseqlogo, gridExtra, ggtree, ggtreeExtra, ggstar, ggnewscale, 
               xlsx, openxlsx, seqinr, ape, rentrez, reutils, colorspace, scales, taxonomizr, phyloR, XML, xml2, reshape2, gtools, data.table, unikn, nord, purrr,
               tidyverse, showtext, extrafont)

list_sample <- ls()

for(i in 1:length(list_sample)){
  enz<-get(list_sample[i])
  names(enz)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore",
                "annotation","species","superkingdom","kingdom","phylum","class","order","family","genus","species",
                "duplicate","unique_read","duplicated_read","total_read")
  assign(list_sample[i], enz)
  enz$duplicate<-gsub("중복", "duplicate", enz$duplicate)
  enz$duplicate<-gsub("??????", "unique", enz$duplicate)
  assign(paste0(list_sample[i],"_filtered"), enz[grep("unique",enz$duplicate),])
  rm(enz)
}

#1.    qseqid    query (e.g., unknown gene) sequence id
#2.    sseqid    subject (e.g., reference genome) sequence id
#3.    pident    percentage of identical matches
#4.    length    alignment length (sequence overlap)
#5.    mismatch    number of mismatches
#6.    gapopen    number of gap openings
#7.    qstart    start of alignment in query
#8.    qend    end of alignment in query
#9.    sstart    start of alignment in subject
#10.    send    end of alignment in subject
#11.    evalue    expect value
#12.    bitscore    bit score

### attatched each environmental info.###
araA_Halo_filtered<-cbind(araA_Halo_filtered, environment ="Halo")
araA_HFD_filtered<-cbind(araA_HFD_filtered, environment ="HFD")
araA_HFD_M_filtered<-cbind(araA_HFD_M_filtered, environment ="HFD_M")
araA_Ind_filtered<-cbind(araA_Ind_filtered, environment ="Ind")
araA_Upo_filtered<-cbind(araA_Upo_filtered, environment ="Upo")

araA_filtered<-rbind(araA_Halo_filtered,
                     araA_HFD_filtered,
                     araA_HFD_M_filtered,
                     araA_Ind_filtered,
                     araA_Upo_filtered)



### KH data ###
`KH_araA_607`<-read.table("araA_607.txt", header=FALSE, sep= "\n", quote="")
name <- grep("^>", `KH_araA_607`$V1, value = TRUE,)
name <- as.data.frame(name)
seq_full<-grep("^>",`KH_araA_607`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,nrow(`KH_araA_607`))
sink("KH_araA_607_seq.csv");for(i in 1:nrow(`name`)){
  w=paste(`KH_araA_607`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(w,"\n")
}
`KH_araA_607_sequence`<-read.csv("KH_araA_607_seq.csv",header=FALSE,sep="\n", stringsAsFactors = FALSE)
`KH_araA_607`<-cbind(name,`KH_araA_607_sequence`)

rm(list=ls(pattern=".{1,}_sequence$")) 
rm(name)

`KH_araA_607`<- separate(`KH_araA_607`, col=name, into=c("sseqid","species"), sep=" ", fill = "right", extra = "merge")


`KH_araA_607`<-separate(`KH_araA_607`, col=species, into=c("annotation","species"), sep=" \\[", fill = "right", extra = "merge")
#`KH_araA_607`<-separate(`KH_araA_607`, col=species, into=c("annotation","species"),sep="\\[(?=[^\\[]+$)",fill="left", extra="merge")
`KH_araA_607`$species<-paste0("[",`KH_araA_607`$species)
names(`KH_araA_607`)[4]<-"sequence"
KH_araA_607$sseqid<-gsub(">","",KH_araA_607$sseqid)
KH_araA_607<-KH_araA_607[,c(1,4)]

tree_info<-merge(`KH_araA_607`, araA_filtered, by ="sseqid",sort=FALSE, all.x = TRUE)


araA_blastx_fasta<-setNames(as.data.frame(cbind(KH_araA_607$sseqid,KH_araA_607$sequence)),c("sseqid","sequence"))
araA_blastx_fasta$sseqid<-substr(araA_blastx_fasta$sseqid, rev(gregexpr("\\|",araA_blastx_fasta$sseqid)[[1]])[2]+1, rev(gregexpr("\\|",araA_blastx_fasta$sseqid)[[1]])[1]-1)
araA_blastx_fasta$sseqid<-paste0(">",araA_blastx_fasta$sseqid)

write.table(araA_blastx_fasta,"araA_blastx_fasta_607.txt", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)

########################################################################################################################################################################
################################################################  phylogeny  ###########################################################################################
########################################################################################################################################################################
araA_tree<-read.tree("T3.raxml.bestTree")

#tip.label<-as.data.frame(araA_tree$tip.label)
#araA_tree$tip.label<-gsub("^[0-9]{1,}_","",araA_tree$tip.label)
#substr(araA_tree$tip.label, (gregexpr("\\_",araA_tree$tip.label)[[1]])[-1], (gregexpr("\\_",araA_tree$tip.label)[[1]])[-1])<-"."
#araA_tree$tip.label<-gsub("\\_[1]{1}$","\\.1",araA_tree$tip.label)

##melt & cast by envs.
tree_info_env<-melt(tree_info, id=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore",
                                    "annotation","species","superkingdom","kingdom","phylum","class","order","family","genus","species",
                                    "duplicate","unique_read","duplicated_read","total_read","sequence","species.1","environment"))
tree_info_env<-dcast(tree_info_env, formula = sseqid ~ environment, fun.aggregate = length)
tree_info<-merge(`KH_araA_607`, araA_filtered, by ="sseqid",sort=FALSE, all.x = FALSE, all.y=FALSE, no.dups = TRUE)
tree_info = tree_info %>% distinct(sseqid,.keep_all = TRUE)
tree_info<-merge(tree_info, tree_info_env, by.x="sseqid")
rm(tree_info_env)

#delete fasta format
tree_info$tip.label<-substr(tree_info$sseqid, rev(gregexpr("\\|",tree_info$sseqid)[[1]])[2]+1, rev(gregexpr("\\|",tree_info$sseqid)[[1]])[1]-1)

character<-setNames(as.data.frame(araA_tree$tip.label),"tip.label")
character<-merge(character,`tree_info`,by.x ="tip.label",sort=FALSE,all.x = TRUE)
character<-character %>% mutate_if(is.character, as.factor)
character<-character %>% mutate_if(is.integer, as.factor)
character<-character %>% mutate_if(is.numeric, as.factor)
character<-character %>% mutate_if(is.double, as.factor)



#manual color
rm(background_color, background_color_list)
background_color_list<-data.table::data.table(phylum=c("Acidobacteria","Actinobacteria","Aquificae","Armatimonadetes","Bacteroidetes","BRC1",
                                           "Caldiserica","Caldithrix_p","Chlamydiae","Chlorobi","Chloroflexi","Chrysiogenetes","Crenarchaeota","Cyanobacteria",
                                           "Deferribacteres","Deinococcus-Thermus","Dictyoglomi","Elusimicrobia","Euryarchaeota",
                                           "Fibrobacteres","Firmicutes","Fusobacteria","Gemmatimonadetes","Lentisphaerae","Nanoarchaeota","Nitrospirae",
                                           "Planctomycetes","Poribacteria_p","Proteobacteria","Spirochaetes","Synergistetes",
                                           "Tenericutes","Thaumarchaeota","Thermodesulfobacteria","Thermotogae","Verrucomicrobia","NA"),
                                  background_color=c("#ff4947", "#f67e7d", "#f69c7d", "#f6b07d", "#AF967D", "#af7a46",
                                                     "#dda46b", "#f2c9a1", "#f8c779", "#ffc661", "#ffb532", "#f59e03", "#f5b247", "#f8c924", 
                                                     "#f3e75b", "#fff098", "#c6dda6", "#bad780", "#919f7f", 
                                                     "#95b46a", "#A3BE8C", "#80d7c6", "#4fc6d0", "#4faad0", "#5c9ce4", "#337cd6", 
                                                     "#5565ca", "#5E81AC", "#AD8CAE", "#bc97ab", "#bb84a1",
                                                     "#cb72a1", "#cc5293", "#a04876",  "#484860", "#184860", "#c0c0c0"))
#LJY_pick
#background_color=c("#ff4947", "#f67e7d", "#f69c7d", "#f6b07d", "#ffd07f", "#fff098", 
#                   "#f3e75b", "#f8c924", "#ffc661", "#ffb532", "#f59e03", "#f5b247", "#f8c779", "#f2c9a1", 
#                   "#dda46b", "#af7a46", "#919f7f", "#95b46a", "#c6dda6", 
#                   "#ddfdb1", "#bad780", "#80d7c6", "#4fc6d0", "#4faad0", "#5c9ce4", "#337cd6", 
#                   "#5565ca", "#6e55ca", "#713b83", "#a04876", "#cc5293", 
#                   "#cb72a1", "#bb84a1", "#bc97ab", "#dcc8d3", "#c0c0c0","#989898")


background_color <- unikn::newpal(col = as.character(background_color_list$background_color),
                                   names = as.character(background_color_list$phylum))
unikn::seecol(background_color)
ggsave("phylum_color.pdf",width=10,height=5,units="cm",dpi=600,limitsize=FALSE)

background_color_list<-background_color_list %>% mutate_if(is.factor, as.character)
background_color_list<-background_color_list %>% mutate_if(is.character, as.factor)

#https://awesomeopensource.com/project/EmilHvitfeldt/r-color-palettes
#polarnight
background_color_list <- data.table(
  phylum = c(paste("polarnight",1:4, sep="_"),
            paste("snowstorm", 1:3, sep="_"),
            paste("frost", 1:4, sep="_"),
            paste("aurora", 1:5, sep="_"),
            paste("lumina", 1:5, sep="_"),
            paste("mountain_forms", 1:5, sep="_"),
            paste("silver_mine", 1:5, sep="_"),
            paste("lake_superior", 1:6, sep="_"),
            paste("victory_bonds", 1:5, sep="_"),
            paste("halifax_harbor",1:6, sep="_"),
            paste("moose_pond",1:8, sep="_"),
            paste("algoma_forest",1:7, sep="_"),
            paste("rocky_mountain",1:6, sep="_"),
            paste("red_mountain",1:8, sep="_"),
            paste("baie_mouton",1:7, sep="_"),
            paste("afternoon_prarie",1:9, sep="_")),
background_color =  c("#2e3440", "#3b4252", "#434c5e", "#4c566a",
                      "#D8DEE9", "#E5E9F0", "#ECEFF4",
                      "#8FBCBB", "#88C0D0", "#81A1C1", "#5E81AC",
                      "#BF616A", "#D08770", "#EBCB8B", "#A3BE8C", "#B48EAD",
                      "#EDDAEB", "#AD8CAE", "#4F93B8", "#306489", "#222B4C",
                      "#184860", "#486078", "#d8d8d8", "#484860", "#181830",
                      "#4B644B", "#647D4B", "#E1E1E1", "#7D96AF", "#647D96",
                      "#7D4B19", "#C89664", "#C87d4B", "#4B647D", "#324B64", "#19324B",
                      "#AF1900", "#C83200", "#E19600", "#193264", "#001964",
                      "#E1C8AF", "#C8AF96", "#AF967D", "#967D7D", "#644B64", "#4B324b",
                      "#4B3232", "#7D4B32", "#966432", "#AF7D32", "#E19632", "#E1AF4B", "#C8C896", "#4B4B4B",
                      "#4B4B4B", "#967D4B", "#AFAF7D", "#C89632", "#647D64", "#96AFAF", "#7D96AF",
                      "#BEBEBE", "#C8C8C8", "#DCD2C8", "#D2C8C8", "#BEBEC8", "#B4B4BE",
                      "#7D3232", "#7D4B4B", "#7D6464", "#AF967D", "#FAC87D", "#E1AF64", "#C8964B", "#32324B",
                      "#304890", "#7890A8", "#90A8C0", "#A8A8A8", "#C0C0A8", "#6A7E4F", "#304848",
                      "#486090", "#6078A8", "#7890A8", "#90A8C0", "#F0D8C0", "#D6BBCF", "#A8C0C0", "#C0D8D8", "#A8A890"))


background_color <- newpal(col = as.character(background_color_list$background_color),
                           names = as.character(background_color_list$phylum))
seecol(background_color)

#palette
par(mfrow=c(8, 2), lheight = 2, mar=rep(1, 4), adj = 0)
walk(names(nord_palettes), nord_show_palette)

#font
#loadfonts()
fonts()
library(extrafont)
loadfonts(device = "win")
windowsFonts("MyriadPro.ttf")

font_import(pattern = "Myriad.*", prompt=FALSE)
font_add("Myriad Pro",
         regular = "MyriadPro-Regular.ttf",
         bold = "MyriadPro-Bold.ttf",
         italic = "MyriadPro-It.ttf")

#font_add_google("Montserrat", "Montserrat")
font_add('Myriad Pro', file.choose())
showtext_auto()

#Phylogeny
OO<-ggtree(araA_tree, layout="circular", size=0.5, linetype=1)+xlim(-0.5,NA)+#geom_nodelab2(aes(label=araA_tree$node.label),size = 2.5, col= "grey")+
  geom_tiplab2(aes(label=gsub("\\[|\\]","",paste(p$data$species))),color="black", fontface="italic", size=2, align = TRUE, linesize=0, linetype = 0, offset = 1, show.legend = FALSE)+
  geom_tiplab2(aes(color=p3$data$phylum, label=""),fontface="italic", size=2, align = TRUE, linesize=0.1, linetype = 3, offset = 1, show.legend = FALSE)+
  #geom_text(size = 2,align =TRUE, aes(family = "Myriad Pro", angle=angle, color=p3$data$phylum, label=gsub("\\[|\\]","",paste(p$data$species))))+
  scale_color_manual(name="phylum", values=c(background_color), na.value = "grey90") + theme(legend.position="right")+
  geom_treescale(x=5, y=0)
#geom_text(aes(label=node))

#as gradient
p<-OO %<+% character
p2<-p %<+% background_color_list
p3<-p2 + aes(color=phylum, show.legend = FALSE)+
  #scale_color_discrete("Phylum")+
  #$guides(color = guide_legend(override.aes = list(size = 6),shape = 21))+
  new_scale("color")+
  
  geom_tippoint(size=1.5, shape=16, aes(color=Halo, x=4.2, na.value = NA)) +
  scale_color_manual(values = c("white","#D2A54A"), na.value = NA)+
  new_scale("fill") +
  new_scale("color")+
  
  geom_tippoint(size=1.5, shape=16, aes(color=HFD, x=4.3, na.value = NA)) +
  scale_color_manual(values = c("white","#7F7F7F"), na.value = NA)+
  new_scale("fill") +
  new_scale("color")+
  
  geom_tippoint(size=1.5, shape=16, aes(color=HFD_M, x=4.4, na.value = NA))+
  scale_color_manual(values = c("white","#262626"), na.value = NA)+
  new_scale("fill") +
  new_scale("color")+
  
  geom_tippoint(size=1.5, shape=16, aes(color=Ind, x=4.5, na.value = NA))+
  scale_color_manual(values = c("white","#C00000"), na.value = NA)+
  new_scale("fill") +
  new_scale("color")+
  
  geom_tippoint(size=1.5, shape=16, aes(color=Upo, x=4.6, na.value = NA))+
  scale_color_manual(values = c("white","#548235"), na.value = NA)+
  new_scale("fill")+
  new_scale("color")

p3

#showtext_end()

#tidytree::MRCA(p, tip=c("WP_014026435.1","WP_023790968.1"))
#p2<-p %>% as.treedata %>% as_data_frame
#parent(p2, 222)
#child(p2,2130)
#viewClade(p2, node=2130)
#p3<-p3 %>% collapse(node=445)
#p3

#C00000 (Ind,hot spring)
#D2A54A (Halo)
#7F7F7F (HFD)
#262626 (HFD-M)
#548235 (Upo)
#916B45 (Soil)
#7CAFDE (SB)
#3B87CD (NewR)

ggsave("araA_tree_blastx.pdf",device = "pdf",width=60,height=60,units="cm",dpi=600, limitsize=FALSE, useDingbats=FALSE)

###############
background_color_list<-data.table(phylum=c("Acidobacteria","Actinobacteria","Apicomplexa","Armatimonadetes","Ascomycota",
                                           "Bacteroidetes","Balneolaeota", "Basidiomycota","Chlamydiae","Chlorobi",
                                           "Chloroflexi", "Chordata", "Cyanobacteria", "Deinococcus-Thermus", "Euglenozoa",
                                           "Euryarchaeota", "Firmicutes", "Lentisphaerae", "Planctomycetes", "Proteobacteria",
                                           "Spirochaetes", "Streptophyta", "Tenericutes", "Thermotogae", "Verrucomicrobia",
                                           "NA"),
                                  background_color=c("#A71B4B", "#B42B49", "#D64612", "#CE473B", "#DA542D", 
                                                     "#EBA47A", "#FEC19E", "#ECD1C0", "#F2DFA4", "#F6BA62", 
                                                     "#DCA958", "#E69C30", "#FFDC00", "#FBD274", "#41AA4A",
                                                     "#607C3C", "#005B00", "#98A5C0", "#7688BB", "#0071BB",
                                                     "#304FAF", "#896378", "#6E4468", "#4B2A53", "#0D2939",
                                                     "#d3d3d3"))



background_color_list<-data.table(phylum=c("Acidobacteria", "Actinobacteria", "Armatimonadetes", "Bacteroidetes", "Balneolaeota",
                                           "Chlorobi", "Chloroflexi", "Deinococcus-Thermus", "Firmicutes", "Lentisphaerae",
                                           "Planctomycetes", "Proteobacteria", "Spirochaetes", "Thermotogae",  "Verrucomicrobia"),
                                  background_color=c("#771218", "#C00001", "#D64612", "#E25B6A", "#F08B83"
                                                     , "#EBA47A", "#FEC19E", "#ECD1C0", "#CCCCCC", "#C0B2B2",
                                                     "#A5878F", "#896378", "#6E4468", "#4B2A53", "#0D2939"))
