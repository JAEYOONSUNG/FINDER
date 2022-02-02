pacman::p_load(ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr,ggtree, magrittr,sjmisc, ggseqlogo, gridExtra, ggtree, ggtreeExtra, ggstar, ggnewscale, 
               xlsx, openxlsx, seqinr, ape, rentrez, reutils, colorspace, scales, taxonomizr, phyloR, XML, xml2, reshape2, gtools, data.table, unikn, nord, purrr,
               tidyverse, showtext, extrafont, myriad, Rttf2pt1)
devtools::install_github("kjhealy/myriad")
devtools::install_github("wch/Rttf2pt1")
##synonym



#ytoP / glutamyl aminopeptidase /3.4.21.62[EC/RN Number=] 
#{
#  `ytoP`<- entrez_search(db="protein", term="glutamyl aminopeptidase",retmode = "xml", retmax = 999999, use_history = TRUE)
#  for(seq_start in seq(1,length(`ytoP`$ids),500)){
#    recs <- entrez_fetch(db="protein", web_history=`ytoP`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
#    cat(recs, file="ytoP_135148.fasta", append=TRUE)
#    cat(seq_start+499, "sequences downloaded\r")
#  }}

{
  `ytoP`<- entrez_search(db="protein", term="M42 family metallopeptidase",retmode = "xml", retmax = 999999, use_history = TRUE)
  for(seq_start in seq(1,length(`ytoP`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`ytoP`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="ytoP_M42_23811.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}

{
  `ytoP`<- entrez_search(db="protein", term="3.4.11.7[EC]",retmode = "xml", retmax = 999999, use_history = TRUE)
  for(seq_start in seq(1,length(`ytoP`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`ytoP`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="ytoP_EC_54959.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}


#data partitioning
`ytoP`<-read.table("ytoP_EC_54981.fasta", header=FALSE, sep= "\t", quote="")
name <- grep("^>", `ytoP`$V1, value = TRUE,)
name <- as.data.frame(name)
seq_full<-grep("^>",`ytoP`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,nrow(`ytoP`))
sink("ytoP_seq.csv");for(i in 1:nrow(`name`)){
  w=paste(`ytoP`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(w,"\n")
}
`ytoP_sequence`<-read.csv("ytoP_seq.csv",header=FALSE,sep="\n")
`ytoP`<-cbind(name,`ytoP_sequence`)

rm(list=ls(pattern=".{1,}_sequence$")) 
rm(name)

`ytoP`<- separate(`ytoP`, col=name, into=c("accession","species"), sep=" ", fill = "right", extra = "merge")

##filter Swiss-Prot/UniProtKB
`ytoP_filtered`<-separate(`ytoP`, col=species, into=c("annotation","species"), sep=" \\[", fill = "right", extra = "merge")
#`ytoP_filtered`<-separate(`ytoP_filtered`, col=species, into=c("annotation","species"),sep="\\[(?=[^\\[]+$)",fill="left", extra="merge")
`ytoP_filtered`$species<-paste0("[",`ytoP_filtered`$species)
names(`ytoP_filtered`)[4]<-"sequence"
ytoP_filtered$sequence<-as.character(ytoP_filtered$sequence)
`ytoP_filtered`<-mutate(`ytoP_filtered`, number =nchar(`ytoP_filtered`$sequence))


#template sequence
template_sequence<-filter(`ytoP_filtered`,grepl("SPY14977.1|ACT28353.1|VHM42516.1|KYD34354.1|ACV58083.1|WP_106450884.1|OLS31414.1|SJN29768.1|CUQ39877.1|WP_092092338.1",accession))

#CAB14964.1=WP_003229261.1	Bacillus subtilis subsp. subtilis str. 168 (putative modified amino acid aminopeptidase) (SPY14977.1)
#YP_006171000.1=WP_001019482.1 =Escherichia coli P12b (aminopeptidase) (ACT28353.1 [Escherichia coli 'BL21-Gold(DE3)pLysS AG'])
#VHM42516.1 [Streptococcus pyogenes] (newly)
#GAD11976.1=Geobacillus kaustophilus GBlys (KYD34354.1 [Geobacillus stearothermophilus])
#ACV58083.1 [Alicyclobacillus acidocaldarius subsp. acidocaldarius DSM 446]
#WP_106450884.1 [Trichococcus alkaliphilus]
#OLS31414.1 [Candidatus Heimdallarchaeota archaeon AB_125]
#SJN29768.1 [Marinilactibacillus psychrotolerans 42ea]
#CUQ39877.1 [[Ruminococcus] torques]
#WP_092092338.1 [Pisciglobus halotolerans]
#EAL89884.1 [Aspergillus fumigatus Af293] (skip)


#AFL66590.1=WP_014767491.1	Desulfurococcus fermentans DSM 16532 (delete)
#YP_005688133.1=WP_003516515.1 Ruminiclostridium thermocellum DSM 1313(Clostridium thermocellum DSM 1313) M42 family metallopeptidase (delete)
#YP_216595.1=WP_000144624.1=Salmonella enterica subsp. enterica serovar Choleraesuis str. SC-B67 (cellulase) M42 family metallopeptidase (delete)
#YP_001739773.1=WP_012311325.1 =Thermotoga sp. RQ2 (AEH50186.1)
#NP_418725.4=WP_000010829.1 Escherichia coli str. K-12 substr. MG1655(delete)
#ADB57229.1=WP_012939565.1 Archaeoglobus profundus DSM 5631(M42 family metallopeptidase) (delete)


`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",accession))
`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]{1,2}",accession))
`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("pdb\\|[0-9][A-Z0-9]{3}",accession))
`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("patent",species))
##RefSeq accession numbers and molecule types.
#AP_	Protein	Annotated on AC_ alternate assembly
#NP_	Protein	Associated with an NM_ or NC_ accession
#YP_c	Protein	Annotated on genomic molecules without an instantiated transcript record
#XP_c	Protein	Predicted model, associated with an XM_ accession
#WP_	Protein	Non-redundant across multiple strains and species

##filter empty species info; pdb, US patent
`ytoP_filtered`<-filter(`ytoP_filtered`, !is.na(species))
`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("hypothetical|unnamed|partial|predicted|putative|LOW QUALITY PROTEIN|Uncharacterised|Uncharacterized|PREDICTED|UNVERIFIED|fragment|precursor|inhibitor|patent|probable|possible",annotation,ignore.case=TRUE))
`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("arginase|fructose",annotation,ignore.case=TRUE))
`ytoP_filtered`<-filter(`ytoP_filtered`, grepl("",annotation,ignore.case=TRUE))
summary<-summary(`ytoP_filtered`$number)

annotation_list<-as.data.frame(unique(ytoP_filtered$annotation))


##sequence length
#`ytoP_filtered`<-filter(`ytoP_filtered`, number >358, ignore.case=TRUE)
`ytoP_filtered`<-filter(`ytoP_filtered`, number >300 & number<400 ,ignore.case=TRUE)


##PRATT motif (expasy prosite; PDOC00129: Neutral zinc metallopeptidases, zinc-binding region signature)
#ZINC_PROTEASE, PS00142; Neutral zinc metallopeptidases, zinc-binding region signature  (PATTERN)
#[GSTALIVN]-{PCHR}-{KND}-H-E-[LIVMFYW]-{DEHRKP}-H-{EKPC}-[LIVMFYWGSPQ]
#`ytoP_filtered`<-filter(`ytoP_filtered`, grepl("([GSTALIVN][^PCHR][^KND]HE[LIVMFYW][^DEHRKP]H[^EKPC][LIVMFYWGSPQ])",sequence))

##delete duplicated sequence
`ytoP_filtered`<-cbind(`ytoP_filtered`, "duplicate"=duplicated(`ytoP_filtered`$sequence))
`ytoP_filtered`<-subset(`ytoP_filtered`, `ytoP_filtered`$duplicate==FALSE) 





#re-confirm contaminant sequences
annotation_list<-as.data.frame(unique(ytoP_filtered$annotation))

##manually filterate ytoP gene (if needed)
#`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("autophagic|",annotation, ignore.case = TRUE))
#`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("catalase/peroxidase",annotation, ignore.case = TRUE))
#`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("SpoS|TPA|similar|spore|UNVERIFIED|probable|domain",annotation, ignore.case = TRUE))
#`ytoP_filtered`<-filter(`ytoP_filtered`, !grepl("(SWPDN).*(VQMGLIYV).*",sequence))





fasta<-as.data.frame(paste(`ytoP_filtered`$accession))
fasta<-setNames(cbind(fasta,`ytoP_filtered`$sequence),c("name","sequence"))

write.table(`fasta`,"ytoP_EC_2935_fasta.txt", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)



########check identity and tree distribution#######
#crawling fasta file
src_dir<-getwd()
src_file <- list.files(src_dir, pattern = "_[0-9]{2}\\_[0-9]{1,6}\\.txt$") # list
src_file_cluster<-as.data.frame(src_file)
src_file_cluster<-substr(src_file_cluster$src_file, rev(gregexpr("\\_",src_file_cluster$src_file)[[1]])[2], rev(gregexpr("\\.",src_file_cluster$src_file)[[1]])[1])
src_file_cluster<-gsub("\\.","",src_file_cluster)

#^(.(?!(some text)))*$
#for(i in 1:length(src_file))
for(i in 1:length(src_file_cluster)){
  assign(paste0("ytoP_clustered",src_file_cluster[i]),read.table(src_file[i], header=FALSE, quote=""))
  assign("name",grep("^>", get(paste0("ytoP_clustered",src_file_cluster[i]))$V1, value = TRUE))
  assign("name",as.data.frame(name))
  assign("seq_full",grep("^>",get(paste0("ytoP_clustered",src_file_cluster[i]))$V1))
  assign("seq_start",seq_full+1)
  assign("seq_finish",c(seq_full[-1]-1,nrow(get(paste0("ytoP_clustered",src_file_cluster[i])))))
  
  sink(paste0("ytoP_clustered",src_file_cluster[i],"_seq.csv"));for(j in 1:nrow(`name`)){
    assign("w",paste(get(paste0("ytoP_clustered",src_file_cluster[i]))$V1[seq_start[j]:seq_finish[j]],collapse = ""))
    cat(w,"\n")
  }
  assign(paste0("ytoP_clustered",src_file_cluster[i]), read.csv(paste0("ytoP_clustered",src_file_cluster[i],"_seq.csv"),header=FALSE,sep="\n"))
  assign(paste0("ytoP_clustered",src_file_cluster[i]), cbind(name,get(paste0("ytoP_clustered",src_file_cluster[i]))))
  assign(paste0("ytoP_clustered",src_file_cluster[i]), setNames(get(paste0("ytoP_clustered",src_file_cluster[i])),c("accession","sequence")))
}


#add reference sequences and delete duplicated reference seq
#template_aligned_sequence<-filter(`fasta`,grepl("AAB47761.1|NP_391784.2|WP_000077872.1|AAC18407.1|WP_014589250.1|CEJ79573.1|TMW79725.1",name))
#add template sequences
template_fasta<-as.data.frame(paste(`template_sequence`$accession))
template_fasta<-setNames(cbind(template_fasta,`template_sequence`$sequence),c("accession","sequence"))

for (i in 1:length(src_file_cluster)){
  assign(paste0("ytoP_clustered",src_file_cluster[i]),rbind(template_fasta, get(paste0("ytoP_clustered",src_file_cluster[i]))))
  assign(paste0("ytoP_clustered",src_file_cluster[i]),cbind(get(paste0("ytoP_clustered",src_file_cluster[i])), "duplicate"=duplicated(get(paste0("ytoP_clustered",src_file_cluster[i]))$sequence)))
  assign(paste0("ytoP_clustered",src_file_cluster[i]),subset(get(paste0("ytoP_clustered",src_file_cluster[i])), get(paste0("ytoP_clustered",src_file_cluster[i]))$duplicate==FALSE))
  assign(paste0("ytoP_clustered",src_file_cluster[i]),get(paste0("ytoP_clustered",src_file_cluster[i]))[,-3])
  assign(paste0("ytoP_clustered",src_file_cluster[i]),merge(`ytoP_filtered`,get(paste0("ytoP_clustered",src_file_cluster[i])),by="accession", sort=FALSE,all.x = FALSE, no.dups = TRUE, how='outer', suffixes=c('', '_y')))
  assign(paste0("ytoP_clustered",src_file_cluster[i]),get(paste0("ytoP_clustered",src_file_cluster[i]))[,-7])
}

for(i in 1:length(src_file_cluster)){
  assign(paste0("ytoP_clustered_fasta",src_file_cluster[i]), cbind(name=get(paste0("ytoP_clustered",src_file_cluster[i]))$accession,get(paste0("ytoP_clustered",src_file_cluster[i]))$sequence))
  write.table(get(paste0("ytoP_clustered_fasta",src_file_cluster[i])),paste0("ytoP_clustered_fasta",src_file_cluster[i]), row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)
}


########################################################################################################################################################################
################################################################   phylogeny  ##########################################################################################
########################################################################################################################################################################
ytoP_tree<-read.tree("ytoP_70_nj.20112913245290.on.nh.txt")
ytoP_tree$tip.label<-gsub("^[0-9]{1,}_","",ytoP_tree$tip.label)
#substr(ytoP_tree$tip.label, (gregexpr("\\_",ytoP_tree$tip.label)[[1]])[-1], (gregexpr("\\_",ytoP_tree$tip.label)[[1]])[-1])<-"."
ytoP_tree$tip.label<-gsub("\\_[1]{1}$","\\.1",ytoP_tree$tip.label)
ytoP_tree$tip.label<-gsub("\\_[2]{1}$","\\.2",ytoP_tree$tip.label)

character<-setNames(as.data.frame(ytoP_tree$tip.label),"tip.label")
character<-separate(character, col =tip.label ,sep = "\\_{2}", into = "tip.label")
character$tip.label<-gsub("\\_[1]{1}$","\\.1",character$tip.label)

ytoP_tree$tip.label<-character$tip.label

#delete fasta format
names(`ytoP_filtered`)[1]<-"Accession"
`ytoP_filtered`$Accession<-gsub(">","",`ytoP_filtered`$Accession)
names(`ytoP_filtered`)[1]<-"tip.label"



character<-setNames(as.data.frame(ytoP_tree$tip.label),"tip.label")
character<-merge(character,`ytoP_filtered`,by.x ="tip.label",sort=FALSE,all.x = TRUE)
character$species<-gsub("\\[|\\]","",character$species)



#BacDive
BacDiveDB<-rbind(archaea_BacDiveDB,bacteria_BacDiveDB)

test<-merge(character, BacDiveDB,by="species",sort=FALSE,all.x = TRUE)

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

test$temp_opt<-round(as.numeric(test$temp_opt),0)
test$pH_opt<-round(as.numeric(test$pH_opt),1)
test$temp<-round(as.numeric(test$temp),0)
test$pH<-round(as.numeric(test$pH),1)

test<-row_count(test, pH,pH_opt,temp,temp_opt, count = NA, append = TRUE)
test2 <- test %>% group_by(tip.label) %>% filter(rowcount == min(rowcount))
test2 <- cbind(test2, "accession_duplicate"=duplicated(test2$tip.label))
test2<-subset(test2, test2$accession_duplicate==FALSE) 

reference_sequences<-setNames(cbind(template_fasta,"reference"),c("tip.label","sequence","reference"))
reference_sequences<-reference_sequences[,c(1,3)]
reference_sequences$tip.label<-gsub(">","",reference_sequences$tip.label)
test2<-merge(test2, reference_sequences, by="tip.label", sort=FALSE, all.x = TRUE)
test2<-test2 %>% mutate_if(is.character, as.factor)
#test2<-test2[,c(2,1,3:ncol(test2))]
#taxa assignment
{
  for(i in 1){
    protein_search<- entrez_search(db="protein", term=test2$tip.label[i],retmode = "xml", retmax = 999999, use_history = TRUE)
    taxonomy_data <- entrez_link(dbfrom = "protein", id = protein_search$ids, db ="taxonomy")
    taxonomy_id <- xmlToDataFrame(entrez_fetch(db="taxonomy", id=taxonomy_data$links$protein_taxonomy[1], rettype="xml"), stringsAsFactors = FALSE) %>% select("TaxId")
    raw_recs <- taxize::classification(taxonomy_id[[1]], db = 'ncbi')
    #raw_recs <- as.data.frame.list(raw_recs)
    raw_recs <- do.call(rbind.data.frame, raw_recs)
    rownames(raw_recs)<-NULL
    raw_recs<-melt(raw_recs, id=c("rank","id"),variable.name="name")
    raw_recs<-dcast(raw_recs, formula = name ~ rank, fun.aggregate = function(x){paste(x,collapse = "_")})
    raw_recs<-cbind(tip.label=test2$tip.label[i], protein_id=protein_search$ids[1], tax_id=taxonomy_data$links$protein_taxonomy[1], raw_recs)
    taxa_assignment<-raw_recs
  }
  
  for(i in (1:nrow(name))[-1]){
    protein_search<- entrez_search(db="protein", term=test2$tip.label[i], retmode = "xml", retmax = 999999, use_history = TRUE)
    taxonomy_data <- entrez_link(dbfrom = "protein", id = protein_search$ids, db ="taxonomy")
    if(!is.null(taxonomy_data$links$protein_taxonomy[1])){
      taxonomy_id <- xmlToDataFrame(entrez_fetch(db="taxonomy", id=taxonomy_data$links$protein_taxonomy[1], rettype="xml"), stringsAsFactors = FALSE) %>% select("TaxId")
      raw_recs <- taxize::classification(taxonomy_id[[1]], db = 'ncbi')
      #raw_recs <- as.data.frame.list(raw_recs)
      raw_recs <- do.call(rbind.data.frame, raw_recs)
      rownames(raw_recs)<-NULL
      raw_recs<-melt(raw_recs, id=c("rank","id"),variable.name="name")
      raw_recs<-dcast(raw_recs, formula = name ~ rank, fun.aggregate = function(x){paste(x,collapse = "_")})
      raw_recs<-cbind(tip.label=test2$tip.label[i], protein_id=protein_search$ids[1], tax_id=taxonomy_data$links$protein_taxonomy[1], raw_recs)
      #raw_recs <- data.frame(matrix(unlist(raw_recs), nrow=length(raw_recs), byrow=T))
      #filter(rank=c("phylum","class","order","family","genus","species","strain"))
    }
    else{raw_recs<-NULL}
    taxa_assignment<-rbind.fill(taxa_assignment,raw_recs)} 
}

for(i in (i+1:nrow(name))){
  protein_search<- entrez_search(db="protein", term=test2$tip.label[i], retmode = "xml", retmax = 999999, use_history = TRUE)
  taxonomy_data <- entrez_link(dbfrom = "protein", id = protein_search$ids, db ="taxonomy")
  if(!is.null(taxonomy_data$links$protein_taxonomy[1])){
    taxonomy_id <- xmlToDataFrame(entrez_fetch(db="taxonomy", id=taxonomy_data$links$protein_taxonomy[1], rettype="xml"), stringsAsFactors = FALSE) %>% select("TaxId")
    raw_recs <- taxize::classification(taxonomy_id[[1]], db = 'ncbi')
    #raw_recs <- as.data.frame.list(raw_recs)
    raw_recs <- do.call(rbind.data.frame, raw_recs)
    rownames(raw_recs)<-NULL
    raw_recs<-melt(raw_recs, id=c("rank","id"),variable.name="name")
    raw_recs<-dcast(raw_recs, formula = name ~ rank, fun.aggregate = function(x){paste(x,collapse = "_")})
    raw_recs<-cbind(tip.label=test2$tip.label[i], protein_id=protein_search$ids[1], tax_id=taxonomy_data$links$protein_taxonomy[1], raw_recs)
    #raw_recs <- data.frame(matrix(unlist(raw_recs), nrow=length(raw_recs), byrow=T))
    #filter(rank=c("phylum","class","order","family","genus","species","strain"))
  }
  else{raw_recs<-NULL}
  taxa_assignment<-rbind.fill(taxa_assignment,raw_recs)} 


taxa_assignment<-taxa_assignment[c("tip.label","protein_id","tax_id","superkingdom","phylum","class","order","family","genus","species","strain")]
names(taxa_assignment)<-paste("ncbi",c("protein_accession","protein_id","tax_id","superkingdom","phylum","class","order","family","genus","species","strain"),sep="_")
taxa_assignment$ncbi_protein_accession<-gsub(">","",taxa_assignment$ncbi_protein_accession)
names(taxa_assignment)[1]<-"tip.label"

########################################################################################################################################################################
################################################################   phylogeny  ##########################################################################################
########################################################################################################################################################################
#test<-test %>% mutate_if(is.character, as.factor)
#dd<-rbind(primordial,modern,uncertain)
#colnames(dd)<-c("taxa","age")
#test<-test %>% mutate_if(is.factor, as.character)
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

background_color <- unikn::newpal(col = as.character(background_color_list$background_color),
                                  names = as.character(background_color_list$phylum))
unikn::seecol(background_color)
ggsave("phylum_color.pdf",width=10,height=5,units="cm",dpi=600,limitsize=FALSE)

background_color_list<-background_color_list %>% mutate_if(is.factor, as.character)
background_color_list<-background_color_list %>% mutate_if(is.character, as.factor)


#font
#loadfonts()
fonts()
library(extrafont)
loadfonts(device = "win")
windowsFonts("MyriadPro.ttf")

font_import(pattern = "Myriad.*", prompt=FALSE)

ttf_import()


fonttable()



#show_text package save typo as vector
font_add('Myriad Pro', file.choose())

font_add("Myriad Pro",
         regular = "MyriadPro-Regular.otf",
         bold = "MyriadPro-Bold.otf",
         italic = "MyriadPro-It.otf")


#showtext_auto()

#showtext_begin()
#phlyogeny
OO<-ggtree(ytoP_tree, layout="circular", size=1.0, linetype=1)+ xlim(-0.5,NA)+
  geom_tiplab2(aes(label=gsub("\\[|\\]","",paste(p$data$species))), color="black", fontface="italic", align = TRUE, linesize=0, linetype = 0, offset = 0.35,  show.legend = FALSE)+
  geom_tiplab2(aes(color=p3$data$phylum, label=""), fontface="italic", align = TRUE, linesize=0.2, linetype = 3, offset = 0.35,  show.legend = FALSE)+
  scale_color_manual(name="phylum", values=c(background_color), na.value = "grey90")+
  scale_fill_manual(name="phylum", values=c(background_color), na.value = "grey90")
#showtext_end()  
#geom_text(size = 2,align =TRUE, aes(family = "Myriad Pro", angle=angle, color=p3$data$phylum, label=gsub("\\[|\\]","",paste(p$data$species))))+

#as gradient
p<-OO %<+% test2
p2<-p %<+% taxa_assignment
p3<-p2 + aes(color=ncbi_phylum)+
         aes(fill=ncbi_phylum)+
         new_scale("fill")+
         geom_fruit(geom=geom_tile, show.legend = FALSE,
                    mapping=aes(y = ncbi_phylum, fill=ncbi_phylum),linesize=0, offset = 1, alpha=0.1, #if need stacked boxes change the x and y axis
                    #width=5.0, 
                    #pwidth=3.0, 
                    #offset=0.0,
                    stat = "identity",
                    position = "auto",
                    axis.params=list(#axis="x", # add axis text of the layer.
                    #text.angle=-45, # the text size of axis.
                    hjust=0,
                    vjust=0# adjust the horizontal position of text of axis.
             ))+

  #scale_color_discrete("phylum")+
  #guides(color = guide_legend(override.aes = list(size = 6),shape = 21))+
  geom_tiplab(aes(label=p$data$pH_opt), color="black", align = TRUE, hjust=-2.5, linesize = 0, linetype = 0, size = 4, offset = 0) +
  geom_tiplab(aes(label=p$data$temp_opt), color="black", align = TRUE, hjust=-3.75, linesize = 0, linetype = 0, size = 4, offset = 0) +
  
  new_scale("color") +
  geom_tippoint(size=4, shape=16, aes(color=temp_opt, x=1.15)) +
  scale_color_distiller(name="Temperature", palette ="YlOrRd", direction = 1, na.value = NA, limits=c(10,70),breaks=c(10,30,50,70), oob=squish)+theme(legend.position="right")+
  
  new_scale("color") +
  geom_tippoint(size=4, shape=16, aes(color=pH_opt, x=1.00))+
  scale_color_distiller(name="pH",palette ="YlGnBu", direction = 1, na.value = NA, limits=c(3,11),breaks=c(3,5,7,9,11), oob=squish)+theme(legend.position="right")

#new_scale("color") +
#geom_tippoint(size=3, shape=21, aes(fill=reference, x=2.5, na.value=NA), na.rm = TRUE)+
#scale_fill_manual(values = "yellow", na.value = NA)
p3

p4<-p3+ new_scale("fill")+
        geom_star(mapping=aes(fill=reference, starshape=reference), size=3, starstroke=0.0)+
        scale_fill_manual(name="reference", values=c("red"), na.value = NA) + theme(legend.position="right")+
        
  #geom_fruit(reference,
  #           geom=geom_star,
  #           mapping=aes(fill=reference, size=5, starshape=26),
  #           position="identity",
  #           starstroke=0.2)
  
        theme(#legend.position=c(0.96, 0.5), # the position of legend.
              legend.background=element_rect(fill=NA), # the background of legend.
              #legend.title=element_text(size=7), # the title size of legend.
              #legend.text=element_text(size=6), # the text size of legend.
              #legend.spacing.y = unit(0.02, "cm")  # the distance of legends (y orientation).
              ) 

p4
#geom_fruit(geom=geom_tile, show.legend = FALSE,
#           mapping=aes(ncbi_phylum, y = 1, fill=ncbi_phylum), alpha=0.3,
#           #width=5.0, 
#           pwidth=1.0, 
#           #offset=0.0,
#           #stat = "identity",
#           axis.params=list(#axis="x", # add axis text of the layer.
#             #text.angle=-45, # the text size of axis.
#             hjust=0,
#             vjust=0# adust the horizontal position of text of axis.
#           )


ggsave("ytoP_tree.pdf",width=60,height=60,units="cm",dpi=600,limitsize=FALSE)
#showtext_auto(FALSE)


########################################################################################################################################################################


additional_template_fasta<-rbind(template_fasta, ytoP_clustered_50_376[grep("KQB87418.1|KGJ95074.1|KXG78355.1" ,ytoP_clustered_50_376$accession),][,c(1,4)])
additional_template_fasta<-merge(ytoP_clustered_50_376, additional_template_fasta, by="accession", sort=FALSE, all.x = FALSE, no.dups = TRUE, how='outer', suffixes=c('', ''))
additional_template_fasta<-additional_template_fasta[,-ncol(additional_template_fasta)]
additional_template_fasta<-cbind(paste(additional_template_fasta$accession, additional_template_fasta$annotation, additional_template_fasta$species), additional_template_fasta$sequence)
write.table(additional_template_fasta,"ytoP_10template.txt", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)

write.table(template_fasta,"ytoP_10template.txt", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)



`ytoP_aligned`<-read.table("ytoP_clustered_fasta_70_121_aligned", header=FALSE, sep= "\t", quote="")
name <- grep("^>", `ytoP_aligned`$V1, value = TRUE,)
name <- as.data.frame(name)
#species <- separate(name, col=name, into=c("protein number","species"), sep=" \\[", fill = "right")
#accession <- separate(name, col=name, into=c("protein number","species"), sep=" ", fill = "right")
#accession<-as.data.frame(accession$`protein number`)
#colnames(accession)<-"accession"
#species$annotation<-gsub(">.{1,15} ","",species$`protein number`)
seq_full<-grep("^>",`ytoP_aligned`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,nrow(`ytoP_aligned`))
sink("ytoP_aligned_seq.csv");for(i in 1:nrow(`name`)){
  w=paste(`ytoP_aligned`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(w,"\n")
}
ytoP_aligned<-read.csv(file = "ytoP_aligned_seq.csv", header=FALSE, sep="\n")
#ytoP_aligned<-ytoP_aligned[-nrow(ytoP_aligned),]
ytoP_aligned<-cbind(name, ytoP_aligned)
################################################################################################
#sequence<-seqinr::read.alignment(file = "ytoP_clustered_fasta_70_562_aligned", format = "fasta", forceToLower = TRUE)
sequence=na.omit(ytoP_aligned$V1)

sequence<-as.character(sequence)

{
  p1 <- ggseqlogo(sequence, seq_type='aa', method = 'bits', facet = "wrap", nrow =5)+theme(legend.position="none")+theme_logo(base_size = 15)+theme(axis.text.x = element_text(angle = 90, hjust = 0))
  p1$coordinates$limits$x<-c(151,250)
  p1$scales$scales[[1]]<- scale_x_continuous(breaks=c(seq(151,250, by=10)),expand = c(0.005,0.005))
  #10000 1260
  p2 <- ggseqlogo(sequence, seq_type='aa', method = 'bits', facet = "wrap", nrow =5)+theme_logo(base_size = 15)+theme(axis.text.x = element_text(angle = 90, hjust = 0))
  p2$coordinates$limits$x<-c(251,350)
  p2$scales$scales[[1]]<- scale_x_continuous(breaks=c(seq(251,350, by=10)),expand = c(0.005,0.005)) 
  #2110 2160
}


gridExtra::grid.arrange(p1, p2, nrow=2, ncol=1)
#ggsave("ytoP_seqlogo.pdf",units="cm",dpi=300,limitsize=FALSE)

ggseqlogo <- arrangeGrob(p1, p2, nrow=2, ncol=1)

#save
ggsave("ytoP_seqlogo.pdf",ggseqlogo,width=50,height=30,units="cm",dpi=300,limitsize=FALSE)







