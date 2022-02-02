pacman::p_load(ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr,ggtree, magrittr,sjmisc, ggseqlogo, gridExtra, ggtree, ggtreeExtra, ggstar, ggnewscale, 
               xlsx, openxlsx, seqinr, ape, rentrez, reutils, colorspace, scales, taxonomizr, phyloR, XML, xml2, reshape2, gtools, data.table, unikn, nord, purrr, inlmisc,
               tidyverse, showtext, extrafont, myriad, Rttf2pt1)

##synonym
#AcpII, Ak.1 protease, Alcalase, Alcalase 0.6L, Alcalase 2.5L, ALE1 subtilase, ALK-enzyme, 
#Alkaline mesentericopeptidase, Alkaline protease, alkaline serine protease, ALP1, Alzwiprase, 
#aprE, AprE51, aqualysin, Arp, AsES, Asp v 13, AtSBT1.9, Bacillopeptidase A, Bacillopeptidase B, Bacillus gibsonii alkaline protease, Bacillus subtilis alkaline proteinase Bioprase, 
#BgAP, Bioprase AL 15, Bioprase APL 30, BLS, BPN', BprB, BprV, C1 subtilase, cold active subtilisin-like serine proteinase, Colistinase, EC 3.4.21.14, EC 3.4.4.16, Esperase, 
#Fe protease, Fe prtS8A, Genenase I, intracellular subtilisin protease, ISP, IvaP, Kazusase, Maxatase, mesenteroicopeptidase, More, Nagarse, Opticlean, ORF2, Orientase 10B, P69 subtilase, 
#PBANKA_1106900, PbSOPT, Peptidase, subtilo-, A, PF3D7_0507300, phytophase, PIMMS2, Protease S, Protease VIII, Protease XXVII, Proteinase K, Proteinase, Bacillus subtilis alkaline, Protin A 3L, 
#PSP-3, psychrophilic subtilisin-like protease, S1P subtilase, SAP, SAS-1, SASP, saspase, Savinase, Savinase 16.0L, Savinase 32.0 L EX, Savinase 4.0T, Savinase 8.0L, savinaseTM, 
#SBc, SBL, SDD1 subtilase, senescence-associated subtilisin protease, SES7, SISBT3 subtilase, SOPT, SP 266, Sspa, SSU0757, SUB1, SUB2, subC, subtilase, subtilase subfamily 1 member 9, 
#subtilase-like protease, subtilisin, subtilisin 72, subtilisin A, Subtilisin amylosacchariticus, Subtilisin BL, subtilisin BPN??, subtilisin C., subtilisin Carlsberg, subtilisin DJ-4, 
#Subtilisin DY, Subtilisin E, subtilisin E-S7, Subtilisin GX, subtilisin JB1, subtilisin Karlsberg, Subtilisin Novo, subtilisin Pr1-like protease, subtilisin protease, subtilisin QK, 
#Subtilisin S41, subtilisin S4I, subtilisin S88, Subtilisin Sendai, subtilisin Sph, subtilisin-like ookinete protein, subtilisin-like ookinete protein important for transmission, 
#subtilisin-like protease, subtilisin-like protease AprV2, subtilisin-like serine protease, Subtilisn J, Subtilopeptidase, Superase, thermitase, thermo-active subtilisin-like serine protease, 
#Thermoase, Thermoase PC 10, thermophilic thermitase, thermostable subtilisin, ThSS45, Tk-subtilisin, trans-cinnamoyl-subtilisin, V. cholerae-secreted serine protease, VC_0157, vPR




#aprE / subtilisin /3.4.21.62[EC/RN Number=] 
{
  `aprE`<- entrez_search(db="protein", term="subtilisin",retmode = "xml", retmax = 999999, use_history = TRUE)
  for(seq_start in seq(1,length(`aprE`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`aprE`$web_history, rettype="fasta", retmax=500, retstart=seq_start)
    cat(recs, file="aprE_135148.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}

{
  `aprE`<- entrez_search(db="protein", term="3.4.21.62[EC]",retmode = "xml", retmax = 999999, use_history = TRUE)
  for(seq_start in seq(1,length(`aprE`$ids),500)){
    recs <- entrez_fetch(db="protein", web_history=`aprE`$web_history, rettype="fasta", retmax=500, rets3www3tart=seq_start)
    cat(recs, file="aprE_9999.fasta", append=TRUE)
    cat(seq_start+499, "sequences downloaded\r")
  }}




#data partitioning
`aprE`<-read.table("aprE_9999_EC.fasta", header=FALSE, sep= "\t", quote="")
name <- grep("^>", `aprE`$V1, value = TRUE,)
name <- as.data.frame(name)
seq_full<-grep("^>",`aprE`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,nrow(`aprE`))
sink("aprE_seq.csv");for(i in 1:nrow(`name`)){
  w=paste(`aprE`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(w,"\n")
}
`aprE_sequence`<-read.csv("aprE_seq.csv",header=FALSE,sep="\n")
`aprE`<-cbind(name,`aprE_sequence`)

rm(list=ls(pattern=".{1,}_sequence$")) 
rm(name)

`aprE`<- separate(`aprE`, col=name, into=c("accession","species"), sep=" ", fill = "right", extra = "merge")

##filter Swiss-Prot/UniProtKB
`aprE_filtered`<-separate(`aprE`, col=species, into=c("annotation","species"), sep=" \\[", fill = "right", extra = "merge")
#`aprE_filtered`<-separate(`aprE_filtered`, col=species, into=c("annotation","species"),sep="\\[(?=[^\\[]+$)",fill="left", extra="merge")
`aprE_filtered`$species<-paste0("[",`aprE_filtered`$species)
names(`aprE_filtered`)[4]<-"sequence"
aprE_filtered$sequence<-as.character(aprE_filtered$sequence)
`aprE_filtered`<-mutate(`aprE_filtered`, number =nchar(`aprE_filtered`$sequence))

#template sequence
template_sequence<-filter(`aprE_filtered`,grepl("NP_388911.2|AEH48379.1|EUB99404.1|CUB36377.1|ADV66372.1|SNQ60367.1|KOX89359.1|GBF57798.1|EOZ98320.1|RIH86974.1",accession))

#CAA24990.1: [Bacillus amyloliquefaciens] ->CUB36377.1
#AAF12593.1: [Deinococcus maricopensis DSM 21211]	-> ADV66372.1
#AFV22290.1: [Methanolobus psychrophilus R15] ->SNQ60367.1
#ADW21000.1: [Thermus scotoductus SA-01] ->KOX89359.1
#GBF57798.1 Candidatus Phycosocius bacilliformis
#EOZ98320.1 [Indibacter alkaliphilus LW1]
#RIH86974.1 [Calidithermus terrae]


`aprE_filtered`<-filter(`aprE_filtered`, !grepl("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}",accession))
`aprE_filtered`<-filter(`aprE_filtered`, !grepl("[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]{1,2}",accession))
`aprE_filtered`<-filter(`aprE_filtered`, !grepl("pdb\\|[0-9][A-Z0-9]{3}",accession))
`aprE_filtered`<-filter(`aprE_filtered`, !grepl("patent",species))
##RefSeq accession numbers and molecule types.
#AP_	Protein	Annotated on AC_ alternate assembly
#NP_	Protein	Associated with an NM_ or NC_ accession
#YP_c	Protein	Annotated on genomic molecules without an instantiated transcript record
#XP_c	Protein	Predicted model, associated with an XM_ accession
#WP_	Protein	Non-redundant across multiple strains and species

##filter empty species info; pdb, US patent
`aprE_filtered`<-filter(`aprE_filtered`, !is.na(species))
`aprE_filtered`<-filter(`aprE_filtered`, !grepl("hypothetical|unnamed|partial|predicted|putative|LOW QUALITY PROTEIN|Uncharacterised|Uncharacterized|PREDICTED|UNVERIFIED|fragment|precursor|inhibitor|patent|probable|possible",annotation,ignore.case=TRUE))
`aprE_filtered`<-filter(`aprE_filtered`, !grepl("phosphatase|reductase|metal|elastase|xylem|meiotic|regulatory|convertase|cyanobactin|similarity|intracellular|toxin|TPA|trypsin|interleukin|phage|sensor|=|processing|anacrylamide|response",annotation,ignore.case=TRUE))
`aprE_filtered`<-filter(`aprE_filtered`, grepl("subtilisin|protease|S8|aprE|serine|alk",annotation,ignore.case=TRUE))
summary<-summary(`aprE_filtered`$number)

annotation_list<-as.data.frame(unique(aprE_filtered$annotation))


##sequence length
#`aprE_filtered`<-filter(`aprE_filtered`, number >358, ignore.case=TRUE)
`aprE_filtered`<-filter(`aprE_filtered`, number >270 & number<900 ,ignore.case=TRUE)


##PRATT motif (expasy prosite; PDOC00125: Serine proteases, subtilase family, active sites signatures and domain profile)
#SUBTILASE_ASP, PS00136; Serine proteases, subtilase family, aspartic acid active site  (PATTERN)
#[STAIV]-{ERDL}-[LIVMF]-[LIVM]-D-[DSTA]-G-[LIVMFC]-x(2,3)-[DNH]
#`aprE_filtered`<-filter(`aprE_filtered`, grepl("([STAIV][^ERDL][LIVMF][LIVM]D[DSTA]G[LIVMFC].{2,3}[DNH])",sequence))
#SUBTILASE_HIS, PS00137; Serine proteases, subtilase family, histidine active site  (PATTERN)
#H-G-[STM]-x-[VIC]-[STAGC]-[GS]-x-[LIVMA]-[STAGCLV]-[SAGM]
#`aprE_filtered`<-filter(`aprE_filtered`, grepl("(HG[STM].[VIC][STAGC][GS].[LIVMA][STAGCLV][SAGM])",sequence))
#SUBTILASE_SER, PS00138; Serine proteases, subtilase family, serine active site  (PATTERN)
#G-T-S-x-[SA]-x-P-x-{L}-[STAVC]-[AG]
#`aprE_filtered`<-filter(`aprE_filtered`, grepl("(GTS.[SA].P.[^L][STAVC][AG])",sequence))
`aprE_filtered`<-filter(`aprE_filtered`, grepl("([STAIV][^ERDL][LIVMF][LIVM]D[DSTA]G[LIVMFC].{2,3}[DNH])|(HG[STM].[VIC][STAGC][GS].[LIVMA][STAGCLV][SAGM])|(GTS.[SA].P.[^L][STAVC][AG])",sequence))

##delete duplicated sequence
aprE_filtered<-rbind(template_sequence, aprE_filtered)
`aprE_filtered`<-cbind(`aprE_filtered`, "duplicate"=duplicated(`aprE_filtered`$sequence))
`aprE_filtered`<-subset(`aprE_filtered`, `aprE_filtered`$duplicate==FALSE) 


##template accession
#CAA24990.1=Bacillus amyloliquefaciens
#NP_388911.2=Bacillus subtilis subsp. subtilis str. 168
#AAF12593.1=Deinococcus radiodurans R1
#AEH48379.1=Geobacillus thermoglucosidasius C56-YS93
#AFV22290.1=Methanolobus psychrophilus R15
#EUB99404.1=Rhizobium sp. CF080
#ADW21000.1=Thermus scotoductus SA-01


#re-confirm contaminant sequences
annotation_list<-as.data.frame(unique(aprE_filtered$annotation))

##manually filterate aprE gene (if needed)
#`aprE_filtered`<-filter(`aprE_filtered`, !grepl("autophagic|",annotation, ignore.case = TRUE))
#`aprE_filtered`<-filter(`aprE_filtered`, !grepl("catalase/peroxidase",annotation, ignore.case = TRUE))
#`aprE_filtered`<-filter(`aprE_filtered`, !grepl("SpoS|TPA|similar|spore|UNVERIFIED|probable|domain",annotation, ignore.case = TRUE))
#`aprE_filtered`<-filter(`aprE_filtered`, !grepl("(SWPDN).*(VQMGLIYV).*",sequence))


fasta<-as.data.frame(paste(`aprE_filtered`$accession))
fasta<-setNames(cbind(fasta,`aprE_filtered`$sequence),c("name","sequence"))
template_fasta<-as.data.frame(paste(`template_sequence`$accession))
template_fasta<-setNames(cbind(template_fasta,`template_sequence`$sequence),c("name","sequence"))
fasta<-rbind(template_fasta, fasta)
`fasta`<-cbind(`fasta`, "duplicate"=duplicated(`fasta`$sequence))
`fasta`<-subset(`fasta`, `fasta`$duplicate==FALSE) 
fasta<-fasta[,-3]
write.table(`fasta`,"aprE_EC_831_fasta.txt", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)


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
  assign(paste0("aprE_clustered",src_file_cluster[i]),read.table(src_file[i], header=FALSE, quote=""))
  assign("name",grep("^>", get(paste0("aprE_clustered",src_file_cluster[i]))$V1, value = TRUE))
  assign("name",as.data.frame(name))
  assign("seq_full",grep("^>",get(paste0("aprE_clustered",src_file_cluster[i]))$V1))
  assign("seq_start",seq_full+1)
  assign("seq_finish",c(seq_full[-1]-1,nrow(get(paste0("aprE_clustered",src_file_cluster[i])))))
  
  sink(paste0("aprE_clustered",src_file_cluster[i],"_seq.csv"));for(j in 1:nrow(`name`)){
    assign("w",paste(get(paste0("aprE_clustered",src_file_cluster[i]))$V1[seq_start[j]:seq_finish[j]],collapse = ""))
    cat(w,"\n")
  }
  assign(paste0("aprE_clustered",src_file_cluster[i]), read.csv(paste0("aprE_clustered",src_file_cluster[i],"_seq.csv"),header=FALSE,sep="\n"))
  assign(paste0("aprE_clustered",src_file_cluster[i]), cbind(name,get(paste0("aprE_clustered",src_file_cluster[i]))))
  assign(paste0("aprE_clustered",src_file_cluster[i]), setNames(get(paste0("aprE_clustered",src_file_cluster[i])),c("accession","sequence")))
}


#add reference sequences and delete duplicated reference seq
#template_aligned_sequence<-filter(`fasta`,grepl("AAB47761.1|NP_391784.2|WP_000077872.1|AAC18407.1|WP_014589250.1|CEJ79573.1|TMW79725.1",name))
#add template sequences
template_fasta<-as.data.frame(paste(`template_sequence`$accession))
template_fasta<-setNames(cbind(template_fasta,`template_sequence`$sequence),c("accession","sequence"))

for (i in 1:length(src_file_cluster)){
  assign(paste0("aprE_clustered",src_file_cluster[i]),rbind(template_fasta, get(paste0("aprE_clustered",src_file_cluster[i]))))
  assign(paste0("aprE_clustered",src_file_cluster[i]),cbind(get(paste0("aprE_clustered",src_file_cluster[i])), "duplicate"=duplicated(get(paste0("aprE_clustered",src_file_cluster[i]))$sequence)))
  assign(paste0("aprE_clustered",src_file_cluster[i]),subset(get(paste0("aprE_clustered",src_file_cluster[i])), get(paste0("aprE_clustered",src_file_cluster[i]))$duplicate==FALSE))
  assign(paste0("aprE_clustered",src_file_cluster[i]),get(paste0("aprE_clustered",src_file_cluster[i]))[,-3])
  assign(paste0("aprE_clustered",src_file_cluster[i]),merge(`aprE_filtered`,get(paste0("aprE_clustered",src_file_cluster[i])),by="accession", sort=FALSE,all.x = FALSE, no.dups = TRUE, how='outer', suffixes=c('', '_y')))
  assign(paste0("aprE_clustered",src_file_cluster[i]),get(paste0("aprE_clustered",src_file_cluster[i]))[,-7])
}

for(i in 1:length(src_file_cluster)){
  assign(paste0("aprE_clustered_fasta",src_file_cluster[i]), cbind(name=get(paste0("aprE_clustered",src_file_cluster[i]))$accession,get(paste0("aprE_clustered",src_file_cluster[i]))$sequence))
  write.table(get(paste0("aprE_clustered_fasta",src_file_cluster[i])),paste0("aprE_clustered_fasta",src_file_cluster[i]), row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)
}

########################################################################################################################################################################
################################################################   phylogeny  ##########################################################################################
########################################################################################################################################################################
#crawling newick file
src_dir<-getwd()
tree_file <- list.files(src_dir, pattern = "\\.on\\.nh$") # list

aprE_tree<-read.tree("aprE_EC50_nj.201116191719286.on.nh")
aprE_tree$tip.label<-gsub("^[0-9]{1,}_","",aprE_tree$tip.label)
#substr(aprE_tree$tip.label, (gregexpr("\\_",aprE_tree$tip.label)[[1]])[-1], (gregexpr("\\_",aprE_tree$tip.label)[[1]])[-1])<-"."
aprE_tree$tip.label<-gsub("\\_[1]{1}$","\\.1",aprE_tree$tip.label)
aprE_tree$tip.label<-gsub("\\_[2]{1}$","\\.2",aprE_tree$tip.label)

character<-setNames(as.data.frame(aprE_tree$tip.label),"tip.label")
character<-separate(character, col =tip.label ,sep = "\\_{2}", into = "tip.label")
character$tip.label<-gsub("\\_[1]{1}$","\\.1",character$tip.label)

aprE_tree$tip.label<-character$tip.label

#delete fasta format
names(`aprE_filtered`)[1]<-"Accession"
`aprE_filtered`$Accession<-gsub(">","",`aprE_filtered`$Accession)
names(`aprE_filtered`)[1]<-"tip.label"

character<-setNames(as.data.frame(aprE_tree$tip.label),"tip.label")
character<-merge(character,`aprE_filtered`,by.x ="tip.label",sort=FALSE,all.x = TRUE)
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
  
  for(i in (1:nrow(test2))[-1]){
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

#manual color
rm(background_color, background_color_list)
length(unique(taxa_assignment$ncbi_phylum))
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

#phlyogeny
OO<-ggtree(aprE_tree, layout="circular", size=0.5, linetype=1)+ xlim(-0.5,NA)+
  geom_tiplab2(aes(label=gsub("\\[|\\]","",paste(p$data$species))), color="black", fontface="italic", align = TRUE, size=1.8, linesize=0, linetype = 0, offset = 0.7,  show.legend = FALSE)+
  geom_tiplab2(aes(color=p3$data$phylum, label=""), fontface="italic", align = TRUE, size=1.8,linesize=0.2, linetype = 3, offset = 0.7,  show.legend = FALSE)+
  scale_color_manual(name="phylum", values=c(background_color), na.value = "grey90")+
  scale_fill_manual(name="phylum", values=c(background_color), na.value = "grey90")

#as gradient
p<-OO %<+% test2
p2<-p %<+% taxa_assignment
p3<-p2 + aes(color=ncbi_phylum)+
         aes(fill=ncbi_phylum)+
         geom_fruit(geom=geom_tile, show.legend = FALSE,
             mapping=aes(y = ncbi_phylum, fill=ncbi_phylum),linesize=0, offset = 1.0, alpha=0.1, #if need stacked boxes change the x and y axis
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
  new_scale("fill") +
  #scale_color_discrete("phylum")+
  #guides(color = guide_legend(override.aes = list(size = 6),shape = 21))+
  geom_tiplab(aes(label=p$data$pH_opt), color="black", align = TRUE, hjust=-9.0, linesize = 0, linetype = 0, size = 2, offset = 0) +
  geom_tiplab(aes(label=p$data$temp_opt), color="black", align = TRUE, hjust=-12.0, linesize = 0, linetype = 0, size = 2, offset = 0) +
  
  new_scale("color") +
  geom_tippoint(size=2.0, shape=16, aes(color=pH_opt, x=1.55))+
  scale_color_distiller(name="pH",palette ="YlGnBu", direction = 1, na.value = NA, limits=c(3,11),breaks=c(3,5,7,9,11), oob=squish)+theme(legend.position="right")+
  
  new_scale("color") +
  geom_tippoint(size=2.0, shape=16, aes(color=temp_opt, x=1.85)) +
  scale_color_distiller(name="Temperature", palette ="YlOrRd", direction = 1, na.value = NA, limits=c(10,70),breaks=c(10,30,50,70), oob=squish)+theme(legend.position="right")+
  new_scale("fill")

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

ggsave("aprE_tree.pdf",width=60,height=60,units="cm",dpi=600,limitsize=FALSE)

#Addition if needed
#KQB87418.1 [Corynebacterium lowii]	
#KGJ95074.1 [Colwellia psychrerythraea]	
#KXG78355.1	[Fervidicola ferrireducens]	
write.table(template_fasta,"aprE_10template.txt", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)

additional_template_fasta<-rbind(template_fasta, aprE_clustered_50_376[grep("KQB87418.1|KGJ95074.1|KXG78355.1" ,aprE_clustered_50_376$accession),][,c(1,4)])
additional_template_fasta<-merge(aprE_clustered_50_376, additional_template_fasta, by="accession", sort=FALSE, all.x = FALSE, no.dups = TRUE, how='outer', suffixes=c('', ''))
additional_template_fasta<-additional_template_fasta[,-ncol(additional_template_fasta)]
additional_template_fasta<-cbind(paste(additional_template_fasta$accession, additional_template_fasta$annotation, additional_template_fasta$species), additional_template_fasta$sequence)
write.table(additional_template_fasta,"aprE_10template.txt", row.names = FALSE, col.names = FALSE, sep = "\n", quote = FALSE)



devtools::install_github("ropenscilabs/ochRe")
library(ochRe)
pal_names <- names(ochre_palettes)
par(mfrow=c(length(ochre_palettes)/2, 2), lheight = 2, mar=rep(1, 4), adj = 0)
for (i in 1:length(ochre_palettes)){
  viz_palette(ochre_palettes[[i]], pal_names[i])
}



#trimmed gapped sequence
msa <- readFasta("aprE_clustered_fasta_70_491_aligned")
msa<-mutate(msa,"aligned_length"=str_length(msa$Sequence))
msa.trimmed <- msaTrim(msa)
msa.trimmed<-mutate(msa.trimmed,"aligned_length"=str_length(msa.trimmed$Sequence))
#gap.end = 0.5, gap.mid = 0.9
#msa.mat <- msa2mat(msa)  # for use with ape::as.DNAbin(msa.mat)
sequence=na.omit(msa.trimmed$Sequence)
sequence<-as.character(sequence)
x<-as.data.frame(str_locate(sequence,"GTSMA"))



`aprE_aligned`<-read.table("aprE_clustered_fasta_70_491_aligned", header=FALSE, sep= "\t", quote="")
name <- grep("^>", `aprE_aligned`$V1, value = TRUE,)
name <- as.data.frame(name)
#species <- separate(name, col=name, into=c("protein number","species"), sep=" \\[", fill = "right")
#accession <- separate(name, col=name, into=c("protein number","species"), sep=" ", fill = "right")
#accession<-as.data.frame(accession$`protein number`)
#colnames(accession)<-"accession"
#species$annotation<-gsub(">.{1,15} ","",species$`protein number`)
seq_full<-grep("^>",`aprE_aligned`$V1) 
seq_start<-seq_full+1
seq_finish<-c(seq_full[-1]-1,nrow(`aprE_aligned`))
sink("aprE_aligned_seq.csv");for(i in 1:nrow(`name`)){
  w=paste(`aprE_aligned`$V1[seq_start[i]:seq_finish[i]],collapse = "")
  cat(w,"\n")
}
aprE_aligned<-read.csv(file = "aprE_aligned_seq.csv", header=FALSE, sep="\n")
#aprE_aligned<-aprE_aligned[-nrow(aprE_aligned),]
aprE_aligned<-cbind(name, aprE_aligned)
################################################################################################
#sequence<-seqinr::read.alignment(file = "aprE_clustered_fasta_70_562_aligned", format = "fasta", forceToLower = TRUE)
sequence=na.omit(aprE_aligned$V1)

sequence<-as.character(sequence)

{
  p1 <- ggseqlogo(sequence, seq_type='aa', method = 'bits', facet = "wrap", nrow =5)+theme(legend.position="none")+theme_logo(base_size = 15)+theme(axis.text.x = element_text(angle = 90, hjust = 0))
  p1$coordinates$limits$x<-c(381,623)
  p1$scales$scales[[1]]<- scale_x_continuous(breaks=c(seq(381,623, by=20)),expand = c(0.005,0.005))
  
  p2 <- ggseqlogo(sequence, seq_type='aa', method = 'bits', facet = "wrap", nrow =5)+theme_logo(base_size = 15)+theme(axis.text.x = element_text(angle = 90, hjust = 0))
  p2$coordinates$limits$x<-c(3901,4900)
  p2$scales$scales[[1]]<- scale_x_continuous(breaks=c(seq(3901,4900, by=100)),expand = c(0.005,0.005)) 
}
gridExtra::grid.arrange(p1, p2, nrow=2, ncol=1)
#ggsave("aprE_seqlogo.pdf",units="cm",dpi=300,limitsize=FALSE)

ggseqlogo <- arrangeGrob(p1, p2, nrow=2, ncol=1)

#save
ggsave("aprE_seqlogo.pdf",ggseqlogo,width=50,height=30,units="cm",dpi=300,limitsize=FALSE)






















