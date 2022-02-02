#for phylum coloroing
pacman::p_load(ggplot2, ggnewscale, GGally, ggtree, tidyr, plyr, dplyr, stringr,ggtree, magrittr, sjmisc, ggseqlogo, gridExtra, ggtree, ggtreeExtra, ggstar, seqinr, ape, rentrez, reutils, colorspace, scales, taxonomizr, phyloR, XML, xml2, reshape2, gtools, BBmisc, magrittr, sjmisc, unikn)

#katE
rm(list=setdiff(ls(), c("katE_filtered",
                        "catA_filtered",
                        "aprE_filtered",
                        "araA_filtered",
                        "eglS_filtered",
                        "ytoP_filtered",
                        "gabT_filtered",
                        "wecB_filtered",
                        "tpiA_filtered",
                        "pyk_filtered")))

#all detected phylum
phylum_filtered<-rbind.fill(`katE_filtered`,
                            `catA_filtered`,
                            `aprE_filtered`,
                            `araA_filtered`,
                            `eglS_filtered`,
                            `ytoP_filtered`,
                            `gabT_filtered`,
                            `wecB_filtered`,
                            `tpiA_filtered`,
                            `pyk_filtered`)


#color assignment
phylum_list<-as.data.frame(unique(phylum_filtered$phylum))
colnames(phylum_list)<-"phylum"
phylum_list$phylum<-as.character(phylum_list$phylum)
phylum_list<-setNames(as.data.frame(phylum_list[order(phylum_list$phylum),]),"phylum")
hcl.pals()
hcl<-c(hcl.colors(n=nrow(phylum_list)-1, palette =c("spectral")),"#d3d3d3")
#brew<-brewer.pal(n=nrow(background_color_list), name="Spectral")
background_color_list<-cbind(phylum_list, hcl)

colnames(background_color_list)<-c("phylum","background_color")
background_color_list<-as.data.frame(background_color_list)
background_color_list$phylum[nrow(phylum_list)]<-"NA"


background_color <- newpal(col = as.character(background_color_list$background_color),
                           names = as.character(background_color_list$phylum))

seecol(background_color)
ggsave("phylum_color.pdf",width=10,height=5,units="cm",dpi=600,limitsize=FALSE)

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

