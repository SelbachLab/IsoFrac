library(seqinr)
library(data.table)
if(!exists("fa3")){
  fa <- read.fasta("/Users/henno/Documents/Skripte/R_Scripts/DataAnalysis/20210129_AK_ProteinIsoforms/20250821_Rebuttal/Fasta_Type/HS_UP000005640_gff_compiled_hz20220603.fa",seqtype = "AA",as.string = T)
  fa2 <- read.fasta("/Users/henno/Documents/Skripte/R_Scripts/DataAnalysis/20210129_AK_ProteinIsoforms/20250821_Rebuttal/Fasta_Type/uniprot-proteome-UP000005640-Human.fasta",seqtype = "AA",as.string = T)
  
  fa2[grep("Q8TC07",names(fa2))]
  fa3 <- append(fa,fa2)
  fa3 <- fa2
  # dt <- rbindlist(lapply(strsplit(names(fa3),"#"),function(x){data.table(ID=x[1],Type=x[2])}) )
  # table(dt$Type)
  G <- names(fa3)
  Attributes_G <- sapply(fa3,function(x){attributes(x)$Annot})
  Attributes_G_gene <- gsub(".*.GN=|.*.#.*.","",Attributes_G)
  Attributes_G_gene <- gsub(" .*.","",Attributes_G_gene)

  
  
  
  
  Gsplit <- strsplit(G,"|",fixed=T)
  Gdt <- as.data.table(Gsplit)
  Gdt <- t(Gdt)
  Gdt <- data.table(Gdt)
  
  Gdt[,AS:=any(grepl("-[1-9]",.BY$V2)),V2]
  Gdt$Sequence <- unlist(fa3)
  table(Gdt$AS)
  Gdt$Gene <- Attributes_G_gene
  
  Gdt[,Clean_Accession:=gsub("_.*.|-.*.","",V2),]
  Gdt[,Gene_2:=paste(sort(unique(Gene[Gene!=""])),collapse=";"),Clean_Accession]
  table(Gdt$Gene_2=="")
  table(is.na(Gdt$Gene_2))
  
  Gdt$Gene_2[Gdt$Gene_2==""]# <- Gdt$Clean_Accession[Gdt$Gene_2==""]
  Gdt[,CountPerGene:=length(unique(MSequence)),.(Gene_2)]
  table(Gdt$CountPerGene)/length(Gdt$CountPerGene)*100
  

  
  Gdt[,c("sptr","Start","Stop"):={
    hm <- unlist(strsplit(.BY$V1,"_"))
    if(length(hm)==3){
      as.list(hm)
    }else{
      list(as.double(NA),as.double(NA),as.double(NA))
    }
    
  },V1]


  
  temp <- GdtBig[grep("Q8TC07",V2)][grep("MAAAGVVSGKIIYEQEGVYIHSSCGKTNDQDGLISGILRVLEKDAEVIVDWRPLDDALDSSSILYARKDSSSVVEWTQAPKERGHRGSEHLNSYEAEWDMVNTVSFKRKPHTNGDAPSHRNGKSKWSFLFSLTDLKSIKQNKEGMGWSYLVFCLKDDVVLPALHFHQGDSKLLIESLEKYVVLCESPQDKRTLLVNCQNKSLSQSFENLLDEPAYGLIQAGLLDRRKLLWAIHHWKKIKKDPYTATMIGFSKVTNYIFDSLRGSDPSTHQRPPSEMADFLSDAIPGLKINQQEEPGFEVITRIDLGERPVVQRREPVSLEEWTKNIDSEGRILNVDNMKQMIFRGGLSHALRKQAWKFLLGYFPWDSTKEERTQLQKQKTDEYFRMKLQWKSISQEQEKRNSRLRDYRSLIEKDVNRTDRTNKFYEGQDNPGLILLHDILMTYCMYDFDLGYVQGMSDLLSPLLYVMENEVDAFWCFASYMDQMHQNFEEQMQGMKTQLIQLSTLLRLLDSGFCSYLESQDSGYLYFCFRWLLIRFKREFSFLDILRLWEVMWTELPCTNFHLLLCCAILESEKQQIMEKHYGFNEILKHINELSMKIDVEDILCKAEAISLQMVKCKELPQAVCEILGLQGSEVTTPDSDVGEDENVVMTPCPTSAFQSNALPTLSASGARNDSPTQIPVSSDVCRLTPA",Sequence)]
  temp2 <- Gdt[grep("Q8TC07",V2)]
  
  GdtBig <- Gdt
  
  GdtBig[,MSequence:=gsub("^M","",Sequence),Sequence]
  

  Gdt <- GdtBig[,{
    temp <- .SD
    temp$ProteinCollect <- paste(sort(unique(temp$V2)),collapse=";")
    
    
    
    gr <- .BY
    if(dim(temp)[1]>1){
      if(any(temp$UniProtType=="sp")){
        temp <- temp[UniProtType=="sp"]
        if(dim(temp)[1]>1&any(V1=="sp")){
          temp <- temp[V1=="sp"]
        }
      }else{
        temp <- temp
        temp <- temp[CountPerGene==max(CountPerGene)]
        # print("Unknown Exception")
        # stop()
      }
    }
    
    if(dim(temp)[1]>1){
      # print("To Big Table")
      # temp1 <<- temp
      temp <- temp[V2==Clean_Accession]
      # stop()
    }
    temp
  },MSequence]
  
  
  Gdt[,FirstAA:={
    tempS <- Sequence#[SeqN==MaxSeqN]
    tempS <- tempS[nchar(tempS)==max(nchar(tempS),na.rm=T)]
    unlist(strsplit(tempS,""))[1]
  },Gene_2]

#  ggplot(Gdt,aes(MaxSeqN/as.numeric(Stop),as.numeric(Start)))+geom_point()
  Gdt[,SeqN:=(nchar(Sequence)),Gene_2]
  Gdt[,MaxSeqN:=max(SeqN,na.rm=T),Gene_2]
  hist(log10(Gdt$MaxSeqN/Gdt$SeqN))
  table(Gdt$SeqN)
 
}

# gff <- fread('/Users/henno/Documents/Skripte/R_Scripts/DataAnalysis/20210129_AK_ProteinIsoforms/fa/uniprot-filtered-proteome%3AUP000005640+AND+organism%3A%22Homo+sapiens+%28Human%--.gff.txt',sep="\t",fill=T,skip=2)
# gff[,V1=="Chain",V1]
# table(gff$V3=="Chain")

# barplot(table(gff$V3),mpg=2,las=2)
Gdt$Type <- "UNDEFINED"
Gdt[,Type:=gsub(".*.#","",.BY$V3),V3]
Gdt$Type[grep("_",Gdt$Type)] <- ""
Gdt$Type[Gdt$Type==""&Gdt$AS==1] <- "AS"
Gdt[,UniProtType:=gsub("_.*.","",.BY$V1),V1]
table(Gdt$Type)
# table(Gdt[,any(Type_2=="Chain"),Gene_2]$V1)
Gdt$CleanCategory <-"HMPF"
ggplot(Gdt,aes(CleanCategory,log10(MaxSeqN/SeqN)))+geom_boxplot()+theme(axis.text.x = element_text(angle=45))

## Function
library(data.table)
library(ggplot2)
library(ggrepel)

Gdt <- Gdt[!is.na(V1)]
gff[grep("P24298",V2)]

Gdt[,Type:={
  temp <- .SD
  
  if(dim(.SD)[1]>1|all(temp$Type=="")){
    if(any(temp$Type==""&temp$UniProtType=="tr"&temp$CountPerGene>1)){
      temp$Type[temp$Type==""&temp$UniProtType=="tr"] <- "Unassigned Variant"
    }
    if(any(temp$Type==""&temp$UniProtType=="sp"&temp$Type!="No Variant")){
      temp$Type[temp$Type==""&temp$UniProtType=="sp"& temp$Type!="No Variant"] <- "Canonical"
    }
  }
  temp$Type
},Gene_2]

temp <- Gdt[grep("P24298",V2)]
temp$Type
table(temp$Type)

# stop()
# Gdt[Start==2&Stop!=MaxSeqN&Type=="Chain"]$Type <- "Methionine"
# Gdt[Type=="Methionine"]
# Gdt[Start==2&Stop!=MaxSeqN&Type_2=="Chain"]$Type_2 <- "Methionine"

ggplot(Gdt,aes(CleanCategory,log10(MaxSeqN/SeqN)))+geom_boxplot()+theme(axis.text.x = element_text(angle=45))


Gdt[CountPerGene==1]$Type <- "No Variant"
Gdt[,Type_2:={
  if(any(.BY$UniProtType=="sp")){
    Type[Type=="Canonical"] <- "AS"
  }
  Type
  
},.(Gene_2,UniProtType)]
table(Gdt$Type_2)


ggplot(Gdt,aes(CleanCategory,log10(MaxSeqN/SeqN)))+geom_boxplot()+theme(axis.text.x = element_text(angle=45))

Gdt[,CleanCategory:={
  other <- c("Peptide","Transit","Propeptide","Natural")
  cate <- .BY$Type
  if(any(cate==other)){
    cate <- "Other"
  }
  if(any(cate=="Alternative")){
    cate <- "Alternative Translation"
  }
  cate
},Type_2]
table(Gdt$CleanCategory)
ggplot(Gdt,aes(CleanCategory,log10(MaxSeqN/SeqN)))+geom_boxplot()+theme(axis.text.x = element_text(angle=45))

### Plot
DT <- Gdt[,{
  ta <- table(CleanCategory)
  list(category=names(ta),value=ta)
},.(UniProtType)]

DT[, fraction := value / sum(value)]
DT[, ymax := cumsum(fraction)]
DT[, ymin := shift(ymax, fill = 0)]
DT[, label_pos := (ymax + ymin) / 2]



DT$CleanCategory  <- factor(DT$category,unique(DT[,sum(value),category][order(V1,decreasing = T)]$category))

# ggplot(DT,aes(category,value,fill=UniProtType))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle=90))#+facet_grid("UniProtType",scale="free_y")+xlab("Count")
ggplot(DT,aes(CleanCategory,value,fill=UniProtType))+geom_bar(stat="identity")+theme(axis.text.x = element_text(angle=90))+ylab("Count")#+facet_grid("UniProtType",scale="free_y")+xlab("Count")

#per Gene Information
GeneSummary<- Gdt[,{
  list(
    SequenceCount=length(unique(Sequence)),
    CleanAccessionCounts=length(unique(V3)),
    AccessionCounts=length(unique(Clean_Accession)),
    Categories=length(unique(CleanCategory)),
    AS=sum(CleanCategory=="AS"),
    Proteolytic_Maturation=sum(CleanCategory=="Chain"),
    Alternative_Translation=sum(CleanCategory=="Alternative Translation"),
    UnAssigned_Variant=sum(CleanCategory=="Unassigned Variant"),
    Methionine=sum(CleanCategory=="Methionine"),
    No_Variant=sum(CleanCategory=="No Variant")
    
  )
},Gene_2]


Counts <- GeneSummary[,(apply(.SD,2,function(x){table(x>0)})),.SDcols=SelVal <- c("AS", "Proteolytic_Maturation", "Alternative_Translation", "UnAssigned_Variant" ,"No_Variant")]
Counts <- data.frame(Counts)
Counts$Boolean <- rownames(Counts)
Counts <- data.table(Counts)
# CountsMelt <- melt(Counts,measure.vars = paste(SelVal,"N",sep="."),id.vars = c("Boolean"))
CountsMelt <- melt(Counts,measure.vars = SelVal,id.vars = c("Boolean"))

CountsMelt[,value_rel:=round(value/sum(value)*100),variable]
ggplot(CountsMelt,aes(variable,value,fill=Boolean,label=value_rel))+geom_bar(stat="identity",position="stack")+geom_text(position = position_stack(vjust = 0.5),color="white")+ggtitle("Protein Product Categories per Gene")+ylab("Gene Count")+xlab("")+
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))


length(unique(Gdt$Sequence))
length(unique(Gdt$Gene_2))

Gdt[Gene_2=="TP53"]
Gdt[Gene_2=="CEBPA"]
Gdt[Gene_2=="CEBPB"]

Gdt[Gene_2=="GNAS"]
Gdt[Gene_2=="ATF4"]
Gdt[Gene_2=="VEGFA"]
Gdt[Gene_2=="POLG"]
Gdt[Gene_2=="EIF4G2"]
Gdt[Gene_2=="GNAS"]


Gdt[,SeqPerGene:=length(unique(MSequence)),.(Gene_2,grepl("_",V1))]
Gdt[SeqPerGene==1&!grepl("_",V1),Type_2:="No Variant"]

gff <- copy(Gdt)


gff[,c("MS_length_referenceproteome","MS_length"):={
temp <- .SD
MS <- temp[grep("_",V2,invert=T),,]$MSequence
  list(length(unique(MS)),length(unique(MSequence)))
},Gene_2]
gff[CleanCategory=="AS"&MS_length_referenceproteome==1]
gff[CleanCategory=="AS"&MS_length==1]

table(gff$MS_length_referenceproteome==1)
table(gff$MS_length==1)

table(gff$MS_length==1)


ggplot(Gdt,aes(CleanCategory,log10(MaxSeqN/SeqN)))+geom_boxplot()+theme(axis.text.x = element_text(angle=45,hjust = 1,vjust=1))
# Extra, 
## if number of Sequences are one per gene in combined fasta and gff files:
gff[MS_length==1]$CleanCategory <- "No Variant"
## if number of Sequences in ReferenceProteome alone is one, but category is AS -> "No VAriant" since it could arise from gff derived ambigous entries.

gff[MS_length_referenceproteome==1&CleanCategory=="AS"]$CleanCategory <- "No Variant"



gff2 <- copy(gff)

gff2[,ReferenceProteomeAccession:=any(!grepl("_",V1)),Clean_Accession]
table(gff2$ReferenceProteomeAccession)

# gff[MS_length_referenceproteome==1&CleanCategory=="AS"]$CleanCategory <- "No Variant"

save(gff,file="compiled_gff_20260224_UniProt_version_2022.rda")
table(Gdt$Type_2)
stop()
ggplot(DT, aes(x="", y=value, fill=category)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(round(value/1000,1),"k")), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL, fill = NULL)+facet_grid("UniProtType")+  # outside labels
  geom_text_repel(
    aes(x = 1.4,
        y = label_pos,
        label = paste0(category, "")),
    size = 4
  ) +labs(x = NULL, y = NULL, fill = NULL) #+
  
  #theme_void() 

