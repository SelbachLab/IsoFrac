#Viewer5
# functions ---------
DF_numericClass <- function(x){
  hum <- sapply(x,function(x){
    x1 <- sum(is.na(x))
    tempx <- as.numeric(x);
    # if(sum(is.na(tempx))>x1){
    #   tempx <- x
    # }
    # tempx
    sum(is.na(tempx))>x1
  })
  # df <- as.data.frame(hum)
  df <- as.data.frame(x)
  for(i in which(!hum)){
    class(df[,i]) <- "numeric"
  }
  
  df
  
}

Traces_GGplot <- function(input,ens="ENSG00000120802",mpwcutoff=0,converttoMW=T){
  hdb <- input$hdbscanTables[EnsG==ens]
  tempValues <- input$TMTSubset_final_singlenorm_zscore[hdb$id,]
  
  tempValues$cl <- hdb$cl_2
  tempValues$id <- hdb$id
  if(any(grepl("^V[0-9]",colnames(tempValues)[1]))){
    me <- melt(tempValues,measure.vars = paste("V",1:28,sep = ""),id.vars = c("cl","id"))
    me$variable <- as.numeric(gsub("V","",me$variable))
  }else{
    me <- melt(tempValues,measure.vars = 2:29,id.vars = c("cl","id"))
    
  }
  
  if(length(input$MW_Predictor)>0&converttoMW){
    md <- input$MW_Predictor
    mds <- 10^predict(md,2:29)
    me$variable <- mds[match(me$variable,2:29)]
  }
  
  g2 <- ggplot(me,aes(as.numeric(variable),value,group=id,color = as.character(cl)))+geom_line()
  g2
}
# Traces_Plotly <- function(input,ens="ENSG00000120802",mpwcutoff=0,converttoMW=T,cbPalette= c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),PurePlot =F,showlegend=T){
#   
#   ens <<- ens
#   hdb <- input$hdbscanTables[EnsG==ens]
#   tempValues <- input$TMTSubset_final_singlenorm_zscore[hdb$id,]
#   
#   tempValues$cl <- hdb$cl_2
#   tempValues$id <- hdb$id
#   if(any(grepl("^V[0-9]",colnames(tempValues)[1]))){
#     me <- melt(tempValues,measure.vars = paste("V",1:28,sep = ""),id.vars = c("cl","id"))
#     me$variable <- as.numeric(gsub("V","",me$variable))
#   }else{
#     
#     me <- melt(tempValues,measure.vars = as.character(2:29),id.vars = c("cl","id"))
#     me$variable <- as.numeric(as.character(me$variable))
#     me$variable <- me$variable-min(me$variable)+1
#     input <<- input
#     me <<- me
#     # plot(me$variable2,as.numeric(as.character()))
#     # abline(0,1)
#     (table(me$variable))
#     # stop()
#   }
#   
#   
#   if(length(input$Fractions_predicted_MW)>0&converttoMW){
#     mds <- input$Fractions_predicted_MW[1:28]
#     # mds <- 10^predict(md,2:29)
#     # mds <- 10^predict(md,1:28)
#     # table(me$variable)
#     # (table(me$variable))
#     me$variable <- 10^mds[match(me$variable,1:28)]
#     # me <- me[order(variable)]
#   }
#   
#   
#   
#   # med <-  dcast(me,variable~paste(cl,id,sep ="_"),value.var = "value")
#   # fig <- plot_ly(med,x=~variable,y=~names(med)[2],type = 'scatter', mode = 'lines')
#   # for(i in unique(me$cl)){
#   #   add <- grep(paste("^",i,"_",sep = ""),names(med),value=T)
#   #   for(u in add){
#   #     fig <- fig%>%add_trace(y=add)
#   #   }
#   # }
#   me <- rbind(me,me[,.(value=max(me$value)+1,variable=NA),.(id,cl)])
#   
#   me$cl <- (factor(me$cl))
#   me$seq <- input$dt_gn$seq[me$id]
#   me <- me[,.SD[order(variable)],id]
#   me <- me[order(id),]
#   me[,idtemp:={.GRP},.(cl,id)]
#   try({me <- me[order(idtemp),]})
#   try({me <- me[order(seq),]})
#   try({me <- me[as.numeric(order(cl)),]})
#   try({me <- me[order(id),]})
#   
#   # me <- me[,.SD[order(as.numeric(id))],id]
#   
#   idCount <- 0
#   idRemove <- c()IFUN
#   for(i in unique(me$cl)){
#     id <- me[cl==i]$id
#     id <- unique(id)
#     id <- (1:length(id))+idCount
#     idCount <- max(id)
#     idRemove <- c(id[-1])
#     
#     
#   }
#   cbp <- colorRampPalette(cbPalette)(max(as.character(me$cl)))
#   cbp <- c("#73FDFF",cbp)
#   if(PurePlot){
#     print("Hu")
#     cbp <- "#00000080"
#     me$color <- "black"
#     me$name <- ""
#     me$legendgroup <- ""
#     # me$cl <- ""
#   }else{
#     me$color <- cbp[as.numeric(as.character(me$cl))+1]
#     
#   }
#   # me <- me[color=="#E69F00"]
#   me <<- me
#   p <- plot_ly(me, 
#                x = ~variable, 
#                y = ~value, 
#                # id = ~id,
#                color=~cl,
#                name=~cl,
#                legendgroup=~cl,
#                type = 'scatter',
#                mode = 'lines',
#                hoverinfo="text",
#                text=me$seq,
#                colors=cbp)%>%
#     layout(showlegend = TRUE) %>%
#     style(showlegend = FALSE, traces = idRemove)
#   # if(PurePlot){
#   #   p <- p%>%layout(title = ens, xaxis = list(title = 'MW converted fraction'), 
#   #                   yaxis = list(title = "zscore"),showlegend=showlegend)
#   # }else{
#   # }
#   p <- p%>%layout(title = ens, xaxis = list(title = 'MW converted fraction'), 
#                   yaxis = list(title = "zscore"), legend = list(title=list(text='hdbscan results')),showlegend=showlegend)
#   
#   
#   list(p=p,data=me,clusterColors=cbp) 
# }
Traces_Plotly <- function(input_data,mpwcutoff=0,converttoMW=T,cbPalette= c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),PurePlot =F,showlegend=T,MeanCluster=F,ShowQuantiles=T,addPolygon=T){
  
  #mpwcutoff=0,converttoMW=T,cbPalette= c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),PurePlot =F,showlegend=T,MeanCluster=F,ShowQuantiles=T,addPolygon=T
  
  input_data <<- input_data
  ens <- paste(unique(input_data$dt_gn$ens, input_data$dt_gn$hgnc_symbol),collapse="")
  hdb <- input_data$hdbscanTables
  tempValues <- input_data$inputtable
  hdb <- hdb[match(input_data$ID,hdb$id),]
  
  tempValues$cl <- as.character(hdb$cl_2)
  vec <- c(0,sort(unique(tempValues$cl[tempValues$cl!=0])))
  tempValues$cl <- factor(tempValues$cl,vec,c("noise traces",vec[-1]))
  tempValues$id <- hdb$id
  
  
  
  if(any(grepl("^V[0-9]",colnames(tempValues)[1]))){
    me <- melt(tempValues,measure.vars = paste("V",1:28,sep = ""),id.vars = c("cl","id"))
    me$variable <- as.numeric(gsub("V","",me$variable))
  }else{
    colnames(tempValues) <- gsub("X","",colnames(tempValues))
    me <- melt(tempValues,measure.vars = as.character(2:29),id.vars = c("cl","id"))
    me$variable <- as.numeric(as.character(me$variable))
    me$variable <- me$variable-min(me$variable)
    # me <<- me
    # plot(me$variable2,as.numeric(as.character()))
    # abline(0,1)
    # length(table(me$variable))
    # stop()
  }
  # me <<- me
  # stop()
  
  me <- me[!is.na(variable)]
  
  if(length(input_data$Fractions_MW$x.Fractions_predicted_MW)>0&converttoMW){
    mds <- input_data$Fractions_MW$x.Fractions_predicted_MW[1:28]
    # mds <- 10^predict(md,2:29)
    # mds <- 10^predict(md,1:28)
    # table(me$variable)
    # (table(me$variable))
    me$variable <- 10^mds[match(me$variable,1:28)]
    # me <- me[order(variable)]
  }
  
  
  # med <-  dcast(me,variable~paste(cl,id,sep ="_"),value.var = "value")
  # fig <- plot_ly(med,x=~variable,y=~names(med)[2],type = 'scatter', mode = 'lines')
  # for(i in unique(me$cl)){
  #   add <- grep(paste("^",i,"_",sep = ""),names(med),value=T)
  #   for(u in add){
  #     fig <- fig%>%add_trace(y=add)
  #   }
  # }
  me <- rbind(me,me[,.(value=max(me$value)+1,variable=NA),.(id,cl)])
  
  me$cl <- (factor(me$cl))
  me$seq <- input_data$dt_gn$seq[match(me$id,input_data$ID)]
  
  me <- me[,.SD[order(variable)],id]
  me <- me[order(id),]
  me[,idtemp:={.GRP},.(cl,id)]
  try({me <- me[order(as.numeric(idtemp)),]})
  try({me <- me[order(as.character(seq)),]})
  try({me <- me[as.numeric(order(as.numeric(cl))),]})
  try({me <- me[order(as.numeric(id)),]})
  
  # me <- me[,.SD[order(as.numeric(id))],id]
  
  idCount <- 0
  idRemove <- c()
  for(i in unique(me$cl)){
    id <- me[cl==i]$id
    id <- unique(id)
    id <- (1:length(id))+idCount
    idCount <- max(id)
    idRemove <- c(id[-1])
    
    
  }
  
  if(MeanCluster){
    me <- me[,{list(value=mean(value,na.rm = T),q90=quantile(value,0.9,na.rm = T),q10=quantile(value,0.1,na.rm = T),n=length(value),seq=paste(unique(seq),collapse="\n"),idtemp=.GRP,id=.GRP)},.(cl,variable)]
  }
  
  
  cbp <- colorRampPalette(cbPalette)(max(as.numeric(me$cl),na.rm=T))
  cbp <- c("#73FDFF",cbp)
  
  if(PurePlot){
    cbp <- "#00000080"
    me$color <- "black"
    me$name <- ""
    me$legendgroup <- ""
    # me$cl <- ""
  }else{
    me$color <- cbp[as.numeric(as.character(me$cl))+1]
    
  }
  
  me <- me
  me$q90a <- me$q90
  me$q10a <- me$q10
  if(ShowQuantiles&MeanCluster==T){
    me$q90 <- abs(me$q90-me$value)
    me$q10 <- abs(me$value-me$q10)
  }else{
    me$q90 <-NA
    me$q10 <-NA
    
  }
  
  
  me <- me[!is.na(me$cl),]
  # me <- me[]
  # me <- me8
  p <- plot_ly(me, 
               x = ~variable, 
               y = ~value, 
               # id = ~id,
               color=~cl,
               name=~cl,
               legendgroup=~cl,
               type = 'scatter',
               mode = 'lines',
               hoverinfo="text",
               text=me$seq,
               colors=cbp,
               line=list(width=3))%>%
    layout(showlegend = TRUE) 
  
  p <- try(p%>%
             style(showlegend = FALSE, traces = idRemove))
  # print("hu")
  if(addPolygon&MeanCluster){
    for(i in unique(me$cl)){
      print(i)
      mei <- me[cl==i]
      
      
      mei <- mei[order(variable,decreasing = T)]
      mei <- mei[!is.na(mei$variable)]
      mei$color[is.na(mei$color)] <- "#73FDFF"
      
      
      p <- p %>% add_polygons(
        x=c(mei$variable,rev(mei$variable),mei$variable[1]),
        y=c(mei$q90a,rev(mei$q10a),mei$q90a[1]),
        line=list(width=2,color="transparent"),
        fillcolor=paste(mei$color[1],20,sep=""), inherit = FALSE,legendgroup=mei$cl[1],name="0.1 to 0.9 Q")
    }
    
  }
  # print("hu")
  
  p <- p%>%layout(title = ens, xaxis = list(title = 'MW converted fraction'), 
                  yaxis = list(title = "zscore"), legend = list(title=list(text='Peptide Clusters')),showlegend=showlegend)
  
  
  list(p=p,data=me,clusterColors=cbp) 
}

# g1 <- Traces_Plotly(dapi,converttoMW = T,MeanCluster = input$MeanCluster,addPolygon  = input$ShowQuantiles)

# g1 <- Traces_Plotly(dapi,converttoMW = T,MeanCluster = T,addPolygon  = T)
# g1
Traces_Data <- function(input,ens="ENSG00000120802",mpwcutoff=0,converttoMW=T){
  hdb <- input$hdbscanTables[EnsG==ens]
  tempValues <- input$TMTSubset_final_singlenorm_zscore[hdb$id,]
  
  tempValues$cl <- hdb$cl_2
  tempValues$id <- hdb$id
  if(any(grepl("^V[0-9]",colnames(tempValues)[1]))){
    me <- melt(tempValues,measure.vars = paste("V",1:28,sep = ""),id.vars = c("cl","id"))
    me$variable <- as.numeric(gsub("V","",me$variable))
  }else{
    me <- melt(tempValues,measure.vars = 2:29,id.vars = c("cl","id"))
    
  }
  
  if(length(input$MW_Predictor)>0&converttoMW){
    md <- input$MW_Predictor
    mds <- 10^predict(md,2:29)
    me$variable <- mds[match(me$variable,2:29)]
  }
  
  
  me
}

ggtraces <- function(input,ens_identifier="ENSG00000120802"){
  
  
  PeakResults_TMPO<- data.table(input$IsoFrac_PeakResultsTable)[ENS==ens_identifier]
  
  PeakResults_TMPO$id <- PeakResults_TMPO$cluster2
  PeakResults_TMPO$cl <- PeakResults_TMPO$cluster2
  
  g2 <- Traces_GGplot(input)
  g3 <- NULL
  try({
    g3 <- g2+geom_errorbar(data=PeakResults_TMPO,aes(x=10^unlist(MW_lm),y=unlist(Height_lm),ymin=unlist(Height_BS_lo),ymax=unlist(Height_BS_hi)),color="black")+
      geom_errorbarh(data=PeakResults_TMPO,aes(x=10^unlist(MW_lm),y=unlist(Height_lm),xmin=10^unlist(MW_lm_lo),xmax=10^unlist(MW_lm_hi)),color="black")+xlab("log10 MW")+ylab("cluster zscore [au]")+
      geom_label_repel(data=PeakResults_TMPO,aes(x=10^unlist(MW_lm),y=unlist(Height_lm),label=round(10^unlist(MW_lm))),size = 2)+
      geom_point(data=PeakResults_TMPO,aes(10^unlist(MW_lm),unlist(Height_lm),label=round(10^unlist(MW_lm)),group=NULL,color=as.character(cluster2)),color="black")
  })
  
  list(Traces=g2,TraceAndPoints=g3)
}

extract_Info <- function(ih,gfu,ai,peaksfun,da,mw_info){
  # detecting Sequence near selected point
  np <- nearPoints(da, ih, maxpoints=1)
  if(length(np)==0){
    return(NULL)
  }
  # return(NULL)
  id <- np$ID-10
  print("Extracting Cluster")
  # extracting other sequences within the cluster
  otherInCluster <- da$ID[da$Color==np$Color]
  otherInCluster <- unique(otherInCluster)-10
  print("Extracting Peaks")
  
  # Peaks of Cluster:
  peaksfun <- peaksfun
  np <- np
  mw_info <- mw_info
  clusterPeak <- peaksfun[sapply(strsplit(peaksfun$co,";"),function(x){any(x==np$Color)}),]$V2
  OBSERVED_MW<- sapply(clusterPeak,FractionToMW,mw_info)#äFractionToMW(clusterPeak,mw_info)
  # OBSERVED_MW <- "test"
  
  SelectedSequence <- ai$SequenceInfo$PeptideSequence[id]
  OtherSequenceInCluster<- ai$SequenceInfo$PeptideSequence[otherInCluster]
  print("Protein Sequence")
  
  dt_ps <- data.table({ai$SequenceInfo$ProteinSequence})
  dt_ps <- dt_ps[,as.list(apply(.SD,2,function(x){paste(unique(x),collapse=";")})),SEQUENCE]
  print("Protein Sequence Match")
  
  SequenceMatch <- lapply(OtherSequenceInCluster,function(x){
    r1 <- regexpr(x,dt_ps$SEQUENCE)
    r1
    
  })
  print("List Output")
  pr <- "not implemented"
  
  list(pe_selected=SelectedSequence,pr=pr,
       cluster_se=OtherSequenceInCluster,
       SequenceMatch=SequenceMatch,
       ProteinInfo= dt_ps,MW=OBSERVED_MW,
       PeaksInfo = peaksfun,
       data=da,
       mw=mw_info)
}

# plot.sequence <- function(IFUN,selected=IFUN$ENSEMBL_PROTEIN[1],ClustColor=NULL){
#   ClustColor <<- ClustColor
#   IFUN <<- IFUN
#   selected <<- selected
#   it <- 1
#   PS <- IFUN[ENSEMBL_PROTEIN==selected,paste(unique(ENSEMBL_PROTEIN),collapse=";"),ProteinSequence]
#   Positions <- IFUN[ENSEMBL_PROTEIN==selected,{
#     gr <- gregexpr(Sequence,ProteinSequence)
#     
#     Posis <- gr[[1]]:(gr[[1]]+attributes(gr[[1]])$match.length-1)
#     Col <- 0
#     if(length(ClustColor)>0){
#       Col <- ClustColor$Color[match(Cluster,ClustColor$Cluster)]
#     }
#     Col <- unique(Col)
#     # if(length(Col)>1){
#     #   Col <- "#888888"
#     # }
#     data.table(Position=Posis,Color=Col,Cluster=Cluster)
#   },.(ProteinSequence,Sequence)]
#   
#   Positions <- Positions[,.(co=paste(unique(Color),collapse = ";"),cl=paste(unique(Cluster),collapse = ";")),Position]
#   
#   Positions <- Positions[Position>0]
#   Positions <- Positions[!is.na(Position)]
#   
#   LE <- strsplit(PS$ProteinSequence,"")
#   maxLe <- 50
#   LEs <- lengths(LE)
#   
#   
#   yrow <- ceiling(LEs/maxLe)[it]
#   maxLe <- 50
#   vec <<- rep(0,(yrow*maxLe))
#   # vec[1:10] <- 1
#   Positions <<- Positions
#   vec[Positions$Position] <- (Positions$cl)
#   # vec[is.na(vec)] <- max(vec,na.rm = T)+1
#   vec <- as.numeric(factor(as.character(vec),ClustColor$Cluster,1:length(ClustColor$Cluster)))
#   vec[is.na(vec)] <- 0
#   # vec[is.na(vec)] <- max(vec,na.rm = T)+1
#   MA <- matrix(rev(vec),nrow = yrow  ,ncol = maxLe,byrow = T)
#   MA <- apply(MA,1,rev)
#   co <- unique(Positions$co)
#   co <- co[grep(";",co,invert =T)]
#   par(mai=c(0.1,0.1,0.5,1.2))
#   image((MA),main="",col=c("#FFFFFF",ClustColor$Color),axes=F)
#   mtext(selected,3,line = 1)
#   legend(par()$usr[2],par()$usr[4],xjust = 0,legend = c(as.character(ClustColor$Cluster),">x"),fill=c(ClustColor$Color,"#33333340"),xpd=NA,bty="n")
#   
#   y <<- 0
#   x <<- 0
#   sapply(1:length(LE[[it]]),function(ite){
#     cat("\r",x)
#     
#     
#     text(x/(maxLe-1),1-y/(yrow-1),LE[[it]][ite],xpd=NA)
#     
#     x <<- x+1
#     if(x%%maxLe==0){
#       
#       y<<-y+1
#       x <<- 0
#     }
#     
#   })
# }
plot.sequence <- function(IFUN,selected=IFUN$ENSEMBL_PROTEIN[1],ClustColor=NULL,ProtSequence=NULL){
  ClustColor <<- ClustColor
  IFUN <<- IFUN
  selected <<- selected
  ProtSequence <<- ProtSequence
  it <- 1
  IFUN$Sequence <- IFUN$Seq
  IFUN$Cluster <- IFUN$Cl
  # tempi <- sapply(IFUN$Seq,function(x){
  #   grep(x,UniSequences_all$SEQUENCE,fixed = T)
  # })
  ClustColorMatch <- ClustColor[match(IFUN$Cluster,ClustColor$cluster2),]
  IFUN$Color <- ClustColorMatch$Color
  IFUN$Peaks <- ClustColorMatch$V1
  
  Positions <- IFUN[,{
    cat("\r",.GRP)
    sequ <- .BY$Sequence
    gr <- gregexpr(sequ,ProtSequence,fixed = T)
    
    Posis <- gr[[1]]:(gr[[1]]+attributes(gr[[1]])$match.length-1)
    
    Col <- 0
    if(length(ClustColor)>0){
      Col <- ClustColor$Color[match(Cluster,ClustColor$Cluster)]
    }
    Col <- unique(Col)
    # if(length(Col)>1){
    #   Col <- "#888888"
    # }
    data.table(Position=Posis)
  },.(Sequence,Color,Peaks)]
  table(Positions$Position)
  Positions$Color[is.na(Positions$Color)] <-   Positions$Peaks[is.na(Positions$Color)]
  
  # Positions <- Positions[Position>0,]
  
  Positions <- Positions[,.(co=paste(unique(Color),collapse = ";"),cl=paste(unique(Peaks),collapse = ";")),Position]
  
  Positions <- Positions[Position>0]
  Positions <- Positions[!is.na(Position)]
  
  LE <- strsplit(ProtSequence,"")
  maxLe <- 50
  LEs <- lengths(LE)
  
  
  yrow <- ceiling(LEs/maxLe)[it]
  maxLe <- 50
  vec <<- rep(0,(yrow*maxLe))
  # vec[1:10] <- 1
  Positions <<- Positions
  vec[Positions$Position] <- (Positions$cl)
  # vec[is.na(vec)] <- max(vec,na.rm = T)+1
  vec <- as.numeric(factor(as.character(vec),ClustColor$V1,1:length(ClustColor$V1)))
  vec[is.na(vec)] <- 0
  # vec[is.na(vec)] <- max(vec,na.rm = T)+1
  MA <- matrix(rev(vec),nrow = yrow  ,ncol = maxLe,byrow = T)
  MA <- apply(MA,1,rev)
  co <- unique(Positions$co)
  co <- co[grep(";",co,invert =T)]
  par(mai=c(0.8,0.1,0.5,2))
  image((MA),main="",col=c("#FFFFFF",ClustColor$Color),axes=F)
  mtext(selected,3,line = 1)
  ClustColorLeg <- ClustColor
  ClustColorLeg$cluster2 <- NULL
  ClustColorLeg <- unique(ClustColorLeg)
  legend(par()$usr[2],par()$usr[4],xjust = 0,legend = c(as.character(ClustColorLeg$V1),">x"),fill=c(ClustColorLeg$Color,"#33333340"),xpd=NA,bty="n",cex = 0.8)
  legend(par()$usr[1],par()$usr[3],xjust = 0,legend = c(as.character(ClustColorLeg$V1),">x"),fill=c(ClustColorLeg$Color,"#33333340"),horiz = T,xpd=NA,bty="n")
  
  y <<- 0
  x <<- 0
  sapply(1:length(LE[[it]]),function(ite){
    cat("\r",x)
    
    
    text(x/(maxLe-1),1-y/(yrow-1),LE[[it]][ite],xpd=NA)
    
    x <<- x+1
    if(x%%maxLe==0){
      
      y<<-y+1
      x <<- 0
    }
    
  })
}
plotly.sequence.3d <- function(IFUN,selected=IFUN$ENSEMBL_PROTEIN[1],ClustColor=NULL,ProtSequence=NULL,plot3d=T,SequenceWidth=20){
  ClustColor <<- ClustColor
  idRemove <- 0
  # ClustColor$V1 <- paste(ClustColor$V1,"Cluster:",ClustColor$cluster2)
  IFUN <<- IFUN
  selected <<- selected
  ProtSequence <<- ProtSequence
  
  IFUN$Sequence <- IFUN$Seq
  IFUN$Cluster <- IFUN$Cl
  # tempi <- sapply(IFUN$Seq,function(x){
  #   grep(x,UniSequences_all$SEQUENCE,fixed = T)
  # })
  ClustColorMatch <- ClustColor[match(IFUN$Cluster,ClustColor$cluster2),]
  IFUN$Color <- ClustColorMatch$Color
  IFUN$Peaks <- ClustColorMatch$V1
  PositionsList <- lapply(unique(ProtSequence$ProteinID),function(Pid){
    ProtSequence_i <<- ProtSequence[ProteinID==Pid]$SEQUENCE
    # ProtSequence_i <- ProtSequence[ProteinID==Pid]$SEQUENCEBlast
    
    Positions <- IFUN[,{
      cat("\r",.GRP)
      sequ <<- .BY$Sequence
      gr <- gregexpr(sequ,ProtSequence_i,fixed = T)
      # stop()
      
      Posis <- gr[[1]]:(gr[[1]]+attributes(gr[[1]])$match.length-1)
      
      Col <- 0
      if(length(ClustColor)>0){
        Col <- ClustColor$Color[match(Cluster,ClustColor$Cluster)]
      }
      Col <- unique(Col)
      # if(length(Col)>1){
      #   Col <- "#888888"
      # }
      data.table(Position=Posis)
    },.(Sequence,Color,Peaks)]
    
    
    
    Positions$ProteinID <- Pid
    # Positions$MaxRow <- MaxRow
    Positions
  })
  Positions <<- rbindlist(PositionsList)
  ProtSequence <<- ProtSequence
  
  # Alignment
  try({
    afu <- align_sequences(ProtSequence$SEQUENCE)
  })
  ProtSequence <- ProtSequence[hclust(dist(t(sapply(afu,function(x){return(x)}))))$order]
  
  Positions <- Positions[order(ProteinID),{.SD[order(as.numeric(Position))]},ProteinID]
  Positions[,c("Group","AlignPosition"):={
    afutemp <<- afu[.BY$ProteinID==ProtSequence$ProteinID][[1]]
    temppo <<- .SD
    posi <- sort(unique(temppo$Position))
    M <- match(posi,afutemp)
    # Msort <- sort(M)
    
    
    groupvec <- rep(0,length(M))
    it <- 1
    newgroup <- F
    dMsort <- diff(M)>1
    dMsort[is.na(dMsort)] <- F
    for(i in 1:length(dMsort)){
      if(dMsort[i]){
        
        newgroup <- T
      }else{
        if(newgroup){
          it <- it+1
          newgroup <- F
        }
        groupvec[i] <- it
      }
      
    }
    # print(length(groupvec))
    # print(length(M))
    po <- match(temppo$Position,posi)
    list(  groupvec[po],M[po])
  },.(ProteinID)]
  # Positions <- Positions[grepl("P42167-3",ProteinID,fixed = T)]
  
  if(plot3d){
    
    library(plotly)
    
    
    fig <- plot_ly() 
    # Positions <- Positions[Positions$Position >0]
    dfall <<- Positions[,{as.list(range(AlignPosition))},.(Sequence,Color,Peaks,ProteinID,Group)]
    dfall <- dfall[order(V1)]
    dfall <- dfall[Group!=0]
    dfall <- dfall[!is.na(Peaks)]
    dfall[,AlignPosition:={
      afutemp <<- afu[.BY$ProteinID==ProtSequence$ProteinID][[1]]
      temppo <<- .SD$Position
      match(temppo,afutemp)
      
    },ProteinID]
    
    dfall$Color <- rainbow(max(dfall$Color,na.rm=T))[dfall$Color]
    # MaxRow <- 1
    MaxRowCollect <- c()
    ProteinOrder <- unique(dfall$ProteinID)
    idCount <- 0
    idRemove <- c()
    PeaksUnique <- unique(dfall$Peaks)
    PeaksUnique <- c(PeaksUnique[!is.na(PeaksUnique)],"THISDOESNTEXIST")
    
    idRemove <- c()
    idColor <- c()
    
    for(pid in ProteinOrder){
      currentRow <- 0
      MaxRow <- 1
      ({
        df <- dfall[ProteinID==pid,]
        # plotly(dfall,y=V1~V2,x=~)
        
        for(i in 1:(nrow(df))){
          
          try({
            if(df$V1[i]<=df$V2[i-1]){
              currentRow <- currentRow+1
            }
            if(df$V1[i]>=df$V2[i-1] & currentRow>0){
              currentRow <- 0
            }
          },silent=T)
          MaxRow <- max(MaxRow,currentRow)
          
          
          
          
          fig <- add_trace(fig,
                           y = c(df$V1[i], df$V2[i]),  # x0, x1
                           z = c(currentRow, currentRow),  # y0, y1
                           x = c(df$ProteinID[i], df$ProteinID[i]),
                           mode = "lines",
                           line = list(color = df$Color[i], width = SequenceWidth),
                           showlegend = F,
                           hoverinfo = "text",
                           type="scatter3d",
                           name=df$Sequence[i],
                           legendgroup=df$Sequence[i],
                           # id=df$Peaks[i],
                           # type="scatter",
                           
                           # Create custom hover text
                           
                           text = paste("ProteinID:",paste(df$ProteinID[i]), "<br>",
                                        "Sequence: ", df$Sequence[i], "<br>",
                                        "Peak: ", df$Peaks[i], "<br>",
                                        "Position: ", paste(df$V1[i],df$V2[i]))
                           
                          
                           # evaluate = T  # needed to avoid lazy loading
          )
        }
        if(MaxRow>5){
          MaxRow <- MaxRow
        }else{
          MaxRow <- 5
        }
        MaxRowCollect <- c(MaxRowCollect,MaxRow)
      })
      
    }
    MaxRow <- max(MaxRowCollect)
    
    # add sequences
    currentRow <- -1
    # afu <<- afu
    # figtest <<- fig
    for(i in 1:length(ProteinOrder)){
      
      a <- afu[[which(ProtSequence$ProteinID==ProteinOrder[i])]]
      # a[a==0] <- NA
      # print(a)
      atemp <- 1:length(a)
      atemp[a==0] <-NA
      a <- atemp
      # stop()
      # a <- as.numeric(a!=0)
      fig <- add_trace(fig,
                       x = rep(ProteinOrder[i],length(a)),
                       y = a,  # x0, x1
                       z = rep(currentRow,length(a)),  # y0, y1
                       mode = "lines",
                       type="scatter3d",
                       line = list(color = "grey", width = SequenceWidth),
                       showlegend = F,
                       legendgroup="NONE",
                       hoverinfo = "text",
                       name="Whole Protein",
                       # type="scatter",
                       
                       # Create custom hover text
                       
                       text = ProtSequence$ProteinID[i]
                       
                       # evaluate = T  # needed to avoid lazy loading
      )
      
      # fig
    }
    
    fig$idRemove <- idRemove
    
    # print("hu")
    # print(MaxRow)
    # print("Hu2")
    # fig <- fig%>%
    #   layout(xaxis = list(range = c(1,nchar(ProtSequence$ProteinID)),title="Protein Sequence Position"),
    #          yaxis = list(range = c(-1,MaxRow+10),title="Sequence Overlap"),title=selected)
    
    # Adding axes names
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'Protein ID',visible=F),
                                       yaxis = list(title = 'Aligned Sequence Position'),
                                       zaxis = list(title = 'Overlapping Peptide Layers',range=c(-1,MaxRow+5))))
    # adding legend by selecting the first entry in the attributes of the legendgroups
    at <- fig$x$attrs
    lg <- sapply(at,function(x){x$legendgroup})
    lg[lengths(lg)==0] <- "init"
    # lg <- lg[lengths(lg)!=0]
    lg <- as.vector(lg)
    # lg <- as.vector(lg)
    
    idremove <- setdiff(1:length(lg),match(unlist(unique(lg)),lg))
    fig %>%
      style(showlegend = FALSE) %>%
      style(showlegend = T, traces = match(unique(unlist(lg)),lg))
    
  }else{
    hm <- sapply(1:length(afu),function(i){
      ProtSequence$ProteinID
      ai <- afu[[i]]
      ai[ai==0] <- -1
      Posis <- Positions[ProteinID==ProtSequence$ProteinID[i]&Position>0]
      # Posis$
      PosisTemp <- Posis[match(ai,Posis$Position)]$Color[ai>=0]
      PosisTemp[is.na(PosisTemp)] <- 0
      ai[ai>=0] <- PosisTemp
      ai
      
      
    })
    colnames(hm) <-  ProtSequence$ProteinID
    
    hmo <- hm[,hclust(dist(t(hm)))$order]
    hmfun <- unique(as.vector(hmo))
    l <- length(hmfun[hmfun>0])
    minus <- any(hmfun==-1)
    ZERO <- any(hmfun==0)
    colors <- c("#FFFFFF"[minus],"#11111140"[ZERO],rainbow(max(hmfun),alpha=0.99)[sort(hmfun[hmfun>0])])
    # plot(1:length(colors),col=colors,cex = 20)
    
    Ticknames <- paste(ClustColor$V1,"Cl:",ClustColor$cluster2)[match(hmfun[hmfun!=0&hmfun!=-1],ClustColor$Color)]
    if(any(hmfun==0)){
      Ticknames <- c("Whole Protein",Ticknames)
    }
    if(any(hmfun==-1)){
      Ticknames <- c("Gap",Ticknames)
    }
    # basic heatmap
    # hmo <- apply(hmo,2,as.character)
    # ProtSequence
    df <- data.frame(z=seq(0,1,length.out=length(colors)),color=colors)
    df$color <- as.character(df$color)
    df$z <- as.character(df$z)
    hmo <- t(hmo)
    Ticknames <<- Ticknames
    colorbar=list(tickmode='array', tickvals=sort(hmfun), ticktext=Ticknames, len=0.2)
    class(hmo[,1])
    hmo <<- hmo 
    p <- plot_ly(x=colnames(hmo), y=rownames(hmo), z = hmo, 
                 xgap=2,
                 ygap=2,
                 type = "heatmap",
                 zmin=min(hmo),
                 zmax=max(hmo),
                 # split=hmo,
                 # mode="markers",
                 colorscale=df,
                 # hoverinfo = "x+y+z",
                 # colorbar=list(tickmode='array', tickvals=c(hmcols), ticktext=levels(hmo), len=0.2)
                 colorbar=colorbar
                 # clusters=length(hmfun)+2
    ) %>%
      layout(margin = list(l=120))
    p
  }
  
}


align_sequences <- function(se,se_specific=NULL,...){
  require(msa)
  if(1)
    # if(1)
    try({
      
      # se <<- se
      if(length(se_specific)>0){
        print("HU")
        
        lu <- lapply(se_specific,function(y){
          y <<- y
          hum <<- gregexec(y,se)
          ml <- sapply(hum,function(x){
            attributes(x)$match.length
          })
          lapply(1:length(hum),function(i){
            hui <- hum[[i]][[1]]
            hui:ml[i]
            
          })
          
        })
        se_new <- lapply(1:length(se),function(i){
          
          keepID <- unique(unlist(lapply(lu,function(x){
            x[[i]]
          })))
          keepID <- keepID[keepID!=-1]
          setemp <- strsplit(se[[i]],"")[[1]]
          seinit <- rep("X",length(setemp))
          seinit[keepID] <- setemp[keepID]
          seinit
        })
        se_new <- lapply(se_new,function(x){x <- x[x!="X"]})
        se_new <- sapply(se_new,paste,collapse="")
        
        se_new <- unlist(se_new)
        
        se_backup <- se
        se <- se_new
        
      }
      # dfAlign <- msa(se,type = "protein",verbose=F,gapOpening = 5,gapExtension = 0,order="input")
      # print(se)
      # se <<- se
      dfAlign <- msa(se,type = "protein",verbose=F,method = "ClustalOmega")
      # stop()
      # dfAlign
      # df
      # df@unmasked@ranges
      dff <- as.data.frame(dfAlign@unmasked)
      M <- matrix(NA,ncol=nchar(dff$x[1]),nrow=dim(dff)[1])
      
      MatchFun <- lapply(strsplit(dff$x,""),function(x){
        x <<- x
        sapply(1:length(x),function(y){
          if(x[y]!="-"){
            s <- sum(x[1:y]!="-")
          }else{
            s <- 0
          }
          
        })
      })
      if(length(names(se))>0){
        names(MatchFun) <- names(se)
      }
      MatchFun
    })
  
}

plotly.sequence <- function(IFUN,selected=IFUN$ENSEMBL_PROTEIN[1],ClustColor=NULL,ProtSequence=NULL){
  ClustColor <<- ClustColor
  IFUN <<- IFUN
  selected <<- selected
  ProtSequence <<- ProtSequence
  it <- 1
  IFUN$Sequence <- IFUN$Seq
  IFUN$Cluster <- IFUN$Cl
  # tempi <- sapply(IFUN$Seq,function(x){
  #   grep(x,UniSequences_all$SEQUENCE,fixed = T)
  # })
  ClustColorMatch <- ClustColor[match(IFUN$Cluster,ClustColor$cluster2),]
  IFUN$Color <- ClustColorMatch$Color
  IFUN$Peaks <- ClustColorMatch$V1
  
  Positions <- IFUN[,{
    cat("\r",.GRP)
    sequ <- .BY$Sequence
    gr <- gregexpr(sequ,ProtSequence,fixed = T)
    
    Posis <- gr[[1]]:(gr[[1]]+attributes(gr[[1]])$match.length-1)
    
    Col <- 0
    if(length(ClustColor)>0){
      Col <- ClustColor$Color[match(Cluster,ClustColor$Cluster)]
    }
    Col <- unique(Col)
    # if(length(Col)>1){
    #   Col <- "#888888"
    # }
    data.table(Position=Posis)
  },.(Sequence,Color,Peaks)]
  
  library(plotly)
  
  
  fig <- plot_ly() 
  Positions <- Positions[Positions$Position >0]
  df <- Positions[,{as.list(range(Position))},.(Sequence,Color,Peaks)]
  df <- df[order(V1)]
  
  df$Color <- rainbow(max(df$Color))[df$Color]
  currentRow <- 0
  MaxRow <- 1
  for(i in 1:(nrow(df))){
    try({
      if(df$V1[i]<=df$V2[i-1]){
        currentRow <- currentRow+1
      }
      if(df$V1[i]>=df$V2[i-1] & currentRow>0){
        currentRow <- 0
      }
    })
    MaxRow <- max(MaxRow,currentRow)
    SequenceWidth <- 10
    fig <- add_trace(fig,
                     x = c(df$V1[i], df$V2[i]),  # x0, x1
                     y = c(currentRow, currentRow),  # y0, y1
                     mode = "lines",
                     line = list(color = df$Color[i], width = SequenceWidth),
                     showlegend = F,
                     hoverinfo = "text",
                     type="scatter",
                     
                     # Create custom hover text
                     
                     text = paste("Sequence: ", df$Sequence[i], "<br>",
                                  "Peak: ", df$Peaks[i], "<br>",
                                  "Position: ", paste(df$V1[i],df$V2[i])),
                     
                     evaluate = T  # needed to avoid lazy loading
    )
  }
  if(MaxRow>5){
    MaxRow <- MaxRow
  }else{
    MaxRow <- 5
  }
  fig%>%
    layout(xaxis = list(range = c(1,nchar(ProtSequence)),title="Protein Sequence Position"),yaxis = list(range = c(-1,MaxRow),title="Sequence Overlap"),title=selected)
  
}
# plotly.sequence(IFUN,sel ,iPEAKSFUN,ProtSequence=PSeq$SEQUENCE)


preparingNetwork <- function(ai,grps,peaksfun,mwinfo){
  ai <<- ai
  grps <<- grps
  peaksfun <<- peaksfun
  mwinfo <<- mwinfo
  
  colo <- grps$Group
  SequenceInfo <- data.table(Sequence=ai$Sequences,Cluster=colo)
  md <- loess(mwinfo$MW~mwinfo$fraction)
  NetWorkFun <- SequenceInfo[,
                             {
                               se <- .BY$Se
                               MappedSeq <- unlist(ai$SequenceInfo$ProteinSequenceMapping[se==ai$SequenceInfo$PeptideSequence])
                               MappedSeq <- unique(MappedSeq)
                               PR <- data.table(ai$SequenceInfo$ProteinSequence[MappedSeq,])
                               PR$ProteinSequence <- PR$SEQUENCE
                               PR$SEQUENCE <- NULL
                               PEAKS <- grps[Sequence==se]
                               if(length(PEAKS)>0){
                                 PEAKS <- as.numeric(unlist(strsplit(PEAKS$Group,";")))
                                 
                                 mw <- predict(md,PEAKS)
                               }else{
                                 mw <- NA
                               }
                              
                               PR <- PR[,2^as.numeric(mw),.(ENSEMBL_PROTEIN,ProteinSequence,MW)]
                               PR
                             },.(Sequence,Cluster)]
  
  NetWorkFun$Diff <- (NetWorkFun$V1-NetWorkFun$MW)
  ab <- (NetWorkFun$Diff)
  # ab <- ab/max(ab)
  diffvec <- rep(0,length(ab))
  diffvec[ab>10000] <- 0.1
  diffvec[ab< -10000] <- -0.1
  
  NetWorkFun$V1 <- round(NetWorkFun$V1)
  NetWorkFun$MW <- round(NetWorkFun$MW)
  nw_edges1 <- data.frame(to=NetWorkFun$Sequence,from = NetWorkFun$V1)
  nw_edges2 <- data.frame(to=NetWorkFun$Sequence,from = NetWorkFun$MW)
  # nw_edges2 <- data.frame(from=NetWorkFun$Sequence,to = NetWorkFun$V1,dashes=NetWorkFun$Diff>0,label=diffvec)
  nw_edges <- rbind(nw_edges1,nw_edges2)
  nw_edges <- unique(nw_edges)
  nw_nodes1 <- NetWorkFun[,.(id=.BY$Sequence,label=Sequence,group=unique(Cluster),value=0),Sequence]
  
    nw_nodes2 <- NetWorkFun[,.(label=V1,id=as.character(unique(V1)),shape="square",group="Observed_MW",value=V1),V1]
  nw_nodes3 <- NetWorkFun[,.(label=paste("U:",MW,sep=""),id=as.character(unique(MW)),shape="square",group="MQ_MW",value=MW),MW]
  
  nw_nodes <- rbind(nw_nodes1,nw_nodes2,nw_nodes3,fill=T)
  li <- list(edges=nw_edges,nodes=nw_nodes,AssembledTable=NetWorkFun)
  class(li) <- c("visnet","list")
  
  return(li)
}
plot.visnet <- function(x){
  vn <- visNetwork(edges = x$edges,nodes=x$nodes)%>% visOptions(manipulation = TRUE)
  vn%>% visOptions(highlightNearest = TRUE)%>% visOptions(manipulation = TRUE)%>% 
    visEdges(arrows = 'from', scaling = list(min = 2, max = 2))
}
FractionToMW <- function(fr,mw,unlog=F){
  if(unlog){
    mw$MW_unlog <- 10^mw$MW
    
  }
  
  mw$MW[mw$fraction==fr]
}

Add_Background <- function(DA,bg){
  background <- data.table(bg)
  colnames(background) <- colnames(DA)
  DA <- rbind(background,DA)
  DA <- DA
  # DA <- DA[rowSums(DA,na.rm = T)!=0,]
}

background_estimator <- function(da){
  di <- dim(da)[2]
  diff <- apply(da[,.SD,.SDcols=c(dim(da)[2],1:dim(da)[2])],1,diff)
  diff <- as.vector(diff)
  diff <- diff[!is.na(diff)]
  
  diff <- t(diff)
  
  background <- matrix(jitter(rep(sd(diff,na.rm = T),28*10),amount = sd(diff,na.rm = T)),nrow=10,ncol=dim(da)[2])
  background <- background + quantile(unlist(da),na.rm = T,0.1)
}

DetectPeak_Wrapper <- function(da,bg,minpeak=0,colo=NULL,background_remove = NULL,repeats=20){
  dagg <- NULL
  if(length(background_remove)>0){
    if(length(colo)==dim(da)[1]){
      colo <- colo[-c(background_remove)]
    }
    da <- da[-c(background_remove),]
  }
  if(length(colo)!=dim(da)[1]){
    colo <- rep(1,dim(da)[1])
  }
  # print("data_agg2")
  # da <<- da
  # colo <<- colo
  
  
  dp_pfr <<- DetectPeak_Compiler(da,bg,colo = colo,minpeak = minpeak)
  try(      showNotification("Calculating False Peaks",id="HUMPE",duration = 30,type="message"))
  decQ <- sapply(1:repeats,function(x){
    # cat("\r",x)
    try({
      dp_pfr_dec <- DetectPeak_Compiler(da,bg,colo = colo,minpeak = minpeak,decoy=T)
      scoredefine(dp_pfr_dec)
    })
    
  })
  try( removeNotification("HUMPE"))
  
  a <- scoredefine(dp_pfr)
  
  CutOff <- quantile(as.numeric(unlist(decQ)),probs = 0.90,na.rm = T)
  
  peakfindResult <- dp_pfr[a>CutOff]
  
  peaks_result <- peakfindResult[,.SD[,.(min=quantile(V3,probs=0.2,na.rm = T),max=quantile(V4,probs=0.8,na.rm = T),max(V1,na.rm = T),co=paste(unique(co),collapse=";")),V2],]
  maxpeakwidth <- 3
  peaksfuntemp <- peaks_result
  peaksfuntemp[,min:={sel <- (V2-min)>maxpeakwidth ;min[sel] <- V2[sel]-maxpeakwidth;min},V2]
  peaksfuntemp[,max:={sel <- (max-V2)>maxpeakwidth ;max[sel] <- V2[sel]+maxpeakwidth;max},V2]
  peaksfuntemp
}
calc_SuSq <- function(df) sum(as.matrix(dist(df)^2)) / (2 * nrow(df))

ggplot_MatrixLine <- function(x,col=1:dim(x)[1],expectedLocation=NULL,background = 1:10,ClusterMean=T,x_replacement=NULL){
  x <<- x
  coli <<- col
  background <<- background
  print("GGPLOT")
  
  # t.test
  x$ID<- 1:dim(x)[1]
  x$Color=col
  if(length(background)>0){
    bg <- x[background,]
    sample <- x[11:dim(x)[1]]
    sigi  <- apply(sample,1,function(x){
      sapply(x,function(y){
        p <- NA
        # try(p <- t.test(unlist(bg),mu=y,alternative = "less")$p.value)
        p
      })
    })
    uniCol <- unique(x$Color[background])
    sel <- !is.na(match(x$Color,uniCol))
    if(any(sel)){
      x$Color[sel] <- paste(x$Color[sel],"background match")
    }
  }else{
    sample <- x#[11:dim(x)[1]]
    
    sigi  <- apply(sample,1,function(x){
      sapply(x,function(y){
        p <- NA
        # try(p <- t.test(unlist(bg),mu=y,alternative = "less")$p.value)
        p
      })
    })
    
  }
  
  #


  # x$ID[background] <- "background"
  if(length(background)>0){
    fu <- melt(x[-background],id.vars=c("ID","Color"))
    fu$p <- 1
    # fu[,p:={
    #   id <<- .BY$ID
    #   ps <- sigi[,id-max(background,na.rm = T)]
    #   ps[match(names(ps),variable)]
    # },ID]
  }else{
    fu <- melt(x,id.vars=c("ID","Color"))
    fu$p <- 1
    
    # fu[,p:=NA,ID]
  }
  
  bg1 <- melt(x[(background)],id.vars=c("ID","Color"))
  bg1$p <- 1
  if(length(x_replacement)>0){
    # x_replacement <<- x_replacement
    x_replacement <- x_replacement[!is.na(x_replacement$fraction)]
    fu[,mw:=FractionToMW(gsub("V","",.BY$variable),x_replacement),variable]
    bg1[,mw:=FractionToMW(gsub("V","",.BY$variable),x_replacement),variable]
    
  }else{
    fu$mw <- fu$variable
    bg1$mw <- bg1$variable
  }
  
  # f
  fu <<- fu
  fu <- fu[!is.na(value)]
  g1 <- ggplot(fu,aes(mw,value,group=(ID),shape=p<0.01,color=(Color)))+ylab("Intensity")+xlab("Fractions (increasing Protein MW)")
  g1 <- g1+ theme(legend.title = element_blank())
  if(ClusterMean){
    g1 <- g1+geom_point()
    g1 <- g1+stat_summary( fun.y=mean,aes(group=Color), geom="line")
  }else{
    # g1 <- g1+geom_point(aes(color=fu$p<0.01))+geom_line()
    g1 <- g1+geom_point()+geom_line()
    
  }
  if(1){
    # g1 <- g1+geom_point(data=bg1,aes(mw,value,group=rep(1,dim(bg1)[1])),color="darkgrey")
    # g1+geom_point(data=bg1,aes(variable,value,group=rep(1,dim(bg1)[1])),color="black")
  }
  g1 <- g1+ guides(fill=guide_legend(title="New Legend Title"))
  print("GGPLOTDONE")
  g1
}
# g1g1 <- ggplot_MatrixLine(DaMo[-c(1:10)],col = factor(colo,rev(unique(colo))),background = NULL,ClusterMean = F)

smoothFun <- function(x,xVec=NULL,type="supsmu",k=3,...){
  co <- colnames(x)
  
  x <- x[,apply(.SD,1,function(y){
    if(length(xVec)==0){
      x <- 2:29
    }else{
      x <- xVec
    }
    y <- y
    
    if(type=="supsmu"){
      ks <- "cv"
      if(k!=1){
        ks=k/50
      }else{
        ks=="cv"
      }
      y <- y[is.na(y)] <- 0
      y <- supsmu(x,y,span=ks,...)$y
    }
    if(type=="runmed"){
      y <- runmed(y,k,...)
    }
    if(type=="smooth.spline"){
      y <- smooth.spline(x,y,spar=k/20)$y
    }
    if(type=="smooth"){
      y <- smooth(y,...)
    }
    
    y
  }),]
  x <- as.data.table(t(x))
  colnames(x) <- co
  x <- x
}
DetectPeak_Compiler <- function(da,bg,colo=rep(1,dim(da)[1]),minpeak=1,decoy=F){
  if(decoy){
    temp <- matrix(sample(unlist(da)),ncol = dim(da)[2],nrow = dim(da)[1],)
    colnames(temp) <- colnames(da)
    da <- data.table(temp)
    dafu <<- da
    colo <- sample(colo)
  }else{
    datru <<- da
  }
  soq <- da[,{
    temp <<- .SD
    soq <- apply(temp,2,function(x){
      hum <- x-median(x,na.rm = T);
      l <- sum(!is.na(x));sum(hum^2,na.rm = T)/l
    })
    as.list(soq)
  },.(co=colo)]
  
  dagg <- da[,{
    temp <<- .SD
    co <- apply(temp,2,median,na.rm = T)
    cn <- colnames(temp)
    if(is.vector(co)){
      co <- matrix(co,ncol=length(co))
    }
    colnames(co) <- cn
    co[is.na(co)] <- sample(bg,sum(is.na(co)))
    co <- data.table(co)
    co$n =dim(temp)[1]
    co
  },.(co=colo)]
  
  
  numbi <- 4
  peakfindResult <- dagg[,{
    temp <- unlist(.SD)
    spreadi <- quantile(temp,probs = 0.1,na.rm = T)
    spreadi <- jitter(rep(spreadi,numbi))
    temp <- c(spreadi,temp,spreadi)
    data.table(findpeaks(temp,
                         minpeakheight=minpeak,
                         nups	=1,ndowns=0,
                         minpeakdistance = 2))},.(co,n)]
  
  peakfindResult$V3 <- peakfindResult$V3-numbi
  peakfindResult$V2 <- peakfindResult$V2-numbi
  peakfindResult$V4 <- peakfindResult$V4-numbi
  peakfindResult$V4[peakfindResult$V4<=dim(da)[2]] <- dim(da)[2]
  peakfindResult <- peakfindResult[V2<=dim(da)[2],]
  # peakfindResult <<- peakfindResult
  
  # peakfindResult[,soq:={
  #   temp <<- .SD
  #   byvec <<- .BY$co
  #   s <- NA
  #   try(s <- soq[soq$co==byvec,.SD,.SDcols=as.character(byvec)])
  #   s
  # },.(V2,co)]
  peakfindResult
}

DetectPeak_Compiler_old <- function(da,colo=rep(1,dim(da)[1]),minpeak=1,decoy=F){
  if(decoy){
    temp <- matrix(sample(unlist(da)),ncol = dim(da)[2],nrow = dim(da)[1],)
    colnames(temp) <- colnames(da)
    da <- data.table(temp)
    dafu <<- da
    
  }
  soq <- da[,{
    temp <<- .SD
    soq <- apply(temp,2,function(x){
      hum <- x-median(x,na.rm = T);
      l <- sum(!is.na(x));sum(hum^2,na.rm = T)/l
    })
    as.list(soq)
  },.(co=colo)]
  
  dagg <- da[,{
    temp <<- .SD
    co <- apply(temp,2,median,na.rm = T)
    cn <- colnames(temp)
    if(is.vector(co)){
      co <- matrix(co,ncol=length(co))
    }
    colnames(co) <- cn
    co[is.na(co)] <- sample(bg,sum(is.na(co)))
    co <- data.table(co)
    co$n =dim(temp)[1]
    co
  },.(co=colo)]
  
  
  peakfindResult <- dagg[,{data.table(findpeaks(unlist(.SD),
                                                minpeakheight=minpeak,
                                                nups	=1,ndowns=1,
                                                minpeakdistance = 3))},.(co,n)]
  # peakfindResult[,humpe:={
  #   temp <<- .SD
  #   s <<- soq[soq$co==.BY$co,.SD,.SDcols=.BY$V2]
  #   s
  # },.(V2,co)]
  peakfindResult
}
CompareVecs <- function(a,b){
  v <- c(a,b)
  c <- c(rep(2,length(a)),rep(1,length(b)))
  o <- order(v)
  return(list(v=v[o],c=c[o]))
}
scoredefine <- function(x){
  x$V1*x$n 
}


ISOFRAC_plot <- function(input){
  rd_c <- input$rdc
  Peaks <- input$Peaks
  
  NewAssignments <- input$Assignmeints
  diff_candis <- input$Diff_Candis
  Quality <- input$Qua
  AveragePeaks <- input$AveragePeaks
  x <- input$GeneName
  
  g2 <- ggplot(rd_c,aes(num,value,group=Clust,col=Clust))+geom_line()+ggtitle(x)
  try({
    for(i in 1:dim(Peaks)[1]){
      if(any(i==diff_candis)){
        RealPeak <- NewAssignments[diff_candis==i]
        
      }else{
        RealPeak <- i
      }
      q <- Quality$predict[i]
      if(q<0.2){
        q <- "grey"
      }else{
        q <- "red"
      }
      temp <- Peaks[i,]
      
      g2 <- g2 + annotate("pointrange", x = temp$V2, y = temp$V1, xmin = temp$V3, xmax = temp$V4,
                          color = q, size = Quality$predict[i]*2)
      g2 <- g2 + annotate("text", x = temp$V2, y = temp$V1,label=RealPeak)
      ap <- AveragePeaks[i, ]
      
      g2 <- g2+annotate("pointrange", x = ap$M, y = ap$H, xmin = ap$M1, xmax =ap$M2,
                        color = "#22880080", size = Quality$predict[i]*2)
      
      g2 <- g2+annotate("text", x = ap$M, y = ap$H, label=RealPeak ,
                        color = "black")
      
    }
    
  })
  g2
}

Neighbour_Mean_imputation <- function(rawda){
  rawdanew <- apply(rawda,1,function(xu){
    xu <- xu
    NAs <- which(!is.na(xu))
    if(any(is.na(xu))){
      repi <- sapply(1:length(xu),function(y){
        v1 <<- xu[y]
        y <<- y
        if(is.na(v1)){
          mi <- max(NAs[NAs<y])
          ma <- min(NAs[NAs>y])
          v1 <- mean(xu[c(mi,ma)],na.rm = T)
        }
        v1
        
      })
      
    }else{
      repi <- xu
    }
  })
  rawdanew <- t(rawdanew)
  rawda <- data.table(rawdanew)
  rawda
}
PreparingData <- function(it,da){
  # library(h2o)
  # library(data.table)
  # library(pracma)
  # library(ggplot2)
  # library(parallel)
  # library(ggplot2)
  # library(cluster)
  # 
  source("Isofrac_Functions.R")
  rawda <- NULL
  # init <- h2o.init()
  # if(!exists("it")){
  #   it <<- 1
  #   le <<- 0
  # }
  nami <- Sys.info()[['nodename']]
  write(it,paste("progress",nami,Sys.getpid(),".txt",sep = ""))
  # x <<- x
  li <- list()
  try({
    x <<- unique(da$dt_gn$GN_ENS)[it]
    GeneName <- x
    sel <- da$dt_gn$GN_ENS==x
    sum(sel)
    rawda <- da$singlenorm_zscore[sel,]
    if(is.vector(rawda)){
      rawda <- t(data.table(da$singlenorm_zscore[sel,]))
      
    }
    if(all(is.na(rawda))){
      print("all NA")
      return(li)
    }
  })
  NAswitch    <- rowSums(as.data.frame(rawda),na.rm = T)!=0
  SequenceID  <- which(NAswitch)
  SequenceID  <- which(sel)[SequenceID]
  
  rawda <- rawda[NAswitch,]
  if(is.vector(rawda)){
    rawda <- t(data.table(rawda))
    
  }
  return(list(rawda=rawda,SequenceID=SequenceID))
}
Feature_Extraction_Wrapper <- function(rawda,SequenceID=NULL,mw=NULL,ReportFeatureTable=F,FalseOutput=F,IncludeFalsePeaks=F,kvec=4){
  # rawda <<- rawda
  
  if(length(SequenceID)==0){
    SequenceID <- 1:dim(rawda)[1]
  }
  # impute NA as median between neighboring points
  rawda <- Neighbour_Mean_imputation(rawda)
  
  
  if(is.vector(rawda)){
    rawda <- t(data.table(rawda))
    rawda_sm <- smoothFun(rawda,xVec =sapply(1:28,FractionToMW,mw),type="runmed",k=3)
    fa <- list(clustering=0)
    Mfa <- list(cluster=0)
    
  }else{
    rawda_sm <- smoothFun(rawda,xVec =sapply(1:28,FractionToMW,mw),type="runmed",k=3)
    rawda_sm[is.na(rawda_sm)] <- jitter(quantile(rawda_sm,na.rm = T,probs = 0.01))
    if(dim(rawda_sm)[1]<kvec){
      kvec <- 2
    }
    
    fa <- list(clustering=1:dim(rawda_sm)[1])
    Mfa <- list(cluster=1:dim(rawda_sm)[1])
    try({Mfa <- kmeans((dist(as.matrix(rawda_sm))),centers = kvec)})
    fa$clustering <- Mfa$cluster
    checkTest  <- try(fa <- fanny((rawda_sm),k = kvec,metric="SqEuclidean",stand=F,memb.exp	=2),silent = F)
    if(class(checkTest)=="try-error"){
      checkTest  <- try(fa <- fanny((rawda_sm),k = 2),silent = F)
    }
    hdbscan_result <- list("cluster" =0, "minPts" =0, "cluster_scores"=0 , "membership_prob"=0, "outlier_scores" =0, "hc"=0)
    try({
      hdbscan_result <- hdbscan(rawda_sm,minPts = 2)
    })
    
    
  }
  SequenceID_Clusters <- data.table(SequenceID,fz=fa$clustering,km=Mfa$cluster,hdbs_clusters=hdbscan_result$cluster,hdbs_membershipscore=hdbscan_result$membership_prob,hdbs_outlierscore=hdbscan_result$outlier_scores,hdbs_clusterscore=hdbscan_result$cluster_scores[hdbscan_result$cluster])
  
  # check SOQ of the different Clusters:
  # PeptideSOQ_fz <- rawda_sm[,.(SOS=calc_SuSq(.SD),n=dim(.SD)[1],(which(fa$clustering==.BY$fa))),fa$clustering]
  # PeptideSOQ_km <- rawda_sm[,.(SOS=calc_SuSq(.SD),n=dim(.SD)[1],(which(Mfa$cluster==.BY$Mfa))),Mfa$cluster]
  # SOQ_Merged <- merge(PeptideSOQ_fz,PeptideSOQ_km,by="V3")
  
  
  # remove noise again:
  # Noisecluster <- fa$clustering[1:20]
  # smoothing
  
  
  rawda_sm$Col <- as.numeric(hdbscan_result$cluster)+1#as.character(fa$clustering)
  rawda_sm$Col_km <- as.character(Mfa$cluster)
  rawda_sm$Col_hdbs <- as.character(hdbscan_result$cluster+1)
  
  # Setting clusters to minus if pairwise correlations is below threshold
  rawda_sm[,Col:= {
    temp <- .SD
    temp$Col_km <- NULL
    temp$Col_hdbs <- NULL
    
    HUM <- apply(cor(t(temp)),2,median,na.rm = T)
    Clust <- rep(.BY$Col,dim(temp)[1])
    if(any(HUM<0.5)){
      Clust[HUM<0.5] <- -as.numeric(.BY$Col)
    }
    Clust
  },Col]
  
  # Preparing Running median for clusters
  Cluster <- rawda_sm[,{
    temp <- .SD
    temp$Col_km <- NULL
    temp$Seq <- NULL
    temp$Col_hdbs <- NULL
    
    M <- apply(temp,2,median,na.rm = T)
    M <- runmed(M,3)
    
    as.list(M)
  },.(Clust=Col)]
  # Preparing Running median for clusters fuzzy
  
  Cluster_km <- rawda_sm[,{
    temp <- .SD
    temp$Col <- NULL
    temp$Seq <- NULL
    temp$Col_hdbs <- NULL
    
    M <- apply(temp,2,median,na.rm = T)
    M <- runmed(M,3)
    
    as.list(M)
  },.(Clust=Col_km)]
  # Preparing Running median for clusters hdbs
  
  # Cluster_hdbs <- rawda_sm[,{
  #   temp <- .SD
  #   temp$Col <- NULL
  #   temp$Seq <- NULL
  #   temp$Col_hdbs <- NULL
  #   
  #   M <- apply(temp,2,median,na.rm = T)
  #   M <- runmed(M,3,na.rm = T)
  #   
  #   as.list(M)
  # },.(Clust=Col_hdbs)]
  
  # rawda_sm(,Clust)
  # calc_SuSq()
  
  # Find best Peak
  SIPEAKS <- rawda_sm[,data.table(findpeaks_Set(unlist(.SD))),.SDcols=grep("^V",names(rawda_sm)),.(Col,ID=1:dim(rawda_sm)[1])]
  SIPEAKS_km <- rawda_sm[,data.table(findpeaks_Set(unlist(.SD))),.SDcols=grep("^V",names(rawda_sm)),.(Col=Col_km,ID=1:dim(rawda_sm)[1])]
  # SIPEAKS$Type <- "fz"
  SIPEAKS$Type <- "hdbs"
  
  SIPEAKS_km$Type <- "km"
  SIPEAKS <- rbind(SIPEAKS_km,SIPEAKS)
  data.table::setcolorder(SIPEAKS,c("V1","V2","V3","V4","Threshold","Type","ID"))
  
  rawda_sm_pure <- copy(rawda_sm)
  rawda_sm_pure$Seq <- NULL
  rawda_sm_pure$Col <- NULL
  rawda_sm_pure$Col_km <- NULL
  rawda_sm_pure$Col_hdbs <- NULL
  # 
  # FeatureTable_SIPEAKS <- apply(data.frame(SIPEAKS)[,1:4],1,Feature_Extraction,tempx=rawda_sm_pure)
  # FeatureTable <- rbindlist(FeatureTable)
  # FeatureTable_h2o <- as.h2o(FeatureTable)
  # Quality <- as.data.frame(h2o.predict(PeakModel2,FeatureTable_h2o))
  # SIPEAKS$Quality <- Quality$predict
  
  print("HUM1")
  
  SOS <- apply(SIPEAKS,1,function(temp){
    temp <<- as.numeric(temp)
    if(temp[7]=="km"){
      Peak_SOS <- calc_SuSq(df <- as.data.frame(rawda_sm_pure)[rawda_sm_pure$Col_km==temp[8],temp[3]:temp[4]])
    }else{
      Peak_SOS <- calc_SuSq(df <- as.data.frame(rawda_sm_pure)[rawda_sm_pure$Col==temp[8],temp[3]:temp[4]])
      
    }
    c(Peak_SOS,dim(df)[1])
  })
  SOS <- t(SOS)
  colnames(SOS) <- c("SOS","Cluster_n")
  print("HUM1.1")
  
  SIPEAKS <- cbind(SIPEAKS,SOS)
  # Position based on all Sequences in Cluster
  print("HUM1.2")
  
  rd_c <- melt(Cluster,id.vars =   "Clust")
  rd_c$num <- as.numeric(gsub("V","",rd_c$variable))
  print("HUM1.3")
  Cluster <<- Cluster
  # running findpeaks on Cluster set
  Peaks <- Cluster[,{
    temp <<- .SD
    
    if(dim(.SD)[1]>1){
      warning("More than one entry")
    }
    
    fi <- findpeaks_Set(unlist(temp[1,]))
    if(length(fi)>0){
      SOS_temp <- apply(fi,1,function(x){
        # xx
        df <- as.data.frame(rawda_sm)[rawda_sm$Col==.BY$Clust,x[3]:x[4]]
        Peak_SOS <- calc_SuSq(df)
        c(Peak_SOS,dim(df)[1])
      })
      SOS_temp <- t(SOS_temp)
      fi$SOS <- as.double(SOS_temp[,1])
      fi$Peptides <- as.double(SOS_temp[,2])
      fi <- data.table(fi)
      
    }else{    
      fi <-  NULL
    }
    
    # print(dim(fi))
  },Clust]  
  
  print("HUM2")
  Cluster_km <- Cluster_km
  # running findpeaks on Cluster km 
  
  Peaks_km <- Cluster_km[,{
    temp <<- .SD
    
    if(dim(.SD)[1]>1){q
      warning("More than one entry")
    }
    
    fi <- findpeaks_Set(unlist(temp[1,]))
    if(length(fi)>0){
      
      SOS_temp <- apply(fi,1,function(x){
        # xx
        df <- as.data.frame(rawda_sm)[rawda_sm$Col_km==.BY$Clust,x[3]:x[4]]
        Peak_SOS <- calc_SuSq(df)
        c(Peak_SOS,dim(df)[1])
      })
      SOS_temp <- t(SOS_temp)
      fi$SOS <- SOS_temp[,1]
      fi$Peptides <- SOS_temp[,2]
      fi <- data.table(fi)
    }else{fi <- NULL}
  },Clust]
  Peaks$Type <- "hdbs"
  Peaks_km$Type <- "km"
  Peaks <<- Peaks
  Peaks_km <<- Peaks_km
  if(dim(Peaks)[1]==0){
    Peaks <- Peaks_km# Cluster Based Analysis
    
  }else{
    if(dim(Peaks_km)[1]==0){
      
    }else{
      Peaks <- rbind(Peaks,Peaks_km)# Cluster Based Analysis
      
    }
  }
  Peaks[Threshold>0.1,DifferentPeaks:={
    AveragePeaks <<- .SD
    if(dim(AveragePeaks)[1]>1){
      pw <- combn(PeaksIndex <- 1:dim(AveragePeaks)[1],2)
      diffis <- apply(pw,2,function(x){
        try({
          a <<- AveragePeaks[x[1],]
          b <<- AveragePeaks[x[2],]
          p <- NA
          a <- c(seq(a$V3,a$V4,by=1),a$M)#as.numeric(unlist(a)[c(2,5,6)])
          # b <- c(as.numeric(unlist(b)[c(2,5,6)]),)
          b <- c(seq(b$V3,b$V4,by=1),b$M)#as.numeric(unlist(a)[c(2,5,6)]) 
        })
        
        p <- 1
        try(p <- t.test(a,b)$p.value)
        p
      })
      Difference_check <- t(rbind(pw,diffis))
      Difference_check <- data.table(Difference_check)
      Difference_check <- Difference_check[diffis>0.01,]
      
      diff_candis <- unique(c(Difference_check$V1,Difference_check$V2))
      NewAssignments <- sapply(PeaksIndex,function(x){
        temp <- Difference_check[V1==x|V2==x ]
        if(dim(temp)[1]==0){
          return(x)
        }
        min(unlist(c(temp$V1,temp$V2)))
      })
    }else{
      NewAssignments <- AveragePeaks$Clust
    }
    
    NewAssignments
  },Type]
  print("HUM3")
  data.table::setcolorder(Peaks,c("V1","V2","V3","V4","Clust","Type","Threshold","SOS"))
  if(IncludeFalsePeaks){
    try({
      # False Peaks Training: 
      Positions <- 1:28
      
      UNI <- apply(cbind(Peaks$V3,Peaks$V4),1,function(x){x[1]:x[2]})
      NoPeakPotential <- setdiff(Positions,unique(unlist(UNI)))
      WithinPeak <- setdiff(unlist(UNI),Peaks$V2)
      
      input_ids <- NoPeakPotential
      
      
      NoPeaks <- Cluster[,{
        temp <<- .SD
        fi <- lapply(c(NoPeakPotential,WithinPeak),function(x){
          data.table(unlist(temp[1,.SD,.SDcols=x]),x,x-2,x+2)
        })
        
        fi <- rbindlist(fi)
        colnames(fi)<- c("V1","V2","V3","V4")
        
        fi$Type <- c(rep("NoPeak_fz",length(NoPeakPotential)),rep("WiPeak_fz",length(WithinPeak)))
        fi$V3[fi$V3<1] <- 1
        fi$V4[fi$V4>28] <- 28
        
        
        if(length(fi)>0){
          
          SOS_temp <- apply(fi,1,function(x){
            x<<-as.numeric(x)
            df <- as.data.frame(rawda_sm)[rawda_sm$Col==.BY$Clust,x[3]:x[4]]
            Peak_SOS <- calc_SuSq(df)
            c(Peak_SOS,dim(df)[1])
          })
          SOS_temp <- t(SOS_temp)
          fi$SOS <- SOS_temp[,1]
          fi$Peptides <- SOS_temp[,2]
          fi <- data.table(fi)
        }else{fi <- NULL}
        
      },Clust]
      NoPeaks_km <- Cluster_km[,{
        temp <<- .SD
        fi <- lapply(c(NoPeakPotential,WithinPeak),function(x){
          data.table(unlist(temp[1,.SD,.SDcols=x]),x,x-2,x+2)
        })
        
        fi <- rbindlist(fi)
        colnames(fi)<- c("V1","V2","V3","V4")
        
        fi$Type <- c(rep("NoPeak_km",length(NoPeakPotential)),rep("WiPeak_km",length(WithinPeak)))
        fi$V3[fi$V3<1] <- 1
        fi$V4[fi$V4>28] <- 28
        
        
        if(length(fi)>0){
          
          SOS_temp <- apply(fi,1,function(x){
            x<<-as.numeric(x)
            df <- as.data.frame(rawda_sm)[rawda_sm$Col_km==.BY$Clust,x[3]:x[4]]
            Peak_SOS <- calc_SuSq(df)
            c(Peak_SOS,dim(df)[1])
          })
          SOS_temp <- t(SOS_temp)
          fi$SOS <- SOS_temp[,1]
          fi$Peptides <- SOS_temp[,2]
          fi <- data.table(fi)
        }else{fi <- NULL}
        
        
      },Clust]
      
      NoPeaks <<- NoPeaks
      NoPeaks_km <<- NoPeaks_km
      if(dim(NoPeaks)[1]==0){
        NoPeaks <- NoPeaks_km# Cluster Based Analysis
        
      }else{
        if(dim(NoPeaks_km)[1]==0){
          
        }else{
          NoPeaks <- rbind(NoPeaks,NoPeaks_km)# Cluster Based Analysis
          
        }
      }
      NoPeaks$Threshold <- 0
      data.table::setcolorder(NoPeaks,c("V1","V2","V3","V4","Clust","Type","Threshold","SOS"))
      NoPeaks$DifferentPeaks <- NoPeaks$Clust
    })
    Peaks <-rbind(Peaks,NoPeaks)
  }
  
  
  
  
  FeatureTable <- apply(data.frame(Peaks)[,1:4],1,Feature_Extraction,tempx=rawda_sm_pure)
  
  FeatureTable <- rbindlist(FeatureTable)
  FeatureTable$SOS <- Peaks$SOS
  FeatureTable$Peptides <- Peaks$Peptides
  
  FeatureTable$SOS[FeatureTable$SOS==0] <-300
  list(Features=FeatureTable,SequenceID_Clusters=SequenceID_Clusters,Peaks=Peaks,rawda_sm=rawda_sm,hdbscan=hdbscan_result)
}
if(0){
  it <- grep("TMPO",unique(da$dt_gn$GN_ENS))
  
  GN_ENS_analysis <- lapply(it,SingleEntry_Isofrac_PeakSearch,da=da,findpeaks_Set=findpeaks_Set,Feature_Extraction=Feature_Extraction,PeakModel2=PeakModel2,mw=mw)
  scl <- GN_ENS_analysis[[1]]$SequenceID_Clusters
  scl$hdbs_clusters
  ggplot(scl,aes(hdbs_membershipscore,hdbs_outlierscore,color=as.character(hdbs_clusters),size=hdbs_clusterscore))+geom_point()+ggtitle("TMPO HDBSCAN Result")
  plot(GN_ENS_analysis[[1]]$hdbscan)
  
  da1 <- GN_ENS_analysis[[1]]$Data
  da1$Col<- NULL
  da1$Col_km <- NULL
  da1$Col_hdbs <- NULL
  ggplot_MatrixLine(da1,col = scl$hdbs_clusters)
  
}

# Feature_Extraction_Wrapper <- function(rawda,SequenceID=NULL,mw=NULL,ReportFeatureTable=F,FalseOutput=F,IncludeFalsePeaks=F,kvec=4){
#   # rawda <<- rawda
#   
#   if(length(SequenceID)==0){
#     SequenceID <- 1:dim(rawda)[1]
#   }
#   # impute NA as median between neighboring points
#   rawda <- Neighbour_Mean_imputation(rawda)
#   
#   
#   if(is.vector(rawda)){
#     rawda <- t(data.table(rawda))
#     rawda_sm <- smoothFun(rawda,xVec =sapply(1:28,FractionToMW,mw),type="runmed",k=3)
#     fa <- list(clustering=0)
#     Mfa <- list(cluster=0)
#     
#   }else{
#     rawda_sm <- smoothFun(rawda,xVec =sapply(1:28,FractionToMW,mw),type="runmed",k=3)
#     rawda_sm[is.na(rawda_sm)] <- jitter(quantile(rawda_sm,na.rm = T,probs = 0.01))
#     if(dim(rawda_sm)[1]<kvec){
#       kvec <- 2
#     }
#     
#     fa <- list(clustering=1:dim(rawda_sm)[1])
#     Mfa <- list(cluster=1:dim(rawda_sm)[1])
#     try({Mfa <- kmeans((dist(as.matrix(rawda_sm))),centers = kvec)})
#     fa$clustering <- Mfa$cluster
#     checkTest  <- try(fa <- fanny((rawda_sm),k = kvec,metric="SqEuclidean",stand=F,memb.exp	=2),silent = F)
#     if(class(checkTest)=="try-error"){
#       checkTest  <- try(fa <- fanny((rawda_sm),k = 2),silent = F)
#     }
#     
#     
#     
#   }
#   SequenceID_Clusters <- data.table(SequenceID,fz=fa$clustering,km=Mfa$cluster)
#   
#   # check SOQ of the different Clusters:
#   PeptideSOQ_fz <- rawda_sm[,.(SOS=calc_SuSq(.SD),n=dim(.SD)[1],(which(fa$clustering==.BY$fa))),fa$clustering]
#   PeptideSOQ_km <- rawda_sm[,.(SOS=calc_SuSq(.SD),n=dim(.SD)[1],(which(Mfa$cluster==.BY$Mfa))),Mfa$cluster]
#   SOQ_Merged <- merge(PeptideSOQ_fz,PeptideSOQ_km,by="V3")
#   
#   
#   # remove noise again:
#   # Noisecluster <- fa$clustering[1:20]
#   # smoothing
#   
#   
#   rawda_sm$Col <- as.character(fa$clustering)
#   rawda_sm$Col_km <- as.character(Mfa$cluster)
#   
#   rawda_sm[,Col:= {
#     temp <- .SD
#     temp$Col_km <- NULL
#     HUM <- apply(cor(t(temp)),2,median,na.rm = T)
#     Clust <- rep(.BY$Col,dim(temp)[1])
#     if(any(HUM<0.5)){
#       Clust[HUM<0.5] <- -as.numeric(.BY$Col)
#     }
#     Clust
#   },Col]
# 
#   Cluster <- rawda_sm[,{
#     temp <- .SD
#     temp$Col_km <- NULL
#     temp$Seq <- NULL
#     M <- apply(temp,2,median,na.rm = T)
#     M <- runmed(M,3)
#     
#     as.list(M)
#   },.(Clust=Col)]
#   Cluster_km <- rawda_sm[,{
#     temp <- .SD
#     temp$Col <- NULL
#     temp$Seq <- NULL
#     M <- apply(temp,2,median,na.rm = T)
#     M <- runmed(M,3)
#     
#     as.list(M)
#   },.(Clust=Col_km)]
#   
#   # rawda_sm(,Clust)
#   # calc_SuSq()
#   
#   # Find best Peak
#   SIPEAKS <- rawda_sm[,data.table(findpeaks_Set(unlist(.SD))),.SDcols=grep("^V",names(rawda_sm)),.(Col,ID=1:dim(rawda_sm)[1])]
#   SIPEAKS_km <- rawda_sm[,data.table(findpeaks_Set(unlist(.SD))),.SDcols=grep("^V",names(rawda_sm)),.(Col=Col_km,ID=1:dim(rawda_sm)[1])]
#   SIPEAKS$Type <- "fz"
#   SIPEAKS_km$Type <- "km"
#   SIPEAKS <- rbind(SIPEAKS_km,SIPEAKS)
#   data.table::setcolorder(SIPEAKS,c("V1","V2","V3","V4","Threshold","Type","ID"))
#   
#   rawda_sm_pure <- copy(rawda_sm)
#   rawda_sm_pure$Seq <- NULL
#   rawda_sm_pure$Col <- NULL
#   rawda_sm_pure$Col_km <- NULL
#   # 
#   # FeatureTable_SIPEAKS <- apply(data.frame(SIPEAKS)[,1:4],1,Feature_Extraction,tempx=rawda_sm_pure)
#   # FeatureTable <- rbindlist(FeatureTable)
#   # FeatureTable_h2o <- as.h2o(FeatureTable)
#   # Quality <- as.data.frame(h2o.predict(PeakModel2,FeatureTable_h2o))
#   # SIPEAKS$Quality <- Quality$predict
#   
#   print("HUM1")
#   
#   SOS <- apply(SIPEAKS,1,function(temp){
#     temp <<- as.numeric(temp)
#     if(temp[7]=="km"){
#       Peak_SOS <- calc_SuSq(df <- as.data.frame(rawda_sm_pure)[rawda_sm_pure$Col_km==temp[8],temp[3]:temp[4]])
#     }else{
#       Peak_SOS <- calc_SuSq(df <- as.data.frame(rawda_sm_pure)[rawda_sm_pure$Col==temp[8],temp[3]:temp[4]])
#       
#     }
#     c(Peak_SOS,dim(df)[1])
#   })
#   SOS <- t(SOS)
#   colnames(SOS) <- c("SOS","Cluster_n")
#   print("HUM1.1")
#   
#   SIPEAKS <- cbind(SIPEAKS,SOS)
#   # Position based on all Sequences in Cluster
#   print("HUM1.2")
#   
#   rd_c <- melt(Cluster,id.vars =   "Clust")
#   rd_c$num <- as.numeric(gsub("V","",rd_c$variable))
#   print("HUM1.3")
#   Cluster <<- Cluster
#   Peaks <- Cluster[,{
#     temp <<- .SD
#     
#     if(dim(.SD)[1]>1){
#       warning("More than one entry")
#     }
#     
#     fi <- findpeaks_Set(unlist(temp[1,]))
#     if(length(fi)>0){
#       SOS_temp <- apply(fi,1,function(x){
#         # xx
#         df <- as.data.frame(rawda_sm)[rawda_sm$Col==.BY$Clust,x[3]:x[4]]
#         Peak_SOS <- calc_SuSq(df)
#         c(Peak_SOS,dim(df)[1])
#       })
#       SOS_temp <- t(SOS_temp)
#       fi$SOS <- as.double(SOS_temp[,1])
#       fi$Peptides <- as.double(SOS_temp[,2])
#       fi <- data.table(fi)
#       
#     }else{    
#       fi <-  NULL
#     }
#     
#     # print(dim(fi))
#   },Clust]  
#   
#   print("HUM2")
#   Cluster_km <- Cluster_km
#   Peaks_km <- Cluster_km[,{
#     temp <<- .SD
#     
#     if(dim(.SD)[1]>1){q
#       warning("More than one entry")
#     }
#     
#     fi <- findpeaks_Set(unlist(temp[1,]))
#     if(length(fi)>0){
#       
#     SOS_temp <- apply(fi,1,function(x){
#       # xx
#       df <- as.data.frame(rawda_sm)[rawda_sm$Col_km==.BY$Clust,x[3]:x[4]]
#       Peak_SOS <- calc_SuSq(df)
#       c(Peak_SOS,dim(df)[1])
#     })
#     SOS_temp <- t(SOS_temp)
#     fi$SOS <- SOS_temp[,1]
#     fi$Peptides <- SOS_temp[,2]
#     fi <- data.table(fi)
#     }else{fi <- NULL}
#   },Clust]
#   Peaks$Type <- "fz"
#   Peaks_km$Type <- "km"
#   Peaks <<- Peaks
#   Peaks_km <<- Peaks_km
#   if(dim(Peaks)[1]==0){
#     Peaks <- Peaks_km# Cluster Based Analysis
#     
#   }else{
#     if(dim(Peaks_km)[1]==0){
#       
#     }else{
#       Peaks <- rbind(Peaks,Peaks_km)# Cluster Based Analysis
#       
#     }
#   }
#   Peaks[Threshold>0.1,DifferentPeaks:={
#     AveragePeaks <<- .SD
#     if(dim(AveragePeaks)[1]>1){
#       pw <- combn(PeaksIndex <- 1:dim(AveragePeaks)[1],2)
#       diffis <- apply(pw,2,function(x){
#         try({
#           a <<- AveragePeaks[x[1],]
#           b <<- AveragePeaks[x[2],]
#           p <- NA
#           a <- c(seq(a$V3,a$V4,by=1),a$M)#as.numeric(unlist(a)[c(2,5,6)])
#           # b <- c(as.numeric(unlist(b)[c(2,5,6)]),)
#           b <- c(seq(b$V3,b$V4,by=1),b$M)#as.numeric(unlist(a)[c(2,5,6)]) 
#         })
#         
#         p <- 1
#         try(p <- t.test(a,b)$p.value)
#         p
#       })
#       Difference_check <- t(rbind(pw,diffis))
#       Difference_check <- data.table(Difference_check)
#       Difference_check <- Difference_check[diffis>0.01,]
#       
#       diff_candis <- unique(c(Difference_check$V1,Difference_check$V2))
#       NewAssignments <- sapply(PeaksIndex,function(x){
#         temp <- Difference_check[V1==x|V2==x ]
#         if(dim(temp)[1]==0){
#           return(x)
#         }
#         min(unlist(c(temp$V1,temp$V2)))
#       })
#     }else{
#       NewAssignments <- AveragePeaks$Clust
#     }
#     
#     NewAssignments
#   },Type]
#   print("HUM3")
#   data.table::setcolorder(Peaks,c("V1","V2","V3","V4","Clust","Type","Threshold","SOS"))
#   if(IncludeFalsePeaks){
#     try({
#       # False Peaks Training: 
#       Positions <- 1:28
#       
#       UNI <- apply(cbind(Peaks$V3,Peaks$V4),1,function(x){x[1]:x[2]})
#       NoPeakPotential <- setdiff(Positions,unique(unlist(UNI)))
#       WithinPeak <- setdiff(unlist(UNI),Peaks$V2)
#       
#       input_ids <- NoPeakPotential
#       
#       
#       NoPeaks <- Cluster[,{
#         temp <<- .SD
#         fi <- lapply(c(NoPeakPotential,WithinPeak),function(x){
#           data.table(unlist(temp[1,.SD,.SDcols=x]),x,x-2,x+2)
#         })
#         
#         fi <- rbindlist(fi)
#         colnames(fi)<- c("V1","V2","V3","V4")
#         
#         fi$Type <- c(rep("NoPeak_fz",length(NoPeakPotential)),rep("WiPeak_fz",length(WithinPeak)))
#         fi$V3[fi$V3<1] <- 1
#         fi$V4[fi$V4>28] <- 28
#         
#         
#         if(length(fi)>0){
#           
#           SOS_temp <- apply(fi,1,function(x){
#             x<<-as.numeric(x)
#             df <- as.data.frame(rawda_sm)[rawda_sm$Col==.BY$Clust,x[3]:x[4]]
#             Peak_SOS <- calc_SuSq(df)
#             c(Peak_SOS,dim(df)[1])
#           })
#           SOS_temp <- t(SOS_temp)
#           fi$SOS <- SOS_temp[,1]
#           fi$Peptides <- SOS_temp[,2]
#           fi <- data.table(fi)
#         }else{fi <- NULL}
#         
#       },Clust]
#       NoPeaks_km <- Cluster_km[,{
#         temp <<- .SD
#         fi <- lapply(c(NoPeakPotential,WithinPeak),function(x){
#           data.table(unlist(temp[1,.SD,.SDcols=x]),x,x-2,x+2)
#         })
#         
#         fi <- rbindlist(fi)
#         colnames(fi)<- c("V1","V2","V3","V4")
#         
#         fi$Type <- c(rep("NoPeak_km",length(NoPeakPotential)),rep("WiPeak_km",length(WithinPeak)))
#         fi$V3[fi$V3<1] <- 1
#         fi$V4[fi$V4>28] <- 28
#         
#         
#         if(length(fi)>0){
#           
#           SOS_temp <- apply(fi,1,function(x){
#             x<<-as.numeric(x)
#             df <- as.data.frame(rawda_sm)[rawda_sm$Col_km==.BY$Clust,x[3]:x[4]]
#             Peak_SOS <- calc_SuSq(df)
#             c(Peak_SOS,dim(df)[1])
#           })
#           SOS_temp <- t(SOS_temp)
#           fi$SOS <- SOS_temp[,1]
#           fi$Peptides <- SOS_temp[,2]
#           fi <- data.table(fi)
#         }else{fi <- NULL}
#         
#         
#       },Clust]
#       
#       NoPeaks <<- NoPeaks
#       NoPeaks_km <<- NoPeaks_km
#       if(dim(NoPeaks)[1]==0){
#         NoPeaks <- NoPeaks_km# Cluster Based Analysis
#         
#       }else{
#         if(dim(NoPeaks_km)[1]==0){
#           
#         }else{
#           NoPeaks <- rbind(NoPeaks,NoPeaks_km)# Cluster Based Analysis
#           
#         }
#       }
#       NoPeaks$Threshold <- 0
#       data.table::setcolorder(NoPeaks,c("V1","V2","V3","V4","Clust","Type","Threshold","SOS"))
#       NoPeaks$DifferentPeaks <- NoPeaks$Clust
#     })
#     Peaks <-rbind(Peaks,NoPeaks)
#   }
# 
#   
#   
#   
#   FeatureTable <- apply(data.frame(Peaks)[,1:4],1,Feature_Extraction,tempx=rawda_sm_pure)
#   
#   FeatureTable <- rbindlist(FeatureTable)
#   FeatureTable$SOS <- Peaks$SOS
#   FeatureTable$Peptides <- Peaks$Peptides
#   
#   FeatureTable$SOS[FeatureTable$SOS==0] <-300
#   list(Features=FeatureTable,SequenceID_Clusters=SequenceID_Clusters,Peaks=Peaks,rawda_sm=rawda_sm)
# }
# SingleEntry_Isofrac_PeakSearch <- function(it,da,findpeaks_Set,Feature_Extraction,PeakModel2,mw,ReportFeatureTable=F){
#   # library(h2o)
#   # library(data.table)
#   # library(pracma)
#   # library(ggplot2)
#   # library(parallel)
#   # library(ggplot2)
#   # library(cluster)
#   
#   h2o.init()
#   source("Isofrac_Functions.R")
#   
#   # init <- h2o.init()
#   # if(!exists("it")){
#   #   it <<- 1
#   #   le <<- 0
#   # }
#   nami <- Sys.info()[['nodename']]
#   write(it,paste("progress",nami,Sys.getpid(),".txt",sep = ""))
#   # x <<- x
#   li <- list()
#   try({
#     x <- unique(da$dt_gn$GN_ENS)[it]
#     GeneName <- x
#     sel <- da$dt_gn$GN_ENS==x
#     sum(sel)
#     rawda <- da$singlenorm_zscore[sel,]
#     dim(rawda)
#     if(is.vector(rawda)){
#       rawda <- t(data.table(da$singlenorm_zscore[sel,]))
#       
#     }
#     if(all(is.na(rawda))){
#       print("all NA")
#       return(li)
#     }
#     
#     NAswitch <- rowSums(as.data.frame(rawda),na.rm = T)!=0
#     SequenceID <- which(NAswitch)
#     SequenceID <- which(sel)[SequenceID]
#     
#     rawda <- rawda[NAswitch,]
#     
#     
#     if(is.vector(rawda)){
#       rawda <- t(data.table(rawda))
#       
#     }
#     
#     FeWrap<- Feature_Extraction_Wrapper(rawda,SequenceID = SequenceID)
#     
#     # temp <- FeWrap$rawda_sm
#     # temp$ID <- 1:dim(temp)[1]
#     # te <- melt(rawda,id.vars = c("Col","Col_km","ID"))
#     # ggplot(te,aes(variable,value,color=Col,group=ID))+geom_line()
#     # ggplot(te,aes(variable,value,color=Col_km,group=ID))+geom_line()
#     # # 
#     Predictions <- h2o.predict(PeakModel2,as.h2o(FeWrap$Features))
#     Predictions <- as.data.table(Predictions)
# 
#     # impute NA as median between neighboring points
#     
#     li <-  list(Data=FeWrap$rawda_sm,
#                 Predictions=Predictions,
#                 Peaks=FeWrap$Peaks,
#                 SequenceID_Clusters=FeWrap$SequenceID_Clusters,
#                 Features=FeWrap$Features,
#                 hdbscan = FeWrap$hdbscan,
#                 GeneName=GeneName
#     )
#   },silent=F)
#   if(length(li)==0){
#     # stop(paste(x,"zero li"))
#   }
#   cat("\r",length(li))
#   li
#   
#   
#   
# }

# findpeaks_Set <- function(x,threshold=0.1){
#   # x <<- x
#   
#   FiPeaks <- lapply(seq(0,1,by = 0.1),function(th){
#     fi <- data.table(findpeaks(x,minpeakdistance = 1,nups =1,ndowns = 1,minpeakheight = 0,threshold = th,zero="+"))
#     fi$Threshold=th
#     fi
#   })
#   
#   FiPeaks <- rbindlist(FiPeaks[lengths(FiPeaks)>2])
#   if(dim(FiPeaks)[1]>0){
#     FiPeaks <- FiPeaks[,.SD[Threshold==max(Threshold)],V2]
#   }else{
#     NULL
#   }
#   
# } 
IsoFracPlot <- function(input,type="fz"){
  # input <<- input
  input$Data$ID <- 1:dim(input$Data)[1]
  Long <-melt(input$Data,id.vars = c("Col","Col_km","ID"))
  Long$variable <- as.numeric(gsub("V","",Long$variable))
  x <- type
  print(x)
  if(x=="km"){
    Long$Clust <- Long$Col_km
  }
  if(x=="fz"){
    Long$Clust <- Long$Col
  }
  g1 <- ggplot(Long,aes(variable,value,group=ID,color=Clust))+geom_line()
  
  P <- input$Peaks
  
  P$V2 <- apply(cbind(P$V2,P$V3,P$V4,P$Clust),1,function(fpi){
    fpi <<- fpi
    o <- NA
    try({
      temp <- input$Data[Col==fpi[4],]
      Vali <- melt(temp,id.vars =c( "Col","Col_km"))
      Vali$variable <- as.numeric(gsub("V","",Vali$variable))
      
      x1 <- fpi[2]:fpi[1]
      v1 <- Vali[!is.na(match(Vali$variable,x1))]
      x1 <- v1$variable
      y1 <- v1$value
      lm_bf <- lm(as.numeric(y1)~as.numeric(x1))
      x2 <- fpi[1]:fpi[3]
      v2 <- Vali[!is.na(match(Vali$variable,x2))]
      x2 <- v2$variable
      y2 <- v2$value
      lm_af <- lm(as.numeric(y2)~as.numeric(x2))
      
      
      
      hu<- lmIntx(lm_bf,lm_af)
      o <- hu$x
    })
    o
  })
  
  P$predictions <- input$Predictions$predict
  sel <- input$Peaks$Type==x
  P <- P[sel,]
  # for(i in 1:dim(input$Peaks)[1]){
  #   g1 <- g1+geom_text_repel(P,aes(x=V))
  #     annotate("Text",x=P$V2[i],y=P$V1[i],label=round(P$predictions[i],2))
  # }
  # geom_text_repel()
  g1 <- g1+geom_point(data=P,aes(x=V2,y=V1,group=Clust,color=Clust,label=round(V2,2),size=P$predictions),box.padding=0.1,label.padding=0.1,label.size =0.1,shape = 8)
  
  g1 <- g1+geom_label_repel(data=P,aes(x=V2,y=V1,group=Clust,fill=Clust,label=round(V2,2),size=P$predictions),color="black",box.padding=0.1,label.padding=0.1,label.size =0.1)
  
  g1
  
  
}
# IsoFracPlot(GN_ENS_analysis[[4]])

# Feature_Extraction_Wrapper(rawda$rawda,SequenceID = NULL)
Feature_Extraction <- function(x,n=1 # clustersize
         ,v=1,assigned=0,tempx=t(as.matrix(x) )# relative cluster deviation
){
  
  x <<- as.numeric(x)
  x[1] <- x[1]-min(tempx)+1
  tempx <<- tempx-min(tempx)+1
  
  el <- list()
  el$q0.01 <- quantile(as.numeric(unlist(tempx)),probs=0.01,na.rm = T) # 0.01 Quantile
  
  # PeakFeatures:
  el$height  <- x[1] # Peak Height
  el$height_rel <- x[1]/el$q0.01 # 
  el$height_rel <- x[1]/el$q0.01 # 
  el$height  <- x[1] # Peak Height
  el$height_rel <- x[1]/el$q0.01 # 
  
  el$pw      <- x[4]-x[3] # Peak Width
  # threepoint Slope before
  p_af <- x[2]:(x[2]+2)
  p_bf <- (x[2]-2):x[2]
  #this part could be improved, handling at the hand of the fractions is not clear, but these regions, may shouldn't be taken into consideration anyhow
  if(is.matrix(tempx)|is.data.frame((tempx))){
    tempx_mean <- apply(tempx,2,mean,na.rm = T)
    
  }else{
    tempx_mean <- tempx
  }
  p_af[p_af>length(tempx_mean)] <- NA
  p_bf[p_bf<=0] <- NA
  cof_af <- lm(tempx_mean[p_af]~p_af)$coefficients
  el$cof_af <- cof_af[grep("Intercept",names(cof_af),invert =T)]#slope before
  
  
  cof_bf <- lm(tempx_mean[p_bf]~p_bf)$coefficients
  el$cof_bf <- cof_bf[grep("Intercept",names(cof_bf),invert =T)]#slope after
  el$cof_bf[is.na(el$cof_bf)] <- 0
  el$cof_af[is.na(el$cof_af)] <- 0
  
  el$height_p1  <- tempx_mean[p_af[2]] # Peak Height before
  el$height_m1  <- tempx_mean[p_bf[2]] # Peak Height before
  el$height_p2  <- tempx_mean[p_af[3]] # Peak Height after
  el$height_m2  <- tempx_mean[p_bf[3]] # Peak Height after
  
  el$height_p1  <- tempx_mean[p_af[2]]/el$q0.01 # relative Peak Height before
  el$height_m1  <- tempx_mean[p_bf[2]]/el$q0.01 # relative Peak Height before
  el$height_p2  <- tempx_mean[p_af[3]]/el$q0.01 # relative Peak Height after
  el$height_m2  <- tempx_mean[p_bf[3]]/el$q0.01 # relative Peak Height after
  
  
  
  
  # surrounding peaks
  # Peak Classes:
  el$n <- n
  el$v <- v
  el$assigned <- assigned
  # el <<- el
  # additional info
  el$X <- x[2]
  el$Xb <- x[3]
  el$Xs<- x[4]
  
  if(all(lengths(el)==1)){
    return(as.data.frame(el))
  }else{
    return(NA)
  }
}
Feature_Extraction <- function(x,n=1 # clustersize
                               ,v=1,assigned=0,tempx=t(as.matrix(x) )# relative cluster deviation
){
  
  x <<- as.numeric(x)
  x[1] <- x[1]-min(tempx)+1
  tempx <<- tempx-min(tempx)+1
  
  el <- list()
  el$q0.01 <- quantile(as.numeric(unlist(tempx)),probs=0.01,na.rm = T) # 0.01 Quantile
  
  # PeakFeatures:
  el$height  <- x[1] # Peak Height
  el$height_rel <- x[1]/el$q0.01 # 
  el$height_rel <- x[1]/el$q0.01 # 
  el$height  <- x[1] # Peak Height
  el$height_rel <- x[1]/el$q0.01 # 
  
  el$pw      <- x[4]-x[3] # Peak Width
  # threepoint Slope before
  p_af <- x[2]:(x[2]+2)
  p_bf <- (x[2]-2):x[2]
  #this part could be improved, handling at the hand of the fractions is not clear, but these regions, may shouldn't be taken into consideration anyhow
  if(is.matrix(tempx)|is.data.frame((tempx))){
    tempx_mean <- apply(tempx,2,mean,na.rm = T)
    
  }else{
    tempx_mean <- tempx
  }
  p_af[p_af>length(tempx_mean)] <- NA
  p_bf[p_bf<=0] <- NA
  cof_af <- lm(tempx_mean[p_af]~p_af)$coefficients
  el$cof_af <- cof_af[grep("Intercept",names(cof_af),invert =T)]#slope before
  
  
  cof_bf <- lm(tempx_mean[p_bf]~p_bf)$coefficients
  el$cof_bf <- cof_bf[grep("Intercept",names(cof_bf),invert =T)]#slope after
  el$cof_bf[is.na(el$cof_bf)] <- 0
  el$cof_af[is.na(el$cof_af)] <- 0
  
  el$height_p1  <- tempx_mean[p_af[2]] # Peak Height before
  el$height_m1  <- tempx_mean[p_bf[2]] # Peak Height before
  el$height_p2  <- tempx_mean[p_af[3]] # Peak Height after
  el$height_m2  <- tempx_mean[p_bf[3]] # Peak Height after
  
  el$height_p1  <- tempx_mean[p_af[2]]/el$q0.01 # relative Peak Height before
  el$height_m1  <- tempx_mean[p_bf[2]]/el$q0.01 # relative Peak Height before
  el$height_p2  <- tempx_mean[p_af[3]]/el$q0.01 # relative Peak Height after
  el$height_m2  <- tempx_mean[p_bf[3]]/el$q0.01 # relative Peak Height after
  
  
  
  
  # surrounding peaks
  # Peak Classes:
  el$n <- n
  el$v <- v
  el$assigned <- assigned
  # el <<- el
  # additional info
  el$X <- x[2]
  el$Xb <- x[3]
  el$Xs<- x[4]
  
  if(all(lengths(el)==1)){
    return(as.data.frame(el))
  }else{
    return(NA)
  }
}
findpeaks_Set <- function(x,threshold=0.1){
  x <<- x
  
  FiPeaks <- lapply(seq(0,1,by = 0.1),function(th){
    fi <- data.table(findpeaks(x,minpeakdistance = 1,nups =1,ndowns = 1,minpeakheight = 0,threshold = th,zero="+"))
    fi$Threshold=th
    fi
  })
  
  FiPeaks <- rbindlist(FiPeaks[lengths(FiPeaks)>2])
  if(dim(FiPeaks)[1]>0){
    FiPeaks <- FiPeaks[,.SD[Threshold==max(Threshold)],V2]
  }else{
    NULL
  }
  
} 
lmIntx <- function(fit1, fit2, rnd=2) {
  #Substitution method
  #https://stackoverflow.com/questions/7114703/finding-where-two-linear-fits-intersect-in-r
  b1<- fit1$coefficient[1]  #y-int for fit1
  m1<- fit1$coefficient[2]  #slope for fit1
  b2<- fit2$coefficient[1]  #y-int for fit2
  m2<- fit2$coefficient[2]  #slope for fit2
  if(m1==m2 & b1==b2) {print("Lines are identical")
  } else if(m1==m2 & b1 != b2) {print("Lines are parallel")
  } else {
    x <- (b2-b1)/(m1-m2)      #solved general equation for x
    y <- m1*x + b1            #plug in the result
    data.frame(x=round(x, rnd), y=round(y, rnd))
  }
}

