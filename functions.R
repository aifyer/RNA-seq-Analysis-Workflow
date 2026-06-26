
Volplot=function(dg1,title="Differential expression",ylimit_vol=c(-1,30)){
  len=str_count(title)
  ggplot(dg1)+
    geom_point(aes(x = logFC, y = -log10(FDR), color = significant), size = 1.2,alpha=0.3)+
    geom_vline(xintercept = c(1,-1),linetype = 2)+
    theme_bw()+labs(title = title)+
    scale_colour_manual(values = c("black", "red"))+
    theme(plot.title = element_text(hjust=0.5,size=(45*20/len)),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          legend.text  = element_text(size=12),
          legend.title = element_text(size=12))+
    xlab(bquote(paste(log[2], "(fold change)", 
                      sep = ""))) + ylab(bquote(paste(-log[10], "FDR", 
                                                      sep = "")))+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    scale_y_continuous(limits = ylimit_vol)
}

require(RColorBrewer)
colList =c(brewer.pal(n=10,name="Paired"),
           brewer.pal(n=9,name="Set1"),
           brewer.pal(n=12,name="Set3"))

Heat_colormap<- colorRampPalette(c("turquoise", "black", "red"))(10000)

PCAplot=function(pcaData=pcaData,pe=pe, title="PCA plot"){
  plot_ly(
    x=pcaData$PC1,
    y=pcaData$PC2,
    z=pcaData$PC3,
    type="scatter3d",
    mode="markers",
    text=pcaData$sampleID,
    showlegend=TRUE,
    legendgroup=pcaData$Group,
    name=factor(pcaData$Group),
    color=factor(pcaData$Group),
    colors=colList[1:4]
  ) %>%
    layout(margin=list(t=50,b=100,l=0,r=20),
           autosize=T,
           title=title,
           scene=list(
             xaxis=list(title=sprintf('PC1 (%.1f%%)',
                                      pe[1]),hoverformat='.2f'),
             yaxis=list(title=sprintf('PC2 (%.1f%%)',
                                      pe[2]),hoverformat='.2f'),
             zaxis=list(title=sprintf('PC3 (%.1f%%)',
                                      pe[3]),hoverformat='.2f')),
           annotations=list(showarrow=F,
                            text=c(" "),
                            x=pcaData$PC1,
                            y=pcaData$PC2,
                            z=pcaData$PC3),
           legend=list(
             orientation="v",
             x=100,y=0.5
           )
    )%>%
    config(displayModeBar=T,modeBarButtons=F,setBackground="transparent")
}

Heatplot=function(m,title = "Top 25 DE-genes Heatmap in pairwise comparison"){
  plot_ly(x=colnames(m),y=rownames(m),
          z=m, type="heatmap",colors = Heat_colormap)%>%
    layout( title = title,
            font= list(size = 11),
            xaxis = list(
              tickfont = list(size = 11),
              tickangle = 35,
              hoverformat = '.2f',
              showticklabels = TRUE),
            yaxis = list(
              tickfont = list(size = 10),
              tickwidth = 0, ticklen = 0.1,
              hoverformat = '.2f',
              showticklabels = TRUE))
  
}

Goplot=function(godat,sample="Comparison"){
  nTerm=min(nrow(godat),20)
  ggplot(godat[1:nTerm,],aes(x=Count,y=reorder(Description,-p.adjust), fill=p.adjust)) +
    geom_bar(stat = "identity") +
    theme(plot.title = element_blank(),
          axis.title.x = element_text(size=18),
          axis.title.y = element_text(size=18),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),
          legend.text  = element_text(size=12),
          legend.title = element_text(size=12))+
    xlab("Gene Counts") + ylab(paste0("Enriched BP of ",sample) )
}

PairComp = function(fit=fit,cpm=cpm,contrast=c(-1,1,0,0,0),path="",
                    sample1="Group1",sample2="Group2",pre=""){
  
  
  qlf <- glmQLFTest(fit, contrast=contrast)
  tr = glmTreat(fit,contrast=contrast,lfc=log2(1))
  
  t=topTags(tr,n=Inf)
  sum = as.data.frame(summary(decideTests(tr))) %>% 
    dplyr::select(Var1,Freq) 
  colnames(sum)=c("Test(FC2 & FDR0.05)",paste0(sample1,"vs",sample2))
  
  write.csv(sum,file=paste0("DE_summary_",sample1,"vs",sample2,".csv"),sep=",",row.names = F)
  
  dt=as.data.frame(decideTests(tr));colnames(dt)="Sig"
  dt = dt %>% mutate(GeneID=rownames(dt),
                     significant_Test=case_when(Sig == 1 ~ "Up",
                                                Sig == 0 ~ "NotSig",
                                                Sig == -1 ~ "Down")
  )
  detable=left_join(t$table[2:7],dt[,2:3],by="GeneID")
  
  de=left_join(cpm,detable,by="GeneID") %>%
    dplyr::select(GeneID,GeneName,contains(sample1),contains(sample2),logFC,logCPM,PValue,FDR,significant_Test)
  
  
  cpm2= cpm %>% dplyr::select(contains(sample1),contains(sample2))
  
  pca=prcomp(t(cpm2),rank. = 10)
  pe = (pca$sdev^2)/sum( (pca$sdev^2) ) *100
  a=unlist(strsplit(rownames(pca$x),split=".",fixed=T))
  group= factor(a[seq(1,length(a),by=2)])
  
  pcaData=as.data.frame(pca$x[,1:3]) %>% mutate(sampleID=rownames(pca$x),Group=group)
  
  pcaPlot = PCAplot(pcaData=pcaData,pe=pe, title=paste0("PCA plot of ",pre,": ", sample1, " vs ",sample2))
  saveWidget(pcaPlot, file=paste0("pcaPlot_(",sample1,"vs",sample2,").html"))
  
  de= transform(de, significant=ifelse(FDR<=0.05,"yes","no"))
  vplot=Volplot(de,title=paste0("Gene expression of ",pre,": ", sample1," vs ",sample2),ylimit_vol=c(0,ceiling(max(-log10(de$FDR)))))
  ggsave(paste0("volPlot_(",sample1,"vs",sample2,").png"),plot=vplot,width = 20, height = 20, units = "cm")
  
  t25=topTags(qlf,n=25)
  de25=right_join(cpm,t25$table[,c(2:7)],by="GeneID") %>%
    dplyr::distinct(GeneName,.keep_all = T)
  row.names(de25)=de25$GeneID
  hdata= de25 %>% arrange(logFC) %>%
    dplyr::select(contains(sample1),contains(sample2)) %>%
    mutate_all(scale)
  
  m=as.matrix(hdata)
  hplot=Heatplot(m,title = paste0("Top 25 DE-genes of ",sample1, " vs " ,sample2 ))
  saveWidget(hplot, file=paste0("heatPlot_(",sample1,"vs",sample2,").html"))
  
  
}
