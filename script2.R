



rda="feature_counts.rda"
glmBase="reference"
dbiAnno=""
attach(rda)

ct= as.data.frame(fCount$counts) %>%
  mutate_all(as.integer)

coldata=read.csv("coldata.csv",header=T,sep=",")

ct <- ct[,coldata$bamID]
if(all(coldata$bamID==colnames(ct))) colnames(ct)=coldata$sampleID

an=as.data.frame(fCount$annotation) 
ct=left_join(mutate(ct,GeneID=rownames(ct)),an,by="GeneID")
gname=read.table("geneName.txt",header=T,sep="\t")
colnames(gname)=c("GeneID","GeneName")

ct=left_join(ct,gname,by="GeneID") %>% 
  dplyr::select(GeneName,colnames(an),everything())

a=unlist(strsplit(colnames(ct[,8:ncol(ct)]),split=".",fixed=T))
group= factor(a[seq(1,length(a),by=2)])

y <- DGEList(counts=ct[,8:ncol(ct)],genes=ct[,1:2],group=group)
rownames(y)=y$genes$GeneID

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

cpm=as.data.frame(cpm(y))
cpm$GeneID=rownames(cpm)
cpm=left_join(cpm,as.data.frame(y$genes),by="GeneID") %>%
  dplyr::select(GeneID,GeneName,everything())

cpm=left_join(cpm,an,by="GeneID") 

fpkm=cpm %>% mutate_at(vars(coldata$sampleID), ~ round(./Length*1000,3))



pca=prcomp(t(cpm[,3:(2+nrow(coldata))]),rank. = 10)
pe = (pca$sdev^2)/sum( (pca$sdev^2) ) *100
pcaData=as.data.frame(pca$x[,1:3]) %>% mutate(sampleID=rownames(pca$x),Group=group)


pcaPlot= PCAplot(pcaData, pe=pe,title="PCA plot of all samples")


