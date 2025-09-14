rm(list=ls())
#------构建phyloseq对象
library(ggClusterNet)
library(phyloseq)
library(tidyverse)
# 从文件读取
metadata = read.delim("./meta_data.tsv")
head(metadata) 
row.names(metadata) = metadata$SampleID
otutab = read.table("./otu_table.txt", header=T, row.names=1, sep="\t", comment.char="", 
                    stringsAsFactors = F,encoding = "utf-8",check.names = F)
head(otutab)
taxonomy = read.table("./taxonomy.tsv", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
head(taxonomy)
library(tidyverse)
library(tidyfst)
tax <- taxonomy%>% separate_dt(Taxon,
                               into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                               sep = ";",
                               remove = F) 
tax = tax[,c("Kingdom","Phylum","Class","Order","Family","Genus","Species")] %>% as.data.frame()

head(tax)
row.names(tax) = row.names(taxonomy)
library(ggtree)
tree = read.tree("./tree.nwk")
# 导入phyloseq(ps)对象
library(phyloseq)
head(metadata)
ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE), 
              tax_table(as.matrix(tax)),
              phy_tree(tree)
)
ps
map = sample_data(ps) %>% as.data.frame()
map$ID= map$SampleID
map$Group = map$Group2
sample_data(ps) = map 
saveRDS(ps,"./ps.rds")
ps0=readRDS("./ps.rds")
library(tidyverse)
library(circlize)
library(phyloseq)
library(ggtree)
library(randomForest)
library(ggClusterNet)
#-主题--颜色等------#-------
package.amp <- function(){
  library(phyloseq)
  library(tidyverse)
  library(ggClusterNet)
  library(EasyStat)
  library(fs)
  library(ggthemes)
  library(RColorBrewer)
  library(magrittr)
  library(MicrobiotaProcess)
  library(ggsignif)
  library(ggtree)
  library(ggtreeExtra)
  library(ggstar)
  library(MicrobiotaProcess)
  library(ggnewscale)
  library(grid)
}

dir.amp <- function(ps0,
                    smart = FALSE){
  result_path <- paste("./","/result_and_plot/",sep = "")
  fs::dir_create(result_path)

  if (smart) {
    tax.1 = c("Fungi",
              "fungi",
              "K__Fungi",
              "k__Fungi",
              "d__Fungi",
              "d__fungi",
              "d__Eukaryota",
              "Eukaryota",
              "K:Fungi",
              "k:Fungi",
              "d:Fungi",
              "d:fungi",
              "d:Eukaryota",
              "Eukaryota",
              "D:Eukaryota"
    )
    TFdir.f <- as.data.frame(table(
      phyloseq::tax_table(ps0)[,1]))[,2][as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,1] %in%
                                           tax.1] > 10
    
    tax.2 = c("Bacteria",
              "K__Bacteria",
              "k__Bacteria",
              "d__Bacteria",
              "k:Bacteria",
              "bacteria",
              "D:bacteria",
              "D__bacteria",
              "d__bacteria",
              "d:prokaryotes",
              "d__prokaryotes",
              "prokaryotes"
    )
    
    TFdir.b <- as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,2][as.data.frame(table(phyloseq::tax_table(ps0)[,1]))[,1] %in%
                                                                        tax.2 ] > 10
    
    if (length(TFdir.f) != 0) {
      print("ITS")
      res1path <- paste(result_path,"/Base_diversity_ITS",sep = "")
      id = tax.1
    }
    
    if (length(TFdir.b) != 0) {
      print("16s")
      res1path <- paste(result_path,"/Base_diversity_16s",sep = "")
      id = tax.2
    }
  }else {
    res1path <- paste(result_path,"/Diversity_Micro",sep = "")
    id = "No check"
  }
  
  
  fs::dir_create(res1path)
  return(list(res1path,id))
}

theme_my = function(ps = ps0) {
  
  mytheme1 = ggplot2::theme_bw() + ggplot2::theme(
    legend.position="right",
    legend.title =  ggplot2::element_blank(),
    legend.background= ggplot2::element_blank(),
    legend.key= ggplot2::element_blank(),
    plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =  ggplot2::element_text(size = 14,colour = "black"),
    axis.title.x = ggplot2::element_text(size = 14,colour = "black"),
    axis.text =  ggplot2::element_text(size = 14),
    axis.text.x =  ggplot2::element_text(colour = "black",size = 14),
    axis.text.y =  ggplot2::element_text(colour = "black",size = 14),
    legend.text =  ggplot2::element_text(size = 15)
  )
  
  mytheme2 = ggplot2::theme_bw() + ggplot2::theme(
    panel.background=  ggplot2::element_blank(),
    panel.grid=  ggplot2::element_blank(),
    legend.position="right",
    legend.title =  ggplot2::element_blank(),
    legend.background=  ggplot2::element_blank(),
    legend.key= ggplot2::element_blank(),
    plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
    axis.title.y =  ggplot2::element_text(size = 24,colour = "black"),
    axis.title.x = ggplot2::element_text(size = 24,colour = "black"),
    axis.text =  ggplot2::element_text(size = 20,colour = "black"),
    axis.text.x =  ggplot2::element_text(colour = "black",size = 14,angle = 90),
    axis.text.y =  ggplot2::element_text(colour = "black",size = 14),
    legend.text =  ggplot2::element_text(size = 15,colour = "black")
  )

  gnum <- unique(phyloseq::sample_data(ps)$Group) %>% length()
  
  if (gnum < 10 ) {
    colset1 <- RColorBrewer::brewer.pal(9,"Set1")
  } else {
    colset1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(gnum)
  }
 
  colset2 <- RColorBrewer::brewer.pal(12,"Paired")
  colset3 <- c(RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(9,"Pastel1"))
  colset4 = colset3
  
  
  return(list(mytheme1,mytheme2,colset1,colset2,colset3,colset4))
  
}

mytheme1.1 = ggplot2::theme_bw() + ggplot2::theme(
  panel.background=  ggplot2::element_blank(),
  panel.grid=  ggplot2::element_blank(),
  legend.position="right",
  legend.title =  ggplot2::element_blank(),
  legend.background= ggplot2::element_blank(),
  legend.key= ggplot2::element_blank(),
  plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =  ggplot2::element_text(colour = "black"),
  axis.title.x = ggplot2::element_text(),
  axis.text =  ggplot2::element_text(),
  axis.text.x =  ggplot2::element_text(),
  axis.text.y =  ggplot2::element_text(),
  legend.text =  ggplot2::element_text()
)

mytheme2.1 = ggplot2::theme_bw() + ggplot2::theme(
  panel.background=  ggplot2::element_blank(),
  panel.grid=  ggplot2::element_blank(),
  legend.position="right",
  
  legend.title =  ggplot2::element_blank(),
  legend.background=  ggplot2::element_blank(),
  legend.key= ggplot2::element_blank(),
  plot.title =  ggplot2::element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =  ggplot2::element_text(),
  axis.title.x = ggplot2::element_text(),
  axis.text =  ggplot2::element_text(),
  axis.text.x =  ggplot2::element_text(angle = 90),
  axis.text.y =  ggplot2::element_text(),
  legend.text =  ggplot2::element_text()
)

res = theme_my()
mytheme1 = res[[1]];mytheme2 = res[[2]]; colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
result<- dir.amp(ps0 = ps0)
res1path = result[[1]];res1path
id = result[[2]];id
ps0 <- ps0 %>%
  phyloseq::filter_taxa(function(x) sum(x ) > 0 , TRUE);ps0
ps = ps0
gnum = phyloseq::sample_data(ps)$Group %>%unique() %>% length()
gnum
axis_order = phyloseq::sample_data(ps)$Group %>%unique()
Top_micro = 300
Top = 10
jj = j = "Phylum"
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)
b = NULL
CK = unique(phyloseq::sample_data(ps)$Group)[1]
heatnum　=　30
lefsenum = 0
ps_lefse <- ps %>%
  phyloseq::subset_taxa(
  )

ps_lefse = ggClusterNet::filter_OTU_ps(ps = ps_lefse,Top = 400)
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)
ROC = FALSE
rfcv = T
optimal = 40
N = 200
zipi = FALSE
otupath = paste(res1path,"/OTU_ps_20250709/",sep = "");otupath
dir.create(otupath)

#--基本表格保存#----------
tabpath = paste(otupath,"/report_table/",sep = "")
dir.create(tabpath)
otu = as.data.frame(t(ggClusterNet::vegan_otu(ps0)))
head(otu)
FileName <- paste(tabpath,"/otutab.csv", sep = "")
write.csv(otu,FileName,sep = "")

tax = as.data.frame((ggClusterNet::vegan_tax(ps0)))
head(tax)
FileName <- paste(tabpath,"/tax.csv", sep = "")
write.csv(otu,FileName,sep = "")

ps0_rela  = phyloseq::transform_sample_counts(ps0, function(x) x / sum(x) );ps0_rela 

otu_norm = as.data.frame(t(ggClusterNet::vegan_otu(ps0_rela)))
FileName <- paste(tabpath,"/otutab_norm.csv", sep = "")
write.csv(otu_norm,FileName,sep = "")

otutax <- cbind(as.data.frame(t(ggClusterNet::vegan_otu(ps0_rela))),as.data.frame((ggClusterNet::vegan_tax(ps0_rela))))
FileName <- paste(tabpath,"/otutax_norm.csv", sep = "")
write.csv(otutax,FileName,sep = "")


for (i in 2: length(phyloseq::rank_names(ps0))) {
  psi  <- ggClusterNet::tax_glom_wt(ps = ps0,ranks = phyloseq::rank_names(ps0)[i])

  otu = as.data.frame(t(ggClusterNet::vegan_otu(psi)))
  FileName <- paste(tabpath,"/otutab",phyloseq::rank_names(ps0)[i],".csv", sep = "")
  write.csv(otu,FileName,sep = "")

  tax = as.data.frame((ggClusterNet::vegan_tax(ps0)))
  FileName <- paste(tabpath,"/tax",phyloseq::rank_names(ps0)[i],".csv", sep = "")
  write.csv(otu,FileName,sep = "")
  
  psi_rela  = phyloseq::transform_sample_counts(psi, function(x) x / sum(x) );psi_rela 
  otu_norm = as.data.frame(t(ggClusterNet::vegan_otu(psi_rela)))
  FileName <- paste(tabpath,"/otutab_norm",phyloseq::rank_names(psi)[i],".csv", sep = "")
  write.csv(otu_norm,FileName,sep = "")
  
  otutax <- cbind(as.data.frame(t(ggClusterNet::vegan_otu(psi_rela))),as.data.frame((ggClusterNet::vegan_tax(psi_rela))))
  FileName <- paste(tabpath,"/otutax_norm",phyloseq::rank_names(ps0)[i],".csv", sep = "")
  write.csv(otutax,FileName,sep = "")
}
#--alpha多样性#---------
alppath = paste(otupath,"/alpha/",sep = "")
dir.create(alppath)
alpha = function(otu = NULL,
                 tax = NULL,
                 map = NULL,
                 ps = NULL,
                 group = "Group",
                 inde="Shannon",
                 sampling = TRUE,
                 Plot = TRUE){
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  if (sampling == TRUE) {
    samplesize = min(phyloseq::sample_sums(ps))
    if (samplesize == 0) {
      print("0 number sequence of some samples")
      print("median number were used")
      ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
    } else{
      ps11  = phyloseq::rarefy_even_depth(ps,sample.size = samplesize)
    }
  } else if(sampling == FALSE){
    ps11 = ps
  }
  
  
  mapping = phyloseq::sample_data(ps11)
  ps11 = phyloseq::filter_taxa(ps11, function(x) sum(x ) >0 , TRUE); ps11
  head(mapping)
  colnames(mapping) = gsub(group,"AA", colnames(mapping))
  
  mapping$Group = mapping$AA
  mapping$Group = as.factor(mapping$Group)
  mapping$Group 
  
  count = as.data.frame(t(ggClusterNet::vegan_otu(ps11)))
  
  alpha=vegan::diversity(count, "shannon")
  
  x = t(count)
  head(x)

  Shannon = vegan::diversity(x)  
  Shannon
  Inv_Simpson <- vegan::diversity(x, index = "invsimpson")
  Inv_Simpson

  S <- vegan::specnumber(x);S  
  S2 = rowSums(x>0)
 Pielou_evenness <- Shannon/log(S)
  Simpson_evenness <- Inv_Simpson/S
  est <- vegan::estimateR(x)
  est <- vegan::estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  report = cbind(Shannon, Inv_Simpson, Pielou_evenness, Simpson_evenness,
                 Richness, Chao1,ACE) 
  head(report)
  index = merge(mapping,report , by="row.names",all=F)
  return(index)
}
index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )

#--多种组合alpha分析和差异分析出图

alp = alpha(ps = ps0,inde="Shannon",group = "Group",Plot = TRUE )
index= alp

head(index)

sel = c(match("Shannon",colnames(index)),match("Richness",colnames(index)),match("Pielou_evenness",colnames(index)))
sel2 = c(match("Inv_Simpson",colnames(index)),match("Chao1",colnames(index)),match("Simpson_evenness",colnames(index)))
sel3=c(match("ACE",colnames(index)),match("Chao1",colnames(index)),match("Simpson_evenness",colnames(index)))
data = cbind(data.frame(ID = 1:length(index$Group),group = index$Group),index[sel])
head(data)

result = EasyStat::MuiKwWlx2(data = data,num = c(3:5))

FileName <- paste(alppath,"/alpha_diversity_different_label.csv", sep = "")
write.csv(result,FileName,sep = "")
FileName <- paste(alppath,"/alpha_diversity_index.csv", sep = "")
write.csv(index,FileName,sep = "")

result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3:5),
                                          result = result,
                                          sig_show ="abc",ncol = 3 )
p1_1 = result1[[1]] + 
  ggplot2::scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)
p1_1

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 3)
p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) +
  mytheme1+ 
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_2

res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:5),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) + 
  mytheme1 + 
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_3


FileName <- paste(alppath,"Alpha_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = ((1 + gnum) * 3), height =4,limitsize = FALSE)


FileName <- paste(alppath,"Alpha_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = ((1 + gnum) * 3), height = 4,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = ((1 + gnum) * 3), height = 4,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_box", ".jpg", sep = "")
ggsave(FileName, p1_1, width = ((1 + gnum) * 3), height =4,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_bar", ".jpg", sep = "")
ggsave(FileName, p1_2, width = ((1 + gnum) * 3), height = 4,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_boxbar", ".jpg", sep = "")
ggsave(FileName, p1_3, width = ((1 + gnum) * 3), height = 4,limitsize = FALSE)

krusk1 = ggpubr::compare_means( Shannon ~ group, data=data, method = "kruskal.test")
krusk2 = ggpubr::compare_means( Richness ~ group, data=data, method = "kruskal.test")
krusk3 = ggpubr::compare_means( Pielou_evenness ~ group, data=data, method = "kruskal.test")
krusk4 = ggpubr::compare_means( Inv_Simpson ~ group, data=data, method = "kruskal.test")
krusk5 = ggpubr::compare_means( Simpson_evenness ~ group, data=data, method = "kruskal.test")
krusk6 = ggpubr::compare_means( Chao1 ~ group, data=data, method = "kruskal.test")
krusk7 = ggpubr::compare_means( ACE ~ group, data=data, method = "kruskal.test")
dat = rbind(krusk1,krusk2,krusk3,krusk4,krusk5,krusk6,krusk7) %>% as.data.frame()
FileName <- paste(alppath,"/alpha_diversity_index_all_p_Kruskal-Wallis.csv", sep = "")
write.csv(dat,FileName,sep = "")
rare <- mean(phyloseq::sample_sums(ps))/10
aov.out = aov(Shannon ~ group, data = data)
summary(aov.out)
#---排序分析beta-diversity#-------
betapath = paste(otupath,"/bray/",sep = "")
dir.create(betapath)
BetaDiv = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,
                   group = "Group", dist = "bray", method ="PCoA",
                   Micromet = "adonis", pvalue.cutoff = 0.05,pair=TRUE){
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )

  if (method == "DCA") {
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    points = ordi$rproj[,1:2]
    colnames(points) = c("x", "y") 
    eig = ordi$evals^2
  }

  if (method == "CCA") {
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    points = ordi$CA$v[,1:2]
    colnames(points) = c("x", "y") 
    eig = ordi$CA$eig^2
  }

  if (method == "RDA") {
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    points = ordi$CA$v[,1:2]
    colnames(points) = c("x", "y") 
    eig = ordi$CA$eig
  }

  if (method == "MDS") {
    ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    points = ordi$vectors[,1:2]
    colnames(points) = c("x", "y") 
    eig = ordi$values[,1]
  }

  if (method == "PCoA") {
    unif = phyloseq::distance(ps1_rela , method=dist, type="samples")
    pcoa = stats::cmdscale(unif, k=2, eig=T) 
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y") 
    eig = pcoa$eig
  }
   otu_table = as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela )))
  if (method == "PCA") {
    otu.pca = stats::prcomp(t(otu_table), scale.default = TRUE)
    points = otu.pca$x[,1:2]
    colnames(points) = c("x", "y") 
    eig=otu.pca$sdev
    eig=eig*eig
  }
  if (method == "LDA") {
    data = t(otu_table)
    data = as.data.frame(data)
    data = scale(data, center = TRUE, scale = TRUE)
    dim(data)
    data1 = data[,1:10]
    map = as.data.frame(sample_data(ps1_rela))
    model = MASS::lda(data, map$Group)
    ord_in = model
    axes = c(1:2)
    points = data.frame(predict(ord_in)$x[, axes])
    colnames(points) = c("x", "y") 
    eig= ord_in$svd^2
  }

  if (method == "NMDS") {
  ordi = phyloseq::ordinate(ps1_rela, method=method, distance=dist)
    points = ordi$points[,1:2]
    colnames(points) = c("x", "y") 
    stress = ordi$stress
    stress= paste("stress",":",round(stress,2),sep = "")
  }
 if (method == "t-sne") {
    data = t(otu_table)
    data = as.data.frame(data)
    data = scale(data, center = TRUE, scale = TRUE)
    
    dim(data)
    map = as.data.frame(sample_data(ps1_rela))
    row.names(map)

    tsne = Rtsne::Rtsne(data,perplexity = 3)

    points = as.data.frame(tsne$Y)
    row.names(points) =  row.names(map)
    colnames(points) = c("x", "y") 
    stress= NULL
  }
  
  g = sample_data(ps)$Group %>% unique() %>% length()
  n = sample_data(ps)$Group%>% length()
  o = n/g
  
  if (o >= 3) {
    title1 = MicroTest(ps = ps1_rela, Micromet = Micromet, dist = dist)
    title1
  } else{
    title1 = NULL
  }
  
  if (pair==TRUE) {
    pairResult = pairMicroTest(ps = ps1_rela, Micromet = Micromet, dist = dist)
    
  } else {
    pairResult = "no result"
  }
  map = as.data.frame(phyloseq::sample_data(ps1_rela))
  map$Group = as.factor(map$Group)
  colbar = length(levels(map$Group))
  
  points = cbind(points, map[match(rownames(points), rownames(map)), ])
  points$ID = row.names(points)
  
  
  if (method %in% c("DCA", "CCA", "RDA",  "MDS", "PCoA","PCA","LDA")) {
    p2 =ggplot(points, aes(x=x, y=y, fill = Group)) +
      geom_point(alpha=.7, size=5, pch = 21) +
      labs(x=paste0(method," 1 (",format(100*eig[1]/sum(eig),digits=4),"%)"),
           y=paste0(method," 2 (",format(100*eig[2]/sum(eig),digits=4),"%)"),
           title=title1) +
      stat_ellipse(linetype=2,level=0.68,aes(group=Group, colour=Group))
    
    p3 = p2+ggrepel::geom_text_repel(aes(label=points$ID),size = 5)
    p3
  }
  
  if (method %in% c("NMDS")) {
    p2 =ggplot(points, aes(x=x, y=y, fill = Group)) +
      geom_point(alpha=.7, size=5, pch = 21) +
      labs(x=paste(method,"1", sep=""),
           y=paste(method,"2",sep=""),
           title=stress)+
      stat_ellipse( linetype = 2,level = 0.65,aes(group  =Group, colour =  Group))
    
    p3 = p2+ggrepel::geom_text_repel( aes(label=points$ID),size=4)
    p3
    if (method %in% c("t-sne")) {
      supp_lab = labs(x=paste(method,"1", sep=""),y=paste(method,"2",sep=""),title=title)
      p2 = p2 + supp_lab
      p3 = p3 + supp_lab
    }
    eig = NULL
  }
  
  if (method %in% c("t-sne")) {
    p2 =ggplot(points, aes(x=x, y=y, fill = Group)) +
      geom_point(alpha=.7, size=5, pch = 21) +
      labs(x=paste(method,"1", sep=""),
           y=paste(method,"2",sep=""))+
      stat_ellipse( linetype = 2,level = 0.65,aes(group  =Group, colour =  Group))
    
    p3 = p2+ggrepel::geom_text_repel( aes(label=points$ID),size=4)
    p3
    if (method %in% c("t-sne")) {
      supp_lab = labs(x=paste(method,"1", sep=""),y=paste(method,"2",sep=""),title=title1)
      p2 = p2 + supp_lab
      p3 = p3 + supp_lab
    }
    eig = NULL
  }
  return(list(p2,points,p3,pairResult,title1,eig))
}

MicroTest = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,
                     group = "Group", Micromet = "MRPP", dist = "bray"){
  
  
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
 ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela
  map = as.data.frame(phyloseq::sample_data(ps1_rela))
  unif = phyloseq::distance(ps, method=dist)
  if (Micromet == "adonis") {
    ado =  vegan::adonis2(unif ~ map$Group,permutations = 999)
    R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
    p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }

  if (Micromet == "MRPP") {
    dat.mrpp = vegan::mrpp(unif, map$Group)
    a = round(dat.mrpp$delta,3)
    R2 = paste("MRPP.delta ",a, sep = "")
    p_v = paste("p: ",round(dat.mrpp$Pvalue,3), sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }

  if (Micromet == "anosim") {
    dat.ano = vegan::anosim(unif, map$Group)
    a = round(dat.ano$statistic,3)
    R2 = paste("ANOSIM.r ",a, sep = "")
    p_v = paste("p: ",round(dat.ano$signif,3), sep = "")
    title1 = paste(R2," ",p_v, sep = "")
    title1
  }
  return(title1)
}
pairMicroTest = function(ps = ps, Micromet = "adonis", dist = "bray"){
  
  ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela
  
  map = as.data.frame(phyloseq::sample_data(ps1_rela)) 
  
  otu = ps1_rela %>% ggClusterNet::vegan_otu() %>%
    as.data.frame()
  map$Group = as.factor(map$Group)
  aa = levels(map$Group)
  aaa = combn(aa,2)
  aaa
  dim(aaa)[2]
  ID = rep("a",dim(aaa)[2])
  R = rep("a",dim(aaa)[2])
  P = rep("a",dim(aaa)[2])
  
  for (i in 1:dim(aaa)[2]) {
    print(i)
    Desep_group = aaa[,i]
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$ID = row.names(map)
    maps <- map[map$Group %in%Desep_group,]
    row.names(maps) = maps$ID
    ps_sub = ps1_rela
    phyloseq::sample_data( ps_sub ) = maps
    ps_sub = phyloseq::filter_taxa(ps_sub, function(x) sum(x ) > 0 , TRUE);ps_sub
    map = as.data.frame(phyloseq::sample_data(ps_sub))
    gg =map$Group
    unif <- phyloseq::distance(ps_sub, method=dist)
    
    if (Micromet == "MRPP") {
      mrpp = vegan::mrpp(unif, map$Group)
      as1 = round(mrpp$delta,3)
      R2 <- paste("MRPP.delta ",as1, sep = "")
      R2
      p_v = paste("p: ",round(mrpp$Pvalue,3), sep = "")
      p_v
    }
    
    if (Micromet == "anosim") {
      dat.ano = vegan::anosim(unif, map$Group)
      a = round(dat.ano$statistic,3)
      R2 <- paste("ANOSIM.r ",a, sep = "")
      R[i] = R2
      p_v = paste("p: ",round(dat.ano$signif,3), sep = "")
      P[i] = p_v
    }
    if (Micromet == "adonis") {
      str(otu)
      str(dune.env)
      ado = vegan:: adonis2( unif~ map$Group,method = "bray", by = NULL)
      R2 = paste("Adonis:R ",round(ado$R2[1],3), sep = "")
      p_v = paste("p: ",ado$`Pr(>F)`[1], sep = "")
      R[i] = R2
      P[i] = p_v
    }
    ID[i] = paste(Desep_group[1],Desep_group[2],sep = "_VS_")
    P[i] = p_v
    R[i] = R2
  }
  
  result = data.frame(ID = ID,stat = R,p = P)
  
  return(result)
}
methodlist = c("NMDS","PCoA", "PCA")

for (method in methodlist) {
  result = BetaDiv(ps = ps, group = "Group", dist = "bray",
                   method = method, Micromet = "anosim", pvalue.cutoff = 0.05,
                   pair = F)
  p3_1 = result[[1]] + 
    mytheme1 + 
    theme(legend.position = c(0.2,0.2))
  p3_1
  p3_2 = result[[3]] +
    mytheme1 + 
    theme(legend.position = c(0.2,0.2))
  p3_2
  
  FileName <- paste(betapath,"/a2_",method,"bray.pdf", sep = "")
  ggsave(FileName, p3_1, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"",method,"bray.jpg", sep = "")
  ggsave(FileName1 , p3_1, width = 12, height = 12)
  
  FileName <- paste(betapath,"/a2_",method,"bray_label.pdf", sep = "")
  ggsave(FileName, p3_2, width = 12, height = 12)
  FileName1 <- paste(betapath,"/a2_",method,"bray_label.jpg", sep = "")
  ggsave(FileName1 , p3_2, width = 12, height = 12)
  plotdata = result[[2]]
  FileName <-  paste(betapath,"/a2_",method,"bray.csv", sep = "")
  write.csv(plotdata,FileName)
  plotdata =result[[2]]
  head(plotdata)
  cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
  cent
  segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
                by = 'Group', sort = FALSE)
  library(ggsci)
  p3_3 = p3_1 +geom_segment(data = segs,
                            mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) +
    geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
    mytheme1 + 
    theme(legend.position = c(0.2,0.2))
  p3_3
  
  FileName <- paste(betapath,"/a2_",method,"bray_star.pdf", sep = "")
  ggsave(FileName, p3_3, width = 8, height = 8)
  FileName1 <- paste(betapath,"/a2_",method,"bray_star.jpg", sep = "")
  ggsave(FileName1 , p3_3, width = 8, height = 8)
  
  
}
p3_3 =ggplot(segs,aes(x = x, y = y))+
  geom_point(aes(shape=segs$geographic_emarcation,color = Group),size = 5)+
  scale_shape_binned()+
  mytheme1 + 
  theme(legend.position = c(0.2,0.2))
p3_3
map
TResult =result[[5]]
head(TResult)
pair = result[[4]]
pair
FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
write.csv(pair,FileName)
FileName <- paste(betapath,"Total_anosim.csv", sep = "")
write.csv(TResult,FileName)
title1 = MicroTest(ps = ps, Micromet = "adonis", dist = "bray")
title1
FileName <- paste(betapath,"Total_adonis.csv", sep = "")
write.csv(title1,FileName)
pairResult = pairMicroTest(ps = ps, Micromet = "adonis", dist = "bray")
FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
write.csv(pair,FileName)
#-------物种组成展示#----------
library(ggalluvial)
barMainplot = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",
                       j = "Phylum",axis_ord = NULL,label = TRUE ,sd = FALSE,Top = 10,tran = TRUE){
  
  
  if (is.null(axis_ord)) {
    axis_order = NA
  }else{
    axis_order = strsplit(basename(axis_ord), "-")[[1]]
  }
  
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  psdata <- ggClusterNet::tax_glom_wt(ps = ps,ranks = j)
  if (tran == TRUE) {
    psdata = psdata%>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x)} )
  }
  
  
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% 
    phyloseq::psmelt()
  
  
  Taxonomies$Abundance = Taxonomies$Abundance * 100
  colnames(Taxonomies) <- gsub(j,"aa",colnames(Taxonomies))
  data = c()
  i = 2
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    
    c <- Taxonomies %>% 
      dplyr::filter(Group == a)
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }

  Taxonomies = table
 by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance), sd(Abundance))
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  
  zhnagxu3 = zhnagxu2
 Taxonomies_x = plyr::ddply(zhnagxu3,"group", summarize,label_sd = cumsum(Abundance),label_y = cumsum(Abundance) - 0.5*Abundance)
  head( Taxonomies_x )
  Taxonomies_x = cbind(as.data.frame(zhnagxu3),as.data.frame(Taxonomies_x)[,-1])
  Taxonomies_x$label = Taxonomies_x$aa
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
 
  p4 <- ggplot(Taxonomies_x , aes(x =  group, y = Abundance, fill = aa, order = aa)) +
    geom_bar(stat = "identity",width = 0.5,color = "black") +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "")
  p4
  if (is.na(axis_order)) {
    p4 = p4
  }else{
    p4 = p4 + scale_x_discrete(limits = axis_order)
  }
  
  
  if (sd == TRUE) {
    p4 =  p4 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }
  
  if (label == TRUE) {
    p4 = p4 +
      geom_text(aes(y = label_y, label = label ))
  }
  
  
  map = as.data.frame(phyloseq::sample_data(ps))
  if (length(unique(map$Group))>3){p4 = p4 + theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}

  cs = Taxonomies_x $aa
  
  lengthfactor <- cs %>%
    levels() %>%
    length()
  cs4 <- cs %>%
    as.factor() %>%
    summary()  %>%
    as.data.frame()
  cs4$id = row.names(cs4)

  df_arrange<- dplyr::arrange(cs4, id)
  Taxonomies_x1<- dplyr::arrange(Taxonomies_x , aa)
  head(Taxonomies_x1)
  Taxonomies_x1$ID = factor(rep(c(1:lengthfactor), cs4$.))
  head(Taxonomies_x1)
  Taxonomies_x1$Abundance
  
  p3 <- ggplot(Taxonomies_x1, aes(x = group, y = Abundance,fill = aa,alluvium = aa,stratum = ID)) +
    ggalluvial::geom_flow(aes(fill = aa, colour = aa),
                          stat = "alluvium", lode.guidance = "rightleft",
                          color = "black",size = 0.2,width = 0.35,alpha = .2)  +
    geom_bar(width = 0.45,stat = "identity") +
    labs(x="",y="Relative abundance (%)",
         title= "") +
    guides(fill = guide_legend(title = j),color = FALSE) +
    scale_y_continuous(expand = c(0,0))
  p3

  p1 <- ggplot(Taxonomies_x1,
               aes(x = group, alluvium = aa, y = Abundance)) +
    ggalluvial::geom_flow(aes(fill = aa, colour = aa), width = 0)  +
    labs(x="",y="Relative abundance (%)",
         title="") +
    guides(fill = guide_legend(title = j),color = FALSE) +
    scale_y_continuous(expand = c(0,0))
  
  
  if (is.na(axis_order)) {
    p1 = p1
    p3 = p3
    
  }else{
    p1 = p1 + scale_x_discrete(limits = axis_order)
    p3 = p3 + scale_x_discrete(limits = axis_order)
    
  }
  if (label == TRUE) {
    p1 = p1 +
      geom_text(aes(y = label_y, label = label ))
    p3 = p3 +
      geom_text(aes(y = label_y, label = label ))
  }
  
  if (sd == TRUE) {
    p1 =  p1 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
    p3 =  p3 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }
  
  if (length(unique(map$Group))>3){	p3=p3+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
  
  return(list(p4,Taxonomies,p3,p1))
}

library(Hmisc)
library(plyr)
phyloseq::rank_names(ps0)

for (j in c("Phylum" , "Class" ,  "Order"  , "Family" , "Genus")) {
  result = barMainplot(ps = ps0,
                       j = j,
                       label = FALSE,
                       sd = FALSE,
                       Top = Top,
  )
  p4_1 <- result[[1]] + 
    scale_fill_manual(values = colset4) +
    scale_x_discrete(limits = axis_order) +
    mytheme1
  
  p4_1
  
  p4_2  <- result[[3]] + 
    scale_fill_manual(values = colset4) +
    scale_x_discrete(limits = axis_order) + 
    mytheme1
  p4_2
  
  databar <- result[[2]] 
  
  FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
  ggsave(FileName1, p4_2, width = (5+ gnum), height =8,limitsize = FALSE )
  FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
  ggsave(FileName2, p4_2, width = (5+ gnum), height =8,limitsize = FALSE )
  
  FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
  ggsave(FileName1, p4_1, width = (5+ gnum), height =8,limitsize = FALSE )
  FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
  ggsave(FileName2, p4_1, width = (5+ gnum), height =8,limitsize = FALSE )
  
  FileName <- paste(barpath,"/a2_",j,"_bar_data",".csv", sep = "")
  write.csv(databar,FileName)
}

detach("package:ggalluvial")
barpath = paste(otupath,"/circle_Micro_strack_bar/",sep = "");print(barpath)
dir.create(barpath)
circle_starc_bar = function(
    ps = ps,
    Top = 15,
    dist = "bray",
    cuttree = 3,
    hcluter_method = "complete"
    
){
  otu = ps %>% 
    ggClusterNet::scale_micro() %>%
    ggClusterNet::vegan_otu() %>% t() %>%
    as.data.frame()
  unif = phyloseq::distance(ps %>% ggClusterNet::scale_micro() , method = dist)
  hc <- stats::hclust(unif, method = hcluter_method )
  clus <- cutree(hc, cuttree )
  
  d = data.frame(label = names(clus), 
                 member = factor(clus))
  map = as.data.frame(phyloseq::sample_data(ps))
  
  dd = merge(d,map,by = "row.names",all = F)
  row.names(dd) = dd$Row.names 
  dd$Row.names = NULL
  p  = ggtree::ggtree(hc, layout='circular') %<+% dd + 
    geom_tippoint(size=5, shape=21, aes(fill= Group, x=x)) + 
    geom_tiplab(aes(color = Group,x=x * 1.2), hjust=1,offset=0.3) + xlim(-0.5,NA)
  p
  
  psdata =  ggClusterNet::tax_glom_wt(ps = ps %>% ggClusterNet::scale_micro(),ranks = "Phylum" )
  psdata = psdata%>% phyloseq::transform_sample_counts(function(x) {x/sum(x)} )
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  head(tax)
 j = "Phylum"
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "Other"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% phyloseq::psmelt()
  
  Taxonomies$Abundance = Taxonomies$Abundance * 100
  
  Taxonomies$OTU = NULL
  colnames(Taxonomies)[1] = "id"
  
  head(Taxonomies)
  
  dat2 = data.frame(id = Taxonomies$id,Abundance = Taxonomies$Abundance,Phylum = Taxonomies$Phylum)
  head(dat2)
  
  p2 <- p + 
    ggnewscale::new_scale_fill() + 
    ggtreeExtra::geom_fruit(
      data=dat2,
      geom=geom_bar,
      mapping=aes(
        x=Abundance, 
        fill= Phylum
      ),
      stat="identity",
      width = 0.4,
      orientation="y",
      offset=0.9,
      pwidth=2,
      axis.params=list(
        axis = "x", 
        text.angle = -45, 
        hjust = 0, 
        vjust = 0.5, 
        nbreak = 4
      )
    ) +
    scale_fill_manual(
      values=c(colset3),
      guide=guide_legend(keywidth=0.5,keyheight=0.5,order=5)
    ) +theme_void()
  
  return(p2)
}
library(ggtree)
j = "Phylum"
p2 = circle_starc_bar(
  ps = ps0,
  Top = 15,
  dist = "bray",
  cuttree = 3,
  hcluter_method = "complete")

FileName2 <- paste(barpath,"/a2_","_bar",".jpg", sep = "")
ggsave(FileName2, p2, width = 10, height =8 )

FileName2 <- paste(barpath,"/a2_","_bar",".pdf", sep = "")
ggsave(FileName2, p2, width = 10, height =8 )
cluMicro.bar <- function(
    dist = "bray",
    otu = NULL,
    tax = NULL,
    map = NULL,
    tree = NULL,
    j = "Phylum", 
    ps = ps,
    rep = 6 ,
    Top = 10, 
    tran = TRUE, 
    hcluter_method = "complete",
    Group = "Group",
    cuttree = 3){
  
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = Group)
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela)))
  unif = phyloseq::distance(ps1_rela , method = dist)
  hc <- stats::hclust(unif, method = hcluter_method )
  clus <- stats::cutree(hc, cuttree )
  d = data.frame(label = names(clus), 
                 member = factor(clus))
  map = as.data.frame(phyloseq::sample_data(ps))
  dd = merge(d,map,by = "row.names",all = F)
  row.names(dd) = dd$Row.names 
  dd$Row.names = NULL
  p  = ggtree::ggtree(hc) %<+% dd + 
    geom_tippoint(size=5, shape=21, aes(fill= Group, x=x)) + 
    geom_tiplab(aes(color = Group,x=x * 1.2), hjust=1)
  p
  psdata =  ggClusterNet::tax_glom_wt(ps = ps1_rela,ranks = j)
   if (tran == TRUE) {
    psdata = psdata %>% phyloseq::transform_sample_counts(function(x) {x/sum(x)} )
  }
   otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
   for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "Other"
    }
  }
  phyloseq::tax_table(psdata)= tax
   Taxonomies <- psdata %>% 
    phyloseq::psmelt()
  head(Taxonomies)
  Taxonomies$Abundance = Taxonomies$Abundance * 100
  
  Taxonomies$OTU = NULL
  colnames(Taxonomies)[1] = "id"
  head(Taxonomies)
  p <- p + ggnewscale::new_scale_fill()
  p
  p1 <- facet_plot(p, panel = 'Stacked Barplot', data = Taxonomies, geom = ggstance::geom_barh,mapping = aes(x = Abundance, fill = !!sym(j)),color = "black",stat='identity' )   
  p1
  grotax <- Taxonomies %>%
    dplyr::group_by(Group,!!sym(j)) %>%
    dplyr::summarise(Abundance = sum(Abundance))
  head(grotax)
  
  data = c()
  i = 2
  for (i in 1:length(unique(phyloseq::sample_data(psdata)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(psdata)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(psdata)$Group))[i,2]
    
    c <- grotax %>% 
      dplyr::filter(Group == a)
    c$Abundance <- c$Abundance/b
    head(c)
    data = c
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  sum( grotax$Abundance)
  head(table)
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) )
  otu = as.data.frame((ggClusterNet::vegan_otu(ps1_rela)))
  iris.split <- split(otu,as.factor(phyloseq::sample_data(ps1_rela)$Group))
  iris.apply <- lapply(iris.split,function(x)colMeans(x))
  iris.combine <- do.call(rbind,iris.apply)
  otuG = t(iris.combine)
  
  ps = phyloseq::phyloseq(phyloseq::otu_table(otuG,taxa_are_rows = T),
                          phyloseq::tax_table(ps1_rela)
                          
                          
  )
  hc = ps %>%
    phyloseq::distance(method = dist) %>%
    stats::hclust( method = hcluter_method )
  clus <- cutree(hc, cuttree )
  d = data.frame(label = names(clus), 
                 member = factor(clus))
  map = data.frame(ID = unique(phyloseq::sample_data(ps1_rela)$Group),row.names = unique(phyloseq::sample_data(ps1_rela)$Group),Group = unique(phyloseq::sample_data(ps1_rela)$Group))
  dd = merge(d,map,by = "row.names",all = F)
  row.names(dd) = dd$Row.names
  dd$Row.names = NULL
  p3  = ggtree(hc) %<+% dd + 
    geom_tippoint(size=5, shape=21, aes(fill= member, x=x)) + 
    geom_tiplab(aes(color = member,x=x * 1.2), hjust=1)
  p3
  p3 <- p3 + ggnewscale::new_scale_fill()
  head(grotax)
  
  p4 <- facet_plot(p3, panel = 'Stacked Barplot', data = table, geom = ggstance::geom_barh,mapping = aes(x = Abundance, fill = !!sym(j)),color = "black",stat='identity' )   
  p4
  
  return(list(p,p1,p3,p4,Taxonomies))
}

max_height <- 50
for (j in c("Phylum" , "Class" ,  "Order"  , "Family" , "Genus")) {
  
  result <-  cluMicro.bar (dist = "bray",
                           ps = ps,
                           j = j,
                           Top = Top, 
                           tran = TRUE,
                           hcluter_method = "complete",
                           Group = "Group",
                           cuttree = length(unique(phyloseq::sample_data(ps)$Group))
  )
  
  p5_1 <- result[[1]]
  p5_2 <- result[[2]]
  p5_3 <- result[[3]]
  p5_4 <- result[[4]]
  clubardata <- result[[5]]
  
  
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_sample",".pdf", sep = "")
  ggsave(FileName1, p5_1, width = 6, height = dim(phyloseq::sample_data(ps))[1]/4,dpi = 150, limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar_sample",".pdf", sep = "")
  ggsave(FileName1, p5_2, width = 12, height = dim(phyloseq::sample_data(ps))[1]/4 ,dpi = 150, limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_sample",".jpg", sep = "")
  ggsave(FileName1, p5_1, width = 6, height =dim(phyloseq::sample_data(ps))[1]/4 ,dpi = 150, limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar_sample",".jpg", sep = "")
  ggsave(FileName1, p5_2, width = 12, height =dim(phyloseq::sample_data(ps))[1]/4 ,dpi = 150, limitsize = FALSE)
  
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_Group",".pdf", sep = "")
  ggsave(FileName1, p5_3, width = 6, height = gnum,dpi = 150, limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar_Group",".pdf", sep = "")
  ggsave(FileName1, p5_4, width = 12, height = gnum ,dpi = 150, limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_Group",".jpg", sep = "")
  ggsave(FileName1, p5_3, width = 6, height = gnum ,dpi = 150, limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar_Group",".jpg", sep = "")
  ggsave(FileName1, p5_4, width = 12, height = gnum ,dpi = 150, limitsize = FALSE)
  
  FileName <- paste(barpath,"/a2_",j,"_cluster_bar_data",".csv", sep = "")
  write.csv(clubardata,FileName)
  
}
library(phyloseq)
library(ggClusterNet)
library(tidyverse)
map = sample_data(ps0) %>% as.data.frame()
map$ID= map$SampleID
map$Group = map$Group2
sample_data(ps0) = map 

library(tidyverse)
library(phyloseq)
res = theme_my(ps0)
mytheme1 = res[[1]]
mytheme2 = res[[2]] 
colset1 = res[[3]];
colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
otupath = paste(res1path,"/OTU_ps_20250709group/",sep = "");otupath
dir.create(otupath)
Micro_tern = function(ps = ps,
                      ternpath = ternpath
){
  ps_rela = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
  otu = ggClusterNet::vegan_otu(ps_rela) %>% as.data.frame()
  iris.split <- split(otu,as.factor(as.factor(phyloseq::sample_data(ps)$Group)))
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine) %>% as.data.frame()
  head(ven2)
  A <- combn(colnames(ven2),3)
  ven2$mean = rowMeans(ven2)
  tax = ggClusterNet::vegan_tax(ps)
  otutax = cbind(ven2,tax)
  head(otutax)
  otutax$Phylum[otutax$Phylum == ""] = "Unknown"
  for (i in 1:dim(A)[2]) {
    x = A[1,i]
    y = A[2,i]
    z = A[3,i]
    
    p <- ggtern::ggtern(data=otutax,aes_string(x = x,y=y,z=z,color = "Phylum",size ="mean" ))+
      geom_point() + 
      theme_void()
    p
    
    
    filename = paste(ternpath,"/",paste(x,y,z,sep = "_"),"_OTU.pdf",sep = "")
    ggsave(filename,p,width = 12,height = 10)
    
    filename = paste(ternpath,"/",paste(x,y,z,sep = "_"),"_OTU.png",sep = "")
    ggsave(filename,p,width = 12,height = 10,dpi = 100)
    filename = paste(ternpath,"/",paste(x,y,z,sep = "_"),"_OTU.csv",sep = "")
    write.csv(otutax,filename,quote = FALSE)
    
  }
  
}

library(ggtern)
ternpath = paste(otupath,"/ggtern/",sep = "")
dir.create(ternpath)

Micro_tern(ps = ps,ternpath = ternpath )
detach("package:ggtern")
ternpath = paste(otupath,"/ggtern/",sep = "")
dir.create(ternpath)

Micro_tern(ps = ps_lefse,ternpath = ternpath )
detach("package:ggtern")
VenSeper =  function(
    N= 0,
    otu = NULL,
    tax = NULL,
    map = NULL,
    tree = NULL,
    ps = NULL,
    group = "Group",
    path = "./",
    j = "Phylum",
    axis_ord = NULL,
    label = FALSE ,
    sd = FALSE,
    Top = 10
){
  Mytheme <- theme_bw()+
    theme(
      
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      
      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),
      legend.text = element_text(size = 15,face = "bold"),
      legend.position = "none"
      
    )
  ps = ggClusterNet::inputMicro(ps = ps, group  = "Group")
  mapping = as.data.frame(phyloseq::sample_data(ps))
  ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela 
  
  aa = ggClusterNet::vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  dim(count)
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  count[count > 0] <- 1
  count2 = as.data.frame(count)
  aa = sub_design[,"Group"]
  colnames(aa) = "Vengroup"
  iris.split <- split(count2,as.factor(aa$Vengroup))
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  for (i in 1:length(unique(phyloseq::sample_data(ps)[,"Group"]))) {
    aa <- as.data.frame(table(phyloseq::sample_data(ps)[,"Group"]))[i,1]
    bb =  as.data.frame(table(phyloseq::sample_data(ps)[,"Group"]))[i,2]
    ven2[,aa] = ven2[,aa]/bb
  }
  for (i in 1:4) {
    aa <- as.data.frame(table(phyloseq::sample_data(ps)[,"Group"]))[i,1]
    bb =  as.data.frame(table(phyloseq::sample_data(ps)[,"Group"]))[i,2]
    ven2[,aa] = ven2[,aa]/bb
  }
  
  ven2[ven2 < N]  = 0
  ven2[ven2 >=N]  = 1
  ven2 = as.data.frame(ven2)
  
  ven3 = as.list(ven2)
  for (i in 1:ncol(ven2)){
    
    ven3[[i]] <-  row.names(ven2[ven2[i] == 1,])
    
  }
  ven_pick = VennDiagram::get.venn.partitions(ven3)
  ven_pick = as.matrix(ven_pick)
  ven_pick = as.data.frame(ven_pick)
  colnames(ven_pick)
  row.names(ven_pick)
  
  
  sub_design <- as.data.frame(phyloseq::sample_data(ps1_rela))
  
  
  mapping = as.data.frame(phyloseq::sample_data(ps))
  
  
  single = as.data.frame(ven_pick$..set..)
  
  name = c()
  n = 1
  for (n in 1:length(ven_pick$..set..) ) {
    name[n] <- sapply(strsplit(basename(single[1,n]), ")"), `[`, 1)
  }
  
  signum = match(unique(mapping$Group),gsub("\\(","",name))
  length(ven_pick$..set..) 
  plot1 = list()
  plot2 = list()
  plot3 = list()
  i = 1
  for (i in 1:length(ven_pick$..set..) ) {
    ab = length(unique(mapping$Group))
    ab = 2+ab
    aab = ven_pick[[i,ab]]
    print(i)
    print(aab)
    if (length(aab) != 0) {
      otu = as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela)))
      dim(otu)
      otu$ID = row.names(otu)
      otu1<- dplyr::filter(as.tibble(otu), ID %in% aab)
      head(otu1)
      otu1 = as.data.frame(otu1)
      row.names(otu1) = otu1$ID
      
      otu1$ID = NULL
      subtab = as.matrix(otu1)
      ps_sub = ps1_rela
      ps_sub <- phyloseq::phyloseq(phyloseq::otu_table(subtab, taxa_are_rows=TRUE),
                                   phyloseq::tax_table(ps_sub),
                                   phyloseq::sample_data(ps_sub)
      )
      ps_sub
      num = ven_pick[i,1:length(unique(mapping$Group))]
      num = as.data.frame(num )
      colnames(num)[num[1,] == FALSE]
      map = as.data.frame(phyloseq::sample_data(ps_sub))
      map$ID = row.names(map)
      head(map)
      maps<- dplyr::filter(as.tibble(map),!Group %in% colnames(num)[num[1,] == FALSE])
      maps = as.data.frame(maps)
      row.names(maps) = maps$ID
      ps_sub = ps_sub
      phyloseq::sample_data( ps_sub ) = maps
      if (i %in%  signum) {
        tax = as.data.frame(ggClusterNet::vegan_tax(ps_sub))
        tax$count  =1
        pp <- ggplot(tax,aes(y = Phylum, x = count,fill = Genus)) + geom_bar(stat = "identity") +  theme_bw()+
          theme(
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.title = element_text(vjust = -8.5,hjust = 0.1),
            axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
            axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
            axis.text = element_text(size = 20,face = "bold"),
            axis.text.x = element_text(colour = "black",size = 14),
            axis.text.y = element_text(colour = "black",size = 14),
            legend.text = element_text(size = 15,face = "bold"),
          )
        FileName <- paste(Venpath,j,gsub("\\(","",name)[i],"OTU_number_of_each_part", ".pdf", sep = "")
        ggsave(FileName, pp, width = 30, height = 8)
      }
    }
    result = barMainplot(ps = ps_sub,j = j,axis_ord = axis_ord,label = label ,sd = sd,Top = Top,tran = FALSE)
    p = result[[1]] + scale_fill_brewer(palette = "Paired")
    result[[2]]
    p3 = result[[3]] + scale_fill_brewer(palette = "Paired")
    filename = paste(path,"/TAX_ven_pick_",ven_pick$..set..[[i]],".csv",sep = "")
    write.csv(result[[2]],filename,quote = FALSE)
vencount = as.data.frame(phyloseq::sample_sums(ps_sub ))
    end = merge(vencount ,sub_design ,by = "row.names",all = FALSE)
data_wt = data.frame(ID= end$Row.names,group = end$Group,count = end$`phyloseq::sample_sums(ps_sub)`)
    print("到这里了")
    head(data_wt)  
    colnames(data_wt)[3] =ven_pick$..set..[[i]]
    
    if (length(unique(data_wt$group)) !=1) {
      result = EasyStat::KwWlx(data = data_wt, i= 3)
      PlotresultBox = EasyStat::aovMuiBoxP(data = data_wt, i= 3,sig_show ="abc",result = result[[1]])
      PlotresultBox[[2]]
      p2 = PlotresultBox[[1]] + Mytheme
      p2
    }
    
    if (length(unique(data_wt$group)) ==1) {
      aa  =colnames(data_wt)[3]
      colnames(data_wt)[3] = "XX"
      
      p2 = ggplot(data_wt) + geom_boxplot(aes(x = group, y =XX,color = group))  + Mytheme  +
        labs(x="",
             y=aa)
      
    }
    filename = paste(path,"/SeqStat_ven_pick_",ven_pick$..set..[[i]],".csv",sep = "")
    write.csv(PlotresultBox[[2]],filename,quote = FALSE)
    
    plot1[[i]] =p
    plot2[[i]] =p2
    plot3[[i]] =p3
  } 

  pa  = ggpubr::ggarrange(plotlist = plot1, common.legend = TRUE, legend="right")
  pa
  
  pb  = ggpubr::ggarrange(plotlist = plot2, common.legend = TRUE, legend="right")
  pb
  
  pc  = ggpubr::ggarrange(plotlist = plot3, common.legend = TRUE, legend="right")
  pc
  return(list(pa,pb,pc))
}

library(ggalluvial)

barMainplot = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",
                       j = "Phylum",axis_ord = NULL,label = TRUE ,sd = FALSE,Top = 10,tran = TRUE){
  
  
  if (is.null(axis_ord)) {
    axis_order = NA
  }else{
    axis_order = strsplit(basename(axis_ord), "-")[[1]]
  }
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  psdata <- ggClusterNet::tax_glom_wt(ps = ps,ranks = j)
  if (tran == TRUE) {
    psdata = psdata%>%
      phyloseq::transform_sample_counts(function(x) {x/sum(x)} )
  }
  otu = phyloseq::otu_table(psdata)
  tax = phyloseq::tax_table(psdata)
  for (i in 1:dim(tax)[1]) {
    if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {
      
      tax[i,j] =tax[i,j]
    } else {
      tax[i,j]= "others"
    }
  }
  phyloseq::tax_table(psdata)= tax
  
  Taxonomies <- psdata %>% 
    phyloseq::psmelt()
  Taxonomies$Abundance = Taxonomies$Abundance * 100
  colnames(Taxonomies) <- gsub(j,"aa",colnames(Taxonomies))
  data = c()
  i = 2
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    a <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,1]
    
    b =  as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,2]
    
    c <- Taxonomies %>% 
      dplyr::filter(Group == a)
    c$Abundance <- c$Abundance/b
    data = data.frame(Sample =c$Sample,Abundance = c$Abundance,aa =c$aa,Group = c$Group)
    
    if (i == 1) {
      table = data
    }
    if (i != 1) {
      table = rbind(table,data)
    }
    
  }
  Taxonomies = table
  by_cyl <- dplyr::group_by(Taxonomies, aa,Group)
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance), sd(Abundance))
  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))
  head(zhnagxu2)
  colnames(zhnagxu2) <- c("aa","group","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = TRUE,levels = cc$aa)
  zhnagxu3 = zhnagxu2
  Taxonomies_x = plyr::ddply(zhnagxu3,"group", summarize,label_sd = cumsum(Abundance),label_y = cumsum(Abundance) - 0.5*Abundance)
  head( Taxonomies_x )
  Taxonomies_x = cbind(as.data.frame(zhnagxu3),as.data.frame(Taxonomies_x)[,-1])
  Taxonomies_x$label = Taxonomies_x$aa
  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = TRUE,levels = c(as.character(cc$aa)))
  p4 <- ggplot(Taxonomies_x , aes(x =  group, y = Abundance, fill = aa, order = aa)) +
    geom_bar(stat = "identity",width = 0.5,color = "black") +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)") +
    guides(fill = guide_legend(title = j)) +
    labs(x="",y="Relative abundance (%)",
         title= "")
  p4
  if (is.na(axis_order)) {
    p4 = p4
  }else{
    p4 = p4 + scale_x_discrete(limits = axis_order)
  }
  
  
  if (sd == TRUE) {
    p4 =  p4 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }
  
  if (label == TRUE) {
    p4 = p4 +
      geom_text(aes(y = label_y, label = label ))
  }
  map = as.data.frame(phyloseq::sample_data(ps))
  if (length(unique(map$Group))>3){p4 = p4 + theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
  cs = Taxonomies_x $aa
  lengthfactor <- cs %>%
    levels() %>%
    length()
  cs4 <- cs %>%
    as.factor() %>%
    summary()  %>%
    as.data.frame()
  cs4$id = row.names(cs4)
  df_arrange<- dplyr::arrange(cs4, id)
  Taxonomies_x1<- dplyr::arrange(Taxonomies_x , aa)
  head(Taxonomies_x1)
  Taxonomies_x1$ID = factor(rep(c(1:lengthfactor), cs4$.))
  head(Taxonomies_x1)
  Taxonomies_x1$Abundance
  
  p3 <- ggplot(Taxonomies_x1, aes(x = group, y = Abundance,fill = aa,alluvium = aa,stratum = ID)) +
    ggalluvial::geom_flow(aes(fill = aa, colour = aa),
                          stat = "alluvium", lode.guidance = "rightleft",
                          color = "black",size = 0.2,width = 0.35,alpha = .2)  +
    geom_bar(width = 0.45,stat = "identity") +
    labs(x="",y="Relative abundance (%)",
         title= "") +
    guides(fill = guide_legend(title = j),color = FALSE) +
    scale_y_continuous(expand = c(0,0))
  p3
  p1 <- ggplot(Taxonomies_x1,
               aes(x = group, alluvium = aa, y = Abundance)) +
    ggalluvial::geom_flow(aes(fill = aa, colour = aa), width = 0)  +
    labs(x="",y="Relative abundance (%)",
         title="") +
    guides(fill = guide_legend(title = j),color = FALSE) +
    scale_y_continuous(expand = c(0,0))
  
  
  if (is.na(axis_order)) {
    p1 = p1
    p3 = p3
    
  }else{
    p1 = p1 + scale_x_discrete(limits = axis_order)
    p3 = p3 + scale_x_discrete(limits = axis_order)
    
  }
  if (label == TRUE) {
    p1 = p1 +
      geom_text(aes(y = label_y, label = label ))
    p3 = p3 +
      geom_text(aes(y = label_y, label = label ))
  }
  
  if (sd == TRUE) {
    p1 =  p1 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
    p3 =  p3 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }
  
  if (length(unique(map$Group))>3){	p3=p3+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}
  
  return(list(p4,Taxonomies,p3,p1))
}

VenUpset = function(otu = NULL,
                    tax = NULL,
                    map = NULL,
                    tree = NULL,
                    ps = NULL,
                    group = "Group",
                    path = Venpath,
                    N = 0.5
){
  ps =  ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  aa =  ggClusterNet::vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  count[count > 0] <- 1
  count2 = as.data.frame(count )
  aa = sub_design[,"Group"]
  colnames(aa) = "Vengroup"
  iris.split <- split(count2,as.factor(aa$Vengroup))
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  for (i in 1:length(unique(phyloseq::sample_data(ps)[,"Group"]))) {
    aa <- as.data.frame(table(phyloseq::sample_data(ps)[,"Group"]))[i,1]
    bb =  as.data.frame(table(phyloseq::sample_data(ps)[,"Group"]))[i,2]
    ven2[,aa] = ven2[,aa]/bb
  }
  ven2[ven2 < N]  = 0
  ven2[ven2 >=N]  = 1
  ven2 = as.data.frame(ven2)
  ven3 = as.list(ven2)
filename = paste(path,"ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    pdf(file=filename3,width = 18, height = 15)
    T<- VennDiagram::venn.diagram(ven3,
                                  filename=NULL,
                                  lwd=2,
                                  lty=1,
                                  fill=c('red',"blue","yellow"), 
                                  col=c('red',"blue","yellow"), 
                                  cat.col=c('red',"blue","yellow"),
                                  cat.cex = 4,
                                  rotation.degree = 0,
                                  main = "",
                                  main.cex = 2,
                                  sub = "",
                                  sub.cex = 1,
                                  cex=3,
                                  alpha = 0.5,
                                  reverse=TRUE,
                                  scaled     = FALSE)
    grid::grid.draw(T)
    dev.off()
    filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    jpeg(file=filename33)
    grid::grid.draw(T)
    dev.off()
    grid::grid.draw(T)
  } 
  return(list(T,ven2))
group = "Group"
ps_Ven = ps
map = as.data.frame(phyloseq::sample_data(ps_Ven))
gnumven <- map[,group] %>% unique() %>% dim()
Venpath = paste(otupath,"/Ven_Upset_super50/",sep = "")
dir.create(Venpath)


result = VenUpset(ps = ps_Ven,
                  group = group,
                  path = Venpath
)

filename3 <- paste(Venpath,"Upset.pdf", sep = "")
pdf(file=filename3,width = 12, height = 12)


heatpath = paste(otupath,"/heapmap_boplot/",sep = "")
dir.create(heatpath)
Microheatmap <- function(ps_rela,
                         id,
                         label =  TRUE,
                         col_cluster = F,
                         row_cluster = TRUE
                         
){
  map = phyloseq::sample_data(ps_rela)
  map$ID = row.names(map)
  phyloseq::sample_data(ps_rela) = map
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  otu = as.matrix(otu[id,])
  ps_heatm = phyloseq::phyloseq(
    phyloseq::otu_table(otu,taxa_are_rows = T),
    phyloseq::tax_table(ps_rela),
    phyloseq::sample_data(ps_rela)
    
  )
  datah <- as.data.frame(t(ggClusterNet::vegan_otu(ps_heatm)))
  head(datah) 
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))
  otutaxh = cbind(datah,tax)
  head(otutaxh)
  otutaxh$id = paste(row.names(otutaxh),otutaxh$Genus,sep = "_")
  row.names(otutaxh) = otutaxh$id 
  data <- otutaxh[,c("id",phyloseq::sample_names(ps))]
  rig <- data[,phyloseq::sample_names(ps)] %>% rowMeans() %>% as.data.frame()
  head(rig)
  colnames(rig) = "MeanAbundance"
  rig$id = row.names(rig)
  rig = rig %>% dplyr::arrange(MeanAbundance)
  rig$id = factor(rig$id,levels = rig$id)
  p_rig = ggplot(rig) + geom_bar(aes(y = id,x = MeanAbundance),
                                 fill = "#A54657",
                                 stat = "identity") + theme_void()
  tem = data[,phyloseq::sample_names(ps)] %>% as.matrix()
  tem = scale(t(tem)) %>% t() %>% 
    as.data.frame()
  data[,phyloseq::sample_names(ps)] = tem
  mat <- data[,-1] 
  if (col_cluster ==  TRUE) {
    clust <- hclust(dist(mat %>% as.matrix())) 
    ggtree_plot <- ggtree::ggtree(clust)
  }
  if (row_cluster ==  TRUE) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }
  
  if (label ==  TRUE) {
    map = as.data.frame(phyloseq::sample_data(ps))
    map$ID = row.names(map)
    labels= ggplot(map, aes(x = ID, y=1, fill=Group)) + geom_tile() +
      theme_void()
  }
  map = phyloseq::sample_data(ps) %>% as.tibble() %>%
    dplyr::arrange(Group) %>% as.data.frame()
  map$ID
  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable,levels = map$ID)
  pcm$id = factor(pcm$id,levels = rig$id)
  p1 = ggplot(pcm, aes(y = id, x = variable)) + 
    geom_tile(aes(size = value,fill = value))+
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    scale_x_discrete(limits = rev(levels(pcm$variable)))  + 
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60))+
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)
    )
  colours = c( "#A54657",  "#582630", "#F7EE7F", "#4DAA57","#F1A66A","#F26157", "#F9ECCC", "#679289", "#33658A",
               "#F6AE2D","#86BBD8")
  p2 = ggplot(pcm, aes(y = id, x = variable)) + 
    geom_point(aes(size = value,fill = value), alpha = 0.75, shape = 21) + 
    scale_size_continuous(limits = c(0.000001, 100), range = c(2,25), breaks = c(0.1,0.5,1)) + 
    labs( y= "", x = "", size = "Relative Abundance (%)", fill = "")  + 
    scale_x_discrete(limits = rev(levels(pcm$variable)))  + 
    scale_y_discrete(position = "right")  +
    scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral")[11:1])(60)) +
    theme(
      panel.background=element_blank(),
      panel.grid=element_blank(),
      axis.text.x = element_text(colour = "black",angle = 90)
    )
  p1 <- p1  %>%
    aplot::insert_right(p_rig, width=.2) 
  p2 <- p2  %>%
    aplot::insert_right(p_rig, width=.2) 
  if (col_cluster ==  T) {
    p1 <- p1  %>%
      aplot::insert_left(ggtree_plot, width=.2)
    p2 <- p2  %>%
      aplot::insert_left(ggtree_plot, width=.2)
  }
  
  if (row_cluster ==  T) {
    p1 <- p1  %>%
      aplot::insert_top(labels, height=.02) 
    p2 <- p2  %>%
      aplot::insert_top(labels, height=.02) 
  }
  
  if (label ==  T) {
    p1 <- p1  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
    p2 <- p2  %>%
      aplot::insert_top(ggtree_plot_col, height=.1)
  }
  return(list(p1,p2))
}

for (j in 2:6) {
  
  ps_tem = ps %>% 
    ggClusterNet::scale_micro(method = "TMM") %>%
    ggClusterNet::tax_glom_wt(ranks = j) 
  rowSD = function(x){
    apply(x,1, sd)
  }
  
  rowCV = function(x){
    rowSD(x)/rowMeans(x)
  }
  
  id <- ps %>% 
    ggClusterNet::scale_micro(method = "TMM") %>%
    ggClusterNet::tax_glom_wt(ranks = j) %>%
    ggClusterNet::filter_OTU_ps(100) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(20) %>%
    names()
  result <- Microheatmap(ps_rela = ps_tem,id = id ,col_cluster = FALSE)
  p1 <- result[[1]] 
  p1
  p2 <- result[[2]]
  p2
  filename = paste(heatpath,"/",phyloseq::rank.names(ps)[j],"Topggheatmap.pdf",sep = "")
  ggsave(filename,p1,width = 14,height = (6 + heatnum/10))
  
  filename = paste(heatpath,phyloseq::rank.names(ps)[j],"Topggbubble.pdf",sep = "")
  ggsave(filename,p2,width = 14,height = (6 + heatnum/10))
  
  filename = paste(heatpath,"/",phyloseq::rank.names(ps)[j],"Topggheatmap.png",sep = "")
  ggsave(filename,p1,width = 14,height = (6 + heatnum/10))
  
  filename = paste(heatpath,phyloseq::rank.names(ps)[j],"Topggbubble.png",sep = "")
  ggsave(filename,p2,width = 14,height = (6 + heatnum/10))
}
#network
netpath = paste(otupath,"/network/",sep = "")
dir.create(netpath)
library(igraph)
library(sna)
library(phyloseq)
ps=ps0
result =ggClusterNet::network.2(ps = ps,
                                N = 500,
                                big = TRUE,
                                select_layout = TRUE,
                                layout_net = "model_maptree",
                                r.threshold=0.5,
                                p.threshold=0.05,
                                label = FALSE,
                                path = netpath,
                                zipi = F,
                                ncol = gnum,
                                nrow = 1,
                                fill = "Phylum"
)

p4_1 = result[[1]] + scale_fill_brewer(palette = "Paired") + mytheme1
data = result[[2]]
plotname1 = paste(netpath,"/network_all.pdf",sep = "")
ggsave(plotname1, p4_1,width = 2*gnum,height = 12,limitsize = FALSE)
tablename <- paste(netpath,"/co-occurrence_Grobel_net",".csv",sep = "")
write.csv(data,tablename)
p4_2 = result[[3]] + 
  scale_fill_brewer(palette = "Paired") +
  mytheme1
plotname1 = paste(netpath,"/network_all_cover.pdf",sep = "")
ggsave(plotname1, p4_2,width = 2*gnum,height = 24,limitsize = FALSE)

library(picante)
library(ape)
library(vegan)
library(FSA)
library(eulerr)
library(grid)
library(gridExtra)
require(minpack.lm)
require(Hmisc)
require(stats4)
library(parallel)

ps=ps0
psphy = filter_taxa(ps, function(x) sum(x ) > 2000 , TRUE);psphy 
ps=psphy
res1path = "./"
phypath = paste(res1path,"/Phylogenetic_analyse/",sep = "")
dir.create(phypath)

map = sample_data(ps)
n = map$Group %>% unique() %>%
  length()
n
library(ggClusterNet)
neutralModel = function(otu = NULL,
                        tax = NULL,
                        map = NULL,
                        tree = NULL,
                        ps = NULL,
                        group  = "Group",
                        ncol = 3,
                        nrow  = 1
                        
){
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  set.seed(72) 
  psrare = rarefy_even_depth(ps)
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))
  map = as.data.frame(sample_data(psrare))
  aa = levels(map$Group)
  aa
  map$ID = row.names(map)
  plots = list()
  dat1 = list()
  dat2 = list()
  i =1
  for (i in 1:length(aa)) {
    maps<- dplyr::filter(as.tibble(map),Group %in%aa[i])
    maps = as.data.frame(maps)
    row.names(maps) = maps$ID
    ps_sub = psrare
    sample_data( ps_sub ) =maps ;ps_sub
    OTU.table = t(otu_table(ps_sub))
    head(OTU.table )
    N <- mean(apply(OTU.table, 1, sum))
    p.m <- apply(OTU.table, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
    p.df = data.frame(p) %>%
      rownames_to_column(var="OTU")
    OTU.table.bi <- 1*(OTU.table>0)
    freq.table <- apply(OTU.table.bi, 2, mean)
    freq.table <- freq.table[freq.table != 0]
    freq.df = data.frame(OTU=names(freq.table), freq=freq.table)
    C <- inner_join(p.df,freq.df, by="OTU") %>%
      arrange(p)
    C.no0 <- C %>%
      filter(freq != 0, p != 0)
    d <- 1/N
    p.list <- C.no0$p
    freq.list <- C.no0$freq
    m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.1))
    m.ci <- confint(m.fit, 'm', level=0.95)
    m.sum <- summary(m.fit)
    m.coef = coef(m.fit)
    freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
    Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))
    fitstats <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2],
                           Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N,
                           Samples=nrow(OTU.table), Richness=length(p.list),
                           Detect=d)
    freq.pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)
    pred.df <- data.frame(metacomm_RA=p.list, frequency=freq.pred,
                          frequency_lowerCI=freq.pred.ci[,2],
                          frequency_upperCI=freq.pred.ci[,3]) %>%
      unique()
    obs.df = C.no0 %>%
      dplyr::rename(metacomm_RA = p, frequency=freq)
    head(obs.df)
    p = ggplot(data=obs.df) +
      geom_point(data=obs.df, aes(x=log10(metacomm_RA), y=frequency),
                 alpha=.3, size=2, color="#8DD3C7") +
      geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency), color="#FFFFB3") +
      geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_lowerCI), linetype=2, color="#FFFFB3") +
      geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_upperCI), linetype=2, color="#FFFFB3") +
      labs(x="Log10 abundance in\nmetacommunity", y="Frequency detected",title = paste(aa[i],paste("R^2 == ", round(fitstats$Rsqr, 3)),paste("italic(m) ==", round(fitstats$m, 3)))) +
      theme_bw() +
      theme(axis.line = element_line(color="black"),
            legend.position = "none",
            axis.title = element_text(size=14),
            axis.text = element_text(size=12))
    labs(
      x="Log10 abundance in\nmetacommunity", 
      y="Frequency detected",
      title = paste(aa[i], 
                    paste("R^2 == ", round(fitstats$Rsqr, 3)),
                    paste("italic(m) ==", round(fitstats$m, 3)))
    ) +
      theme_bw() +
      theme(
        axis.line = element_line(color="black"),
        legend.position = "none",
        axis.title = element_text(size=14),
        axis.text = element_text(size=12)
      )    
    
    p
    plots[[aa[i]]] = p
    dat1[[aa[i]]] = obs.df
    dat2[[aa[i]]] = pred.df
  }
  p  = ggpubr::ggarrange(plotlist = plots,common.legend = TRUE, legend="right",ncol = ncol,nrow = nrow)
  p
  return(list(p,plots,dat1,dat2))
}
result = neutralModel(ps = psphy,group  = "Group",ncol = 4)
p1 =  result[[1]]
p1
FileName <- paste(phypath,"1_neutral_modelCul", ".pdf", sep = "")
ggsave(FileName, p1,width = 12,height = 4)
FileName <- paste(phypath,"1_neutral_modelCul", ".png", sep = "")
ggsave(FileName, p1,width = 12,height = 4)
phyloSignal = function(otu = NULL,
                       tax = NULL,
                       map = NULL,
                       tree = NULL ,
                       ps = NULL,
                       env = env,
                       group  = "Group",
                       path = "./"){
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  set.seed(72)
  psrare = rarefy_even_depth(ps)
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))
  map = as.data.frame(sample_data(psrare))
  mapE =merge(map,env,by = "row.names",all= TRUE)
  row.names(mapE) = mapE$Row.names
  mapE$Row.names = NULL
  mapE$ID = row.names(mapE)
  sample_data(ps.norm) = mapE
  aa = levels(mapE$Group)
  dir.create(path)
  eco = "Endosp."
  for (eco in as.character(unique(mapE$Group))){
    print(paste("Now running", eco))
    sub.physeq = ps.norm
    otu = as.data.frame(vegan_otu(ps.norm))
    head(otu)
    map = as.data.frame(sample_data(ps.norm))
    mapsub <- map[map$Group == eco,]
    sample_data(sub.physeq) = mapsub
    OTU.table = otu_table(sub.physeq)
    OTU.table[OTU.table > 0] = 1
    OTU.freq = rowSums(OTU.table)
    OTU.freq = OTU.freq[OTU.freq > 2]
    sub.physeq = prune_taxa(names(OTU.freq), sub.physeq)
    sub.physeq
    tree = phy_tree(sub.physeq)
    phylo.dist = cophenetic(tree)
    sample_OTUs = tree$tip.label
    sam.phylo.dist = phylo.dist[sample_OTUs, sample_OTUs]
    sam.phylo.dist[upper.tri(sam.phylo.dist, diag=TRUE)] = NA
    site.chem.mat =  env[row.names(env) %in% row.names(mapsub),]
    site.chem.mat = as.matrix(site.chem.mat)
    otu.table = t(otu_table(sub.physeq))
    match(row.names(otu.table),row.names(site.chem.mat))
    OTU.niche = wascores(site.chem.mat, otu.table)
    OTU.niche.df = data.frame(OTU.niche)
    head( OTU.niche.df)
    for (i in 1:dim(OTU.niche.df)[2]) {
      pH.pref = OTU.niche.df[[i]]
      names(pH.pref) = rownames(OTU.niche.df)
      pH.dist = as.matrix(dist(pH.pref), labels=TRUE)
      sam.pH.dist = pH.dist[sample_OTUs, sample_OTUs]
      sam.pH.dist[upper.tri(sam.pH.dist, diag=TRUE)] = NA
      sam.pH.crlg = mantel.correlog(sam.pH.dist, sam.phylo.dist)
      filename = paste(path,eco,colnames(OTU.niche.df[i]), "_crlg.rds", sep="_")
      saveRDS(sam.pH.crlg, file=filename)
    }
  }
}
phySigPlot = function(otu = NULL,
                      tax = NULL,
                      map = NULL,
                      tree = NULL,
                      ps = NULL,
                      group  = "Group",
                      env = env,
                      path = "./"){
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  mapE = as.data.frame(sample_data(ps))
  for (eco in levels(mapE$Group)) {
    for (i in 1:length(colnames(env))) {
      ag.pH.crlg = data.frame(readRDS(file=paste(path,eco,colnames(env[i]), "_crlg.rds", sep="_"))$mantel.res) %>%
        mutate(Group = eco, property = colnames(env)[i])
      if (i == 1) {
        data = ag.pH.crlg
      }
      if (i != 1) {
        data = rbind(data,ag.pH.crlg )
      }
    }
    if (eco == levels(mapE$Group)[1]) {
      data2 =  data
    }
    if (eco != levels(mapE$Group)[1]) {
      data2 = rbind(data2, data)
    }
  }
  
  dim(data2)
  eco.crlg = data2 %>%
    mutate(sig = ifelse(Pr.corrected. <= 0.05, "significant", "non-significant")) %>%
    filter(!(is.na(Pr.corrected.)))
  eco.crlg$Group= factor(eco.crlg$Group)
  
  p = ggplot(data=eco.crlg, aes(x=class.index, y=Mantel.cor)) +
    geom_point(data=eco.crlg[eco.crlg$sig=="significant",], color = "black", size=2, shape=16) +
    geom_point(data=eco.crlg[eco.crlg$sig=="non-significant",], color = "black",size=2, shape=1) +
    geom_line(data=eco.crlg, aes(color=property)) +
    geom_hline(yintercept = 0, linetype=2) +
    labs(x = "Phylogenetic distance class", y="Mantel correlation", color="property") +
    facet_wrap(~Group,scales="free_y",ncol  = 4)
  return(list(p,eco.crlg,data2))
}
nullModel <- function(otu = NULL,
                      tax = NULL,
                      map = NULL,
                      tree = NULL ,
                      ps = NULL,
                      group  = "Group",
                      dist.method =  "bray",
                      gamma.method = "total",
                      transfer = "none",
                      null.model = "ecosphere"){
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  
  map = as.data.frame(sample_data(ps))
  grp1 = unique(map$Group)
  grp=list()
  for (i in 1:length(grp1)) {
    grp[[i]]=rownames(map)[which(map$Group==grp1[i])]
  }
  names(grp) = grp1
  
  report = c()
  dat4anova = c()
  grp4anova = c()
  report.ES = c()
  report.SES = c()
  otu = as.data.frame(t(vegan_otu(ps)))
  otu = as.matrix(otu)
  dataCK1 = otu[,grp[[x]]]
  if(gamma.method == "group"){
    rsum1 = rowSums(dataCK1)
    tempCK1 = which(rsum1==0)
    if(length(tempCK1)!=0) {dataCK1 = dataCK1[-tempCK1,]}
  }
   beta.dist = vegdist(t(dataCK1),method = dist.method)
  similarity.ob = 1 - beta.dist
  gamma = nrow(dataCK1)
  alpha = colSums(dataCK1>0)
  if(gamma.method == "group"){
    occur = apply(dataCK1, MARGIN=1, FUN=sum)
  }else{
    occur = apply(otu, MARGIN=1, FUN=sum)  
  }
 r = 100
  similarity.pm = matrix(0, nrow=ncol(dataCK1), ncol=ncol(dataCK1))
  similarity.pm = as.dist(similarity.pm)
  set.seed(123)  
  if(length(similarity.pm) >= 3) {
    sample_size <- min(5000, length(similarity.pm))
    sample_data <- sample(similarity.pm, sample_size)
    shapiro.test(sample_data)
  } else {
    print("no")
  }

  for(i in 1:r){
    PRM1 = matrix(0, ncol= ncol(dataCK1), nrow = nrow(dataCK1))
    if(null.model == "ecosphere"){
     for(j in 1:ncol(dataCK1)){
        aa = dataCK1[dataCK1[,j]>0,j]
        PRM1[sample(1:gamma, alpha[j], replace=FALSE, prob=occur), j] = aa
      }

    }else if(null.model == "ecosim"){
      PRM1 = randomizeMatrix(dataCK1, null.model="independentswap")
    }else if(null.model == "frequency"){
      PRM1 = randomizeMatrix(dataCK1, null.model="frequency")
    }
    dist_pm = vegdist(t(PRM1),method = dist.method)
    similarity.pm = similarity.pm + (1- dist_pm)
  }
  similarity.pm = similarity.pm/r
  normality = shapiro.test(similarity.pm)
  nor.p = normality$p.value
  ttest = t.test(similarity.pm, similarity.ob, alternative="two.sided", paired = TRUE, conf.level = 0.95)
  tt.p = ttest$p.value
  conf.int = ttest$conf.int
  pm.mean = mean(similarity.pm)
  pm.sd = sd(similarity.pm)
  ES = log(similarity.ob) - log(similarity.pm)
  effect.size = mean(ES)
  effect.size.sd = sd(ES)
  SES = (similarity.ob - similarity.pm)/pm.sd
  sd.effect.size = mean(SES)
  sd.effect.size.sd = sd(SES)
  ratio = 1 - similarity.pm / similarity.ob
  ratio.mean = mean(ratio)
  ratio.sd = sd(ratio)
  dat4anova = c(dat4anova, as.vector(ratio))
  grp4anova = c(grp4anova, rep(names(grp)[x], length(ratio)))
  
  conf.int.str = paste("[",paste(signif(conf.int,digits=3),collapse="~"),"]",sep="")
  report = rbind(report, c(mean(similarity.ob),sd(similarity.ob), pm.mean,  pm.sd, conf.int.str, nor.p, tt.p , effect.size, effect.size.sd, sd.effect.size, sd.effect.size.sd, ratio.mean, ratio.sd))
  report.ES = c(report.ES, effect.size)
  report.SES = c(report.SES, sd.effect.size)
}
for(x in c(1:length(grp))){  
  dataCK1 = otu[,grp[[x]]]
  if(gamma.method == "group"){
    rsum1 = rowSums(dataCK1)
    tempCK1 = which(rsum1==0)
    if(length(tempCK1)!=0) {dataCK1 = dataCK1[-tempCK1,]}
  }
  beta.dist = vegdist(t(dataCK1),method = dist.method)
  similarity.ob = 1 - beta.dist
  gamma = nrow(dataCK1)
  alpha = colSums(dataCK1>0)
  if(gamma.method == "group"){
    occur = apply(dataCK1, MARGIN=1, FUN=sum)
  }else{
    occur = apply(otu, MARGIN=1, FUN=sum)  
  }
  r = 100
  similarity.pm = matrix(0, nrow=ncol(dataCK1), ncol=ncol(dataCK1))
  similarity.pm = as.dist(similarity.pm)
  for(i in 1:r){
     PRM1 = matrix(0, ncol= ncol(dataCK1), nrow = nrow(dataCK1))
    if(null.model == "ecosphere"){
      for(j in 1:ncol(dataCK1)){
       aa = dataCK1[dataCK1[,j]>0,j]
        PRM1[sample(1:gamma, alpha[j], replace=FALSE, prob=occur), j] = aa
      }
    }else if(null.model == "ecosim"){
      PRM1 = randomizeMatrix(dataCK1, null.model="independentswap")
    }else if(null.model == "frequency"){
      PRM1 = randomizeMatrix(dataCK1, null.model="frequency")
    }
     dist_pm = vegdist(t(PRM1),method = dist.method)
     similarity.pm = similarity.pm + (1- dist_pm)
  }
  similarity.pm = similarity.pm/r
   if(length(similarity.pm) >= 3) {
     if(length(unique(similarity.pm)) < 3) {
      print("警告：所有相似度值都相同，无法进行正态性检验")
      nor.p = NA
    } else {
      sample_size <- min(5000, length(similarity.pm))
      sample_data <- sample(similarity.pm, sample_size)
      normality <- tryCatch({
        shapiro.test(sample_data)
      }, error = function(e) {
        print(paste("正态性检验失败:", e$message))
        return(list(p.value = NA))
      })
      nor.p <- normality$p.value
    }
  } else {
    print("警告：样本量不足，无法进行正态性检验")
    nor.p <- NA
  }
  ttest <- tryCatch({
    t.test(similarity.pm, similarity.ob, alternative="two.sided", paired = TRUE, conf.level = 0.95)
  }, error = function(e) {
    print(paste("t检验失败:", e$message))
    return(list(p.value = NA, conf.int = c(NA, NA)))
  })
  tt.p = ttest$p.value
  conf.int = ttest$conf.int
  pm.mean = mean(similarity.pm)
  pm.sd = sd(similarity.pm)
  tryCatch({
    ES = log(similarity.ob) - log(similarity.pm)
    effect.size = mean(ES)
    effect.size.sd = sd(ES)
    SES = (similarity.ob - similarity.pm)/pm.sd
    sd.effect.size = mean(SES)
    sd.effect.size.sd = sd(SES)
    ratio = 1 - similarity.pm / similarity.ob
    ratio.mean = mean(ratio)
    ratio.sd = sd(ratio)
  }, error = function(e) {
    print(paste("效应量计算失败:", e$message))
    ES = NA
    effect.size = NA
    effect.size.sd = NA
    SES = NA
    sd.effect.size = NA
    sd.effect.size.sd = NA
    ratio = NA
    ratio.mean = NA
    ratio.sd = NA
  })
  dat4anova = c(dat4anova, as.vector(ratio))
  grp4anova = c(grp4anova, rep(names(grp)[x], length(ratio)))
  conf.int.str = ifelse(all(!is.na(conf.int)), 
                        paste("[", paste(signif(conf.int, digits=3), collapse="~"), "]", sep=""),
                        "NA")
  tryCatch({
    report = rbind(report, c(mean(similarity.ob), sd(similarity.ob), pm.mean, pm.sd, 
                             conf.int.str, nor.p, tt.p, effect.size, effect.size.sd, 
                             sd.effect.size, sd.effect.size.sd, ratio.mean, ratio.sd))
    report.ES = c(report.ES, effect.size)
    report.SES = c(report.SES, sd.effect.size)
  }, error = function(e) {
    print(paste("更新报告失败:", e$message))
    report = rbind(report, rep(NA, 13))
    report.ES = c(report.ES, NA)
    report.SES = c(report.SES, NA)
  })
}
rownames(report) = grp1
colnames(report) = c("Mean of observed similarity", "Standard deviation of observed similarity",
                     "Mean of permutated similarity", "Standard deviation of permutated similarity",
                     "95% Conf int of perm similarity", "Normality test (p) on Perm similarity",
                     "T test on Ob and Perm similarity", "Effect size (ES)", "SD of ES",
                     "Standardized effect size (SES)", "SD of SES", "Mean of Ratio", "SD of Ratio")
head(report)

rep = t(report)
head(rep)
if (length(unique(grp4anova)) > 1) {
  aov.re = aov(dat4anova ~ grp4anova)
} else {
  aov.re = NULL
}
ratio = data.frame(ratio = dat4anova,group = grp4anova)
return(list(rep,ratio,aov.re))
result <- nullModel(ps = psphy,
                    group="Group",
                    dist.method =  "bray",
                    gamma.method = "total",
                    transfer = "none",
                    null.model = "ecosphere"
)
nullModeltab <- result[[1]]
ratiotab <- result[[2]]
aovtab <- result[[3]]
FileName <- paste(phypath,"3_nullModeltab", ".csv", sep = "")
write.csv(nullModeltab,FileName)
FileName <- paste(phypath,"3_ratiotab", ".csv", sep = "")
write.csv(ratiotab,FileName)
bNTICul = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",num = 99,thread = 1){
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  ps_sub <- ps
  map = as.data.frame(sample_data(ps_sub))
  map$ID = row.names(map)
  sample_data(ps) = map
  set.seed(72)  
  psrare = rarefy_even_depth(ps_sub)
  sample_sums(psrare)
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))
  bMNTD_null_func <- function(i, OTU.table, tree){
    tree$tip.label = sample(tree$tip.label)
    bMNTD_s = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
    A <- attr(bMNTD_s, "Size")
    B <- if (is.null(attr(bMNTD_s, "Labels"))) sequence(A) else attr(bMNTD_s, "Labels")
    if (isTRUE(attr(bMNTD_s, "Diag"))) attr(bMNTD_s, "Diag") <- FALSE
    if (isTRUE(attr(bMNTD_s, "Upper"))) attr(bMNTD_s, "Upper") <- FALSE
    bMNTD_s.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                            Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                            bMNTD = as.vector(bMNTD_s),
                            rep=i)
    return(bMNTD_s.df)
  }
  Phylo_turnover <- function(physeq, reps, nproc){
    OTU.table = t(otu_table(physeq))
    tree = phy_tree(physeq)
   bMNTD_o = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
    A <- attr(bMNTD_o, "Size")
    B <- if (is.null(attr(bMNTD_o, "Labels"))) sequence(A) else attr(bMNTD_o, "Labels")
    if (isTRUE(attr(bMNTD_o, "Diag"))) attr(bMNTD_o, "Diag") <- FALSE
    if (isTRUE(attr(bMNTD_o, "Upper"))) attr(bMNTD_o, "Upper") <- FALSE
    bMNTD_o.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                            Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                            bMNTD = as.vector(bMNTD_o))
    rep.list = seq(1, reps)
    bMNTD_s.df.list = mclapply(rep.list, bMNTD_null_func, OTU.table=OTU.table, tree=tree, mc.cores=nproc)
    bMNTD_s.df <- do.call("rbind", bMNTD_s.df.list)
    bMNTD_s.means.df = bMNTD_s.df %>%
      group_by(Sample_1, Sample_2) %>%
      dplyr::summarize(mean_bMNTD = mean(bMNTD),
                       sd_bMNTD = sd(bMNTD))
    bMNTD_o.df = inner_join(bMNTD_o.df, bMNTD_s.means.df, by=c("Sample_1", "Sample_2")) %>%
      mutate(bNTI = (bMNTD - mean_bMNTD)/sd_bMNTD)
    
    return(bMNTD_o.df)
  }
  bNTI = Phylo_turnover(psrare, num, thread)
  return(list(bNTI))
}
result = bNTICul(ps = psphy,group  = "Group",num = 99,thread = 1) 
bNTI = result[[1]]
head(bNTI)
filename = paste(phypath,"/4_bNTI.csv",sep = "")
write.csv(bNTI, filename)
RCbary = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",num = 99,thread = 1){
  ps_sub <- ps
  map = as.data.frame(sample_data(ps_sub))
  map$ID = row.names(map)
  sample_data(ps) = map
  set.seed(72)  
  psrare = rarefy_even_depth(ps_sub )
  sample_sums(psrare)
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))
  RCbray_null_func <- function(i, freq.abd.df, alpha1, alpha2, N){
    simcom1 = data.frame(table(sample(freq.abd.df$OTU, size=alpha1, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
    colnames(simcom1) = c("OTU","simcom1")
    simcom1$OTU = as.character(simcom1$OTU)
    simcom1 = inner_join(simcom1, freq.abd.df, by="OTU")
    simcom2 = data.frame(table(sample(freq.abd.df$OTU, size=alpha2, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
    colnames(simcom2) = c("OTU","simcom2")
    simcom2$OTU = as.character(simcom2$OTU)
    simcom2 = inner_join(simcom2, freq.abd.df, by="OTU")
    simcom1.abd = data.frame(table(sample(simcom1$OTU, size=N-alpha1, replace=T, prob=simcom1$p)), stringsAsFactors = F)
    colnames(simcom1.abd) = c("OTU","simcom1.abd")
    simcom1.abd$OTU = as.character(simcom1.abd$OTU)
    simcom1 = full_join(simcom1, simcom1.abd, by="OTU") %>%
      mutate(simcom1.abd = ifelse(is.na(simcom1.abd), 1, simcom1.abd)) %>%
      select(OTU, simcom1.abd)
    simcom2.abd = data.frame(table(sample(simcom2$OTU, size=N-alpha2, replace=T, prob=simcom2$p)), stringsAsFactors = F)
    colnames(simcom2.abd) = c("OTU","simcom2.abd")
    simcom2.abd$OTU = as.character(simcom2.abd$OTU)
    simcom2 = full_join(simcom2, simcom2.abd, by="OTU") %>%
      mutate(simcom2.abd = ifelse(is.na(simcom2.abd), 1, simcom2.abd)) %>%
      select(OTU, simcom2.abd)
    simcom = full_join(simcom1, simcom2, by="OTU")
    simcom[is.na(simcom)] = 0
    rownames(simcom) = simcom$OTU
    simcom$OTU = NULL
    
    null.dist = vegdist(t(simcom), method="bray")[1]
    return(null.dist)
  }
  Calc_RCbray <- function(physeq, reps, nproc){
    otu.table = otu_table(physeq)
    otu.PA.table = otu.table
    otu.PA.table[otu.PA.table > 0] = 1
    alpha.df = data.frame(Sample_ID = colnames(otu.PA.table), OTU.n = colSums(otu.PA.table), stringsAsFactors = F)
    beta.table = as.matrix(vegdist(t(otu.PA.table), method="bray", diag=TRUE, upper=TRUE))
    N <- mean(apply(t(otu.table), 1, sum))
    p.m <- apply(t(otu.table), 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
    otu.table.bi <- 1*(t(otu.table)>0)
    freq <- apply(otu.table.bi, 2, mean)
    freq <- freq[freq != 0]
    freq.abd.df = data.frame(p=p, freq=freq) %>%
      tibble::rownames_to_column(var="OTU") %>%
      filter(p != 0, freq != 0) %>%
      arrange(p)
    comps = combn(alpha.df$Sample_ID, m=2, simplify = F)
    RCb.df = data.frame(Site1 = character(), Site2 = character(), RCb = numeric(), stringsAsFactors = F)
    for (j in seq(1, length(comps))){
      sam = comps[[j]]
      alpha1 = alpha.df[alpha.df$Sample_ID == sam[1],]$OTU.n
      alpha2 = alpha.df[alpha.df$Sample_ID == sam[2],]$OTU.n
      rep.list = seq(1, reps)
      null.list = mclapply(rep.list, RCbray_null_func, freq.abd.df=freq.abd.df, alpha1=alpha1, alpha2=alpha2, N=N, mc.cores=nproc)
      
      RCb = (length(null.list[null.list > beta.table[sam[1], sam[2]]]) + (0.5*length(null.list[null.list == beta.table[sam[1], sam[2]]])))/reps
      RCb = (RCb - 0.5)*2
      
      RCb.df = rbind(RCb.df, data.frame(Site1=sam[1], Site2=sam[2], RCb=RCb, stringsAsFactors = F))
    }
    RCb.df
    return(RCb.df)
  }
  RCb = Calc_RCbray(psrare, num, thread)
  head(RCb)
  return(list(RCb))
}
result = RCbary(ps = psphy ,group  = "Group",num = 99,thread = 1)

RCbary = result[[1]]
head(RCbary)
filename = paste(phypath,"/5_RCb.csv",sep = "")
write.csv(RCbary,filename)
bNTIRCPlot = function(otu = NULL,tax = NULL,
                      map = NULL,tree = NULL ,
                      ps = NULL,
                      RCb  = RCb,bNTI = bNTI,group  = "Group"){
  
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  
  psrare <- ps
  map = as.data.frame(sample_data(psrare))
  map$ID = row.names(map)
  sample_data(psrare) = map
  eco.meta1 = data.frame(sample_data(psrare)) %>%
    select(ID, Group) %>%
    dplyr::rename(Sample_1 = ID, Group_1 = Group)
  
  eco.meta2=data.frame(sample_data(psrare)) %>%
    select(ID, Group) %>%
    dplyr::rename(Sample_2 = ID, Group_2 = Group)
  bNTI.df = inner_join(bNTI, eco.meta1) %>%
    inner_join(eco.meta2)
  turnover.df = inner_join(bNTI.df, RCb)
  head(turnover.df)
  dim(turnover.df)
  dim(bNTI.df)
  within.bNTI.df = bNTI.df %>%
    filter(Group_1 == Group_2) %>%
    mutate(Group = Group_1)
  
  head(within.bNTI.df )
  eco.bNTI.plot <- ggplot(within.bNTI.df, aes(x=Group, y=bNTI)) +
    geom_jitter(alpha = 0.1,color ="#984EA3") +
    geom_boxplot(outlier.shape=1,outlier.alpha = 0,fill = "#984EA3") +
    
    geom_hline(yintercept = 2, linetype=2, size=0.5) +
    geom_hline(yintercept = -2, linetype=2, size=0.5) +
    labs(x="", y="bNTI") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title = element_text(size=14))
  eco.turnover.df = turnover.df %>%
    filter(Group_1 == Group_2) %>%
    mutate(Group = Group_1)
  head(eco.turnover.df )
  eco.turnover.df = eco.turnover.df %>%
    mutate(process = ifelse(abs(bNTI) < 2,
                            ifelse(abs(RCb) < 0.95, "Drift",
                                   ifelse(RCb >= 0.95, "Dispersal Limited",
                                          ifelse(RCb <= -0.95, "Homogenizing Dispersal", "ERROR"))),
                            ifelse(bNTI >= 2, "Variable Selection",
                                   ifelse(bNTI <= -2, "Homogeneous Selection", "ERROR"))))
  eco.turnover.df$process = factor(eco.turnover.df$process, levels = c("Drift",
                                                                       "Dispersal Limited", "Homogenizing Dispersal",
                                                                       "Variable Selection", "Homogeneous Selection"))
  
  head(eco.turnover.df)
   pre = eco.turnover.df %>%
    dplyr::group_by(Group, process) %>%
    dplyr::summarize(n_sites = n(),
                     perc=(n()/45)*100) %>%
    as.data.frame
  numeco <- pre %>%  dplyr::group_by(Group) %>% 
    dplyr::summarise(num = sum(n_sites))
  alleco <- pre %>% dplyr::left_join(numeco,by = "Group")
  alleco$perc =  alleco$n_sites/ alleco$num * 100
  sum.eco.turnover.df = alleco
  eco.turnover.plot = ggplot(sum.eco.turnover.df, aes(x=Group, y=perc, fill=process)) +
    geom_bar(stat="identity", color="black") +
    labs(x="", y="Percent of pairs (%)", fill="Process") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=12),
          axis.text.x = element_text(angle=45, hjust=1),
          axis.title = element_text(size=14),
          legend.key.size = unit(10, "mm"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14))
  eco.turnover.plot
  eco.plot = cowplot::plot_grid(eco.bNTI.plot, eco.turnover.plot,
                                rel_widths=c(0.6, 1), labels=c("A", "B"))
  eco.plot
  
  
  return(list( eco.bNTI.plot, eco.turnover.plot,eco.plot,turnover.df,sum.eco.turnover.df))
}

bNTI = read.csv(paste(phypath,"/4_bNTI.csv",sep = ""),row.names = 1)
head(bNTI)
RCb = read.csv(paste(phypath,"/5_RCb.csv",sep = ""),row.names = 1) %>%
  dplyr::mutate(Sample_1 = Site2, Sample_2 = Site1)
head(RCb)
result = bNTIRCPlot(ps = psphy ,RCb  = RCb,bNTI = bNTI,group  = "Group")
p3 <- result[[1]] + mytheme1
p3
p4 <- result[[2]] + mytheme1
p4
p5 <- result[[3]]
p5
plotdata = result[[4]]
head(plotdata)

filename = paste(phypath,"/6_bNTI_RCbray.csv",sep = "")
write.csv(plotdata,filename)

FileName <- paste(phypath,"6_bNTI", ".pdf", sep = "")
ggsave(FileName, p3,width =8,height = 6)

FileName <- paste(phypath,"6_RCbary", ".pdf", sep = "")
ggsave(FileName, p4,width = 6,height = 6)

FileName <- paste(phypath,"6_BNTI_RCbray", ".pdf", sep = "")
ggsave(FileName, p5,width = 12,height = 8)

FileName <- paste(phypath,"6_bNTI", ".png", sep = "")
ggsave(FileName, p3,width =8,height = 6)

FileName <- paste(phypath,"6_RCbary", ".png", sep = "")
ggsave(FileName, p4,width = 6,height = 6)

FileName <- paste(phypath,"6_BNTI_RCbray", ".png", sep = "")
ggsave(FileName, p5,width = 12,height = 8)

#upset
library(VennDiagram)
library(UpSetR)
library(dplyr) 
library(ggplot2)
library(RColorBrewer) 
Data = read.delim("20250807.csv", header=T, row.names= 1, sep=",", stringsAsFactors = T, check.names = FALSE) 
design = read.delim("design.csv", header=T, row.names= 1, sep=",", stringsAsFactors = T, check.names = FALSE)
index = rownames(design) %in% colnames(Data) 
Data = Data[,rownames(design)] 
Data_t = t(Data)
Data_t2 = merge(design, Data_t, by="row.names")
Data_t2 = Data_t2[,c(-1,-3)]
cols_to_convert <- setdiff(names(Data_t2), names(Data_t2[1])) 
Data_t2[cols_to_convert] <- lapply(Data_t2[cols_to_convert], function(x) {
  as.numeric(as.character(x)) 
})
Data_t2[is.na(Data_t2)] <- 0
Data_mean = aggregate(Data_t2[,-1], by=Data_t2[1], FUN=mean)
Data4Pic = as.data.frame(do.call(rbind, Data_mean)[-1,])
group = Data_mean$Group
colnames(Data4Pic) = group
Data4Pic[Data4Pic>0]=1
write.table(Data4Pic,"datavenn.txt",sep = "\t")
pdf(file="Genus_venn2025008.pdf", height = 4, width = 6)
p1 <- venn.diagram(
  x=list(
    A=row.names(Data4Pic[Data4Pic$DE==1,]),
    B=row.names(Data4Pic[Data4Pic$ME==1,]),
    C=row.names(Data4Pic[Data4Pic$GE==1,]),
    D=row.names(Data4Pic[Data4Pic$GP==1,])),
  filename = NULL, 
  lwd = 3,
  alpha = 0.6,
  label.col = "white",
  cex = 1.5,
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
  fontfamily = "serif",
  fontface = "bold",
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  margin = 0.05)
p1
grid.draw(p1)
dev.off()
pdf(file="Genus_Upsetplot202508.pdf", height = 4, width = 6)
p2<- upset(Data4Pic, sets =colnames(Data4Pic) , order.by = "freq")
p2
dev.off()
pdf(file="Genus_Upsetplot202508.pdf", height = 10, width = 10)

pdf(file="Genus_Upsetplot_indiv202508.pdf",height = 4, width = 10)
p3<-upset(Data4Pic, sets = colnames(Data4Pic), mb.ratio = c(0.55, 0.45), order.by = "freq",
          queries = list(list(query=intersects, params=list("A", "B"), color="purple", active=T), 
                         list(query=intersects, params=list("C", "D", "A"), color="green", active=T), 
                         list(query=intersects, params=list("B", "C", "A", "D"), color="blue", active=T)), 
          nsets = 3, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Shared Genus",
          sets.x.label = "General Number in Each Group", text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5))
p3
dev.off()
library(UpSetR)
library(ggplot2)
library(grid)
library(plyr)
pdf(file="Genus_Upsetplot_attribute2.pdf",height = 8, width = 10)
p4<-upset(Data4Pic, sets = colnames(Data4Pic), mb.ratio = c(0.55, 0.45), order.by = "freq",
          queries = list(list(query=intersects, params=list("A", "B"), color="purple", active=T), 
                         list(query=intersects, params=list("C", "D", "A"), color="green", active=T), 
                         list(query=intersects, params=list("B", "C", "A", "D"), color="blue", active=T)), 
          nsets = 3, number.angles = 0, point.size = 4, line.size = 1, mainbar.y.label = "Number of Shared Genus",
          sets.x.label = "General Number in Each Group", text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5), 
          attribute.plots = list(gridrows = 45, plots = list(list(plot = scatter_plot,  x = "ReleaseDate", y = "AvgRating", queries = T), 
                                                             list(plot = scatter_plot,  x = "AvgRating", y = "Watches", queries = F)), ncols = 2), 
          query.legend = "bottom")
p4
dev.off()

##mantel
library(linkET)
library(ggplot2)
library(dplyr)
library(vegan)
library(openxlsx)
library(corrplot)
library(RColorBrewer)
bacteria_data= read.csv("./otu2025.csv", header=T, row.names=1, sep=",", comment.char="", stringsAsFactors = F)
env = read.csv("./2025EVN.csv", header=T, row.names=1, sep=",", comment.char="", stringsAsFactors = F)
otu_beta<-as.data.frame(t(bacteria_data))
mantel = mantel_test(
  spec = otu_beta, env = env,
  spec_dist =  dist_func(.FUN = "vegdist", method = "bray"), 
  env_dist = dist_func(.FUN = "vegdist", method = "euclidean"),
  mantel_fun = 'mantel',
  na_omit=TRUE,
)
cor_plot<-cor(varechem)
corrplot(cor_plot,method="color",type="upper",addrect=1,insig="blank",rect.col = "blue",rect.lwd = 2)
mt<-function(varespec,varechem)
  library(vegan)
library(dplyr)
vars <- colnames(varechem)
models<-list()
dev.new(
  
  title = "mantel test",
  
  width = 16,
  
  height = 8,
  
  noRStudioGD = TRUE
  
)
corrplot(cor(varechem),method = "color",type = "upper", addrect = 1,insig="blank",rect.col = "blue",rect.lwd = 2)

segments(df_segment$subx, df_segment$suby, df_segment$x, df_segment$y, lty = 'solid', lwd = df_segment$lwd, col = df_segment$lcol, xpd = TRUE)
points(subx, suby, pch = 24, col = 'blue', bg = 'blue', cex = 2, xpd = TRUE)
text(subx - 0.5, suby, labels = grp, adj = c(0.8, 0.5), cex = 1.2, xpd = TRUE)

#SEM
library("sem")
library("dplyr")
library("lavaan")
library(semPlot)
library(lavaan)
library(simsem)
library(devtools)
SEM <- read.csv("sem.csv", header=T, sep = ",")
SEM2 <- select(SEM,-c(Relationship2,Relationship,name1,name2))
fit <- lm(formula = dis~longitude+latitude+hei+prec+temp+host,data=SEM2)
summary(fit)
emPaths(fit,what = "paths",whatLabels = "std",nCharNodes = 0,style = "lisrel",
        layout = "tree2",rotation = 2,edge.color = "blue",color = "lightblue")
sem.model <- "dis~longitude+latitude+hei+prec+temp+host"
semmodel.fit <- sem(model=sem.model,data = SEM2)
summary(semmodel.fit, standardized = TRUE) 
cfi <- cfa(sem.model,data = SEM2,std.lv=T)
summary(semmodel.fit,standardized=T,rsq=T,fit.measures=T)
summary(semmodel.fit, std.nox = TRUE, rsquare = TRUE)
options(max.print = 10000)
result<-parameterestimates(semmodel.fit,boot.ci.type = "bca.simple",standardized = T)
semPaths(semmodel.fit,what = "paths",whatLabels = "est",nCharNodes = 0,label.cex=1.5,
         layout = "circle2",rotation = 1,style = "lisrel",
         edge.label.cex=1.2,esize=2,asize=2.5,label.prop=0.8,edge.label.bg=T,
         edge.color = "black")

#CCA
RDApath = paste(res1path,"/CCAenv/",sep = "")
dir.create(RDApath)
RDA_CCA = function(otu = NULL,tax = NULL,map = NULL,ps = NULL,env = env,group = "Group",path = "./",chose.env = F){
  dir.create(path)
  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  ps
  ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  mapping = as.data.frame( phyloseq::sample_data(ps_rela))
  env.dat = env
  samp.fg = colnames(otu)
  env.st = vegan::decostand(env.dat, method="standardize", MARGIN=2,na.rm = TRUE)
  samp.env= rownames(env.st)
  my.env = match(samp.fg, samp.env)
  env.st3 = na.omit(env.st[my.env, ])  
  samp.env= rownames(env.st3)
  my.fg = match(samp.env, samp.fg)
  otu = otu[, my.fg]
  env.st2=env.st3
  otu = t(otu)
  rowSums(otu)
  DCA= vegan::decorana(otu)  
  xxa = as.data.frame(DCA$rproj)
  if(max(xxa$DCA1)<=4 & max(xxa$DCA1)>=3){twochiose = "T"}else{twochiose = "F"}
  if(max(xxa$DCA1) > 4 | twochiose == "T") {
    C.whole = vegan::cca(otu, env.st3) 
    C.whole
    inf_factor = vegan::vif.cca(C.whole)
    na_env = which(is.na(inf_factor))
    if(isTRUE(length(na_env) > "0") ){
      inf_factor = inf_factor[-na_env]
    }
    
    max_env = which(inf_factor == max(inf_factor,na.rm = TRUE))
    env.st4 = env.st3
    while ( inf_factor[max_env] > 20){
      env.st4 = env.st4[,-max_env]
      C.reduced = vegan::cca(otu, env.st4)
      inf_factor = vegan::vif.cca(C.reduced)
      max_env = which(inf_factor == max(inf_factor,na.rm = TRUE))
    }
    output2 = inf_factor ;output2
    env.st4
    ind.p = array(0,dim=c(1,ncol(env.st4)))
    ind.F = array(0,dim=c(1,ncol(env.st4)))
    for(j in 1:ncol(env.st4)){
      ind.cca = vegan::cca(otu, env.st4[,j]) 
      ind.sig = anova(ind.cca,step=1000)
      ind.p[1,j] = ind.sig$Pr[1]
      ind.F[1,j] = ind.sig$F[1]
    }
    
    colnames(ind.p) = colnames(env.st4)
    inf_Fp=rbind(output2,ind.F,ind.p)
    row.names(inf_Fp)=c("inf_factor","F","p")
    C.whole = vegan::cca(otu, env.st4) 
    x.sig = anova(C.whole)
    x.p = x.sig$Pr[1] ;x.p
    x.F = x.sig$F[1]  ;x.F
    F1 <- paste("anova F: ",round(x.F, 2), sep = "")
    pp1 = paste("p: ",round(x.p, 2), sep = "")
    title = paste(F1," ",pp1, sep = "")
    
    output1 = summary(C.whole)
    
    str(output1)
    a=output1$sites;a 
    b=output1$cont$importance;b 
    c=output1$biplot*5;c  
    filenamea = paste(path,"cca_inf_Fp.txt",sep = "")
    write.table(inf_Fp,file=filenamea,sep="\t",col.names=NA)
    ca1=round(b[2,1],2);ca1
    ca2=round(b[2,2],2);ca2
    
    aa = merge(a,mapping, by="row.names",all=F);aa
    c = as.data.frame(c)
    
    library("ggplot2")
    library(ggrepel)
     p=ggplot()+
      geom_point(data=aa,aes(x=CCA1,y=CCA2, fill = Group),pch = 21,colour = "black",size = 4)+
      geom_segment(data = c,aes(x = 0, y = 0, xend = CCA1, yend = CCA2), 
                   arrow = arrow(angle=22.5,length = unit(0.2,"cm")),linetype=1, size=0.6,colour = "black")+
      geom_text_repel(data = c,aes(x=CCA1,y=CCA2, label = row.names(c)))+
      
      labs(x=paste("CCA 1 (", ca1*100, "%)", sep=""),
           y=paste("CCA 2 (", ca2*100, "%)", sep=""),
           title=title)
    
    p
    p = p+theme_bw()+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      theme(
       panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title = element_text(vjust = -8.5,hjust = 0.1),
        axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
        axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
        axis.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 14),
        legend.text = element_text(size = 15,face = "bold")
      ) 
    p
    plotnamea = paste(path,"/CCA_env.pdf",sep = "")
    ggsave(plotnamea, p, width = 8, height = 6, device = cairo_pdf, family = "Song")
    plotnamea2 = paste(path,"/CCA_env.jpg",sep = "")
    ggsave(plotnamea2, p, width = 8, height = 6, device = cairo_pdf, family = "Song")
  } else if (max(xxa$DCA1) < 3| twochiose == "T" ){
    C.whole = vegan::rda(otu, env.st2)
    if (chose.env == T) {
       inf_factor = vegan::vif.cca(C.whole)
      inf_factor
      na_env = which(is.na(inf_factor))
      if(isTRUE(length(na_env) > "0") ){
        inf_factor = inf_factor[-na_env]
      }
      max_env = which(inf_factor == max(inf_factor,na.rm = TRUE))
      env.st4 = env.st3
      while ( inf_factor[max_env] > 20){
        env.st4 = env.st4[,-max_env]
        C.reduced = vegan::cca(otu, env.st4)
        inf_factor = vegan::vif.cca(C.reduced)
        max_env = which(inf_factor == max(inf_factor,na.rm = TRUE))
      }
      output2 = inf_factor ;output2
     C.whole = vegan::rda(otu, env.st4)
    } else if (chose.env == F) {
      C.whole = vegan::rda(otu, env.st2)
      C.whole
      inf_factor = vegan::vif.cca(C.whole)
      inf_factor
      output2 = inf_factor ;output2
      ind.p = array(0,dim=c(1,ncol(env.st2)))
      ind.F = array(0,dim=c(1,ncol(env.st2)))
      for(j in 1:ncol(env.st2)){
        ind.cca = vegan::cca(otu, env.st2[,j]) 
        ind.sig = anova(ind.cca,step=1000)
        ind.p[1,j] = ind.sig$Pr[1]
        ind.F[1,j] = ind.sig$F[1]
      }
      colnames(ind.p) = colnames(env.st2)
      inf_Fp=rbind(output2,ind.F,ind.p)
      row.names(inf_Fp)=c("inf_factor","F","p")
    }
  x.sig = anova(C.whole)
    x.p = x.sig$Pr[1] ;x.p
    x.F = x.sig$F[1]  ;x.F
    F1 <- paste("anova F: ",round(x.F, 2), sep = "")
    pp1 = paste("p: ",round(x.p, 2), sep = "")
    title = paste(F1," ",pp1, sep = "")
    output1 = summary(C.whole)
    str(output1)
    a=output1$sites;a  
    b=output1$cont$importance;b 
    c=output1$biplot;c  
    ca1=round(b[2,1],2);ca1
    ca2=round(b[2,2],2);ca2
    exp = c(ca1,ca2)
    aa = merge(a,mapping, by="row.names",all=F);aa
    c = as.data.frame(c)
    p=ggplot()+
      geom_point(data=aa,aes(x=RDA1,y=RDA2, fill = Group),pch = 21,colour = "black",size = 4)+
      geom_segment(data = c,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
                   arrow = arrow(angle=22.5,length = unit(0.2,"cm")),linetype=1, size=0.6,colour = "black")+
      ggrepel::geom_text_repel(data = c,aes(x=RDA1,y=RDA2, label = row.names(c)))+
      labs(x=paste("RDA 1 (", ca1*100, "%)", sep=""),
           y=paste("RDA 2 (", ca2*100, "%)", sep=""),
           title=title)
    p
    p = p+theme_bw()+
      geom_hline(aes(yintercept=0), colour="black", linetype=2) +
      geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
      theme(
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.title = element_text(vjust = -8.5,hjust = 0.1),
        axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
        axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
        axis.text = element_text(size = 20,face = "bold"),
        axis.text.x = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 14),
        legend.text = element_text(size = 15,face = "bold")
      ) 
    p
  }
  p3 = p+ ggrepel::geom_text_repel(data=aa, aes_string(x=colnames(aa)[2],y=colnames(aa)[3],label=colnames(aa)[1]),size=4)#?stat_el
  p3
  return(list(p,aa,p3,inf_Fp,c,title,exp))
}
CCA_explain_percent = function(ps = ps,env.dat = envRDA){
  ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  mapping = as.data.frame( phyloseq::sample_data(ps_rela))
  samp.fg = colnames(otu)
  env.st = vegan::decostand(env.dat, method="standardize", MARGIN=2,na.rm = TRUE)
  samp.env= rownames(env.st)
  my.env = match(samp.fg, samp.env)
  env.st2 = na.omit(env.st[my.env, ])  
  samp.env= rownames(env.st2)
  my.fg = match(samp.env, samp.fg)
  otu = otu[, my.fg]
  otu = t(otu)
  C.whole = vegan::rda(otu, env.st2)
  C.whole
  inf_factor = vegan::vif.cca(C.whole)
  inf_factor
  na_env = which(is.na(inf_factor))
  if(isTRUE(length(na_env) > "0") ){
    inf_factor = inf_factor[-na_env]
  }
  max_env = which(inf_factor == max(inf_factor,na.rm = T))
  env.st4 = env.st2
  while ( inf_factor[max_env] > 20){
    env.st4 = env.st4[,-max_env]
    C.reduced = vegan::cca(otu, env.st4)
    inf_factor = vegan::vif.cca(C.reduced)
    max_env = which(inf_factor == max(inf_factor,na.rm = T ))
  }
  output2 = inf_factor ;output2
  C.whole = vegan::rda(otu, env.st4)  
  total.chi = C.whole$tot.chi
  ind.p = array(0,dim=c(1,ncol(env.st4)))
  for(j in 1:ncol(env.st4)){
    ind.par = vegan::rda(otu, env.st4[,j], env.st4[,-j])
    ind.chi = ind.par$CCA$tot.chi
    ind.per = ind.chi/total.chi
    ind.p[j] = ind.per
  }
  ind.p
  rowname = colnames(env.st4);rowname
  out = matrix(data=NA,ncol=length(colnames(env.st4)),nrow=1);out
  out = ind.p
  rownames(out) = "percent"
  colnames(out) = rowname
  out
  total.chi = C.whole$tot.chi;total.chi
  total.constrained = C.whole$CCA$tot.chi ; total.constrained
  explained.percent = (total.constrained) / total.chi;explained.percent
  unexplained.percent = (total.chi - total.constrained) / total.chi;unexplained.percent
  exp = data.frame(ID = c("explained.percent","unexplained.percent"),count = c(explained.percent,unexplained.percent))
  exp
  return(list(out,exp))
}
env = read.csv("./env.csv")
result = RDA_CCA(ps = ps_Rlefse,
                 env = envRDA,
                 path = RDApath,
                 chose.env = F
)
p1 = result[[1]] + mytheme1 + scale_fill_manual(values = colset1)
p1
dataplot = result[[2]]
p2 = result[[3]]+ mytheme1  + scale_fill_manual(values = colset1)
aov = result[[4]]
plotnamea = paste(RDApath,"/CCA_envnew.pdf",sep = "")
ggsave(plotnamea, p1, width = 8, height = 6)
plotnamea4 = paste(RDApath,"/CCA_envnew.jpg",sep = "")
ggsave(plotnamea4, p1, width = 8, height = 6)
filenamea = paste(RDApath,"dataplotnew.txt",sep = "")
write.table(dataplot ,file=filenamea,sep="\t",col.names=NA)
filenamea = paste(RDApath,"aovnew.txt",sep = "")
write.table(aov,file=filenamea,sep="\t",col.names=NA)
plotnamea = paste(RDApath,"/RDA_envlabelnew.pdf",sep = "")
ggsave(plotnamea, p2, width = 18, height = 12)
plotnamea4 = paste(RDApath,"/RDA_envlabelnew.png",sep = "")
ggsave(plotnamea4, p2, width = 18, height = 12)
result = RDA_CCA_explain_percent(ps = ps_Rlefse,
                                 env.dat = envRDA)
out = result[[1]]
wxp = result[[2]]
filenamea = paste(RDApath,"each_env_exp_percent.csv",sep = "")
write.csv(out,filenamea)
filenamea = paste(RDApath,"all_index_explain_percent.csv",sep = "")
write.csv(exp,filenamea)










