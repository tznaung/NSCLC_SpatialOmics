# install.packages("anndata")
# reticulate::install_miniconda()
# anndata::install_anndata()
library(anndata)
library(umap)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(RColorBrewer)

ad <- read_h5ad("4301_ytma_nsclc_CTfinal.h5ad")
expr = ad$X
clinInfo = ad$obs

sel = which(clinInfo$TMA %in% c('YTMA_471'))
pats = unique(clinInfo[sel,]$Image)
nPats = length(pats)
tm_idxs = which(clinInfo$CT_final=="Tumour")
m2_idxs = which(clinInfo$CT_final=="M2_mac")

# test for correlation
tb = matrix(0,nPats,2)
for(i in 1:nPats){
  print(i)
  idxs = which(clinInfo$Image==pats[i])
  idxs_tm = intersect(idxs,tm_idxs)
  idxs_m2 = intersect(idxs,m2_idxs)
  
  #tb[i,1] = mean(expr[idxs_tm,21])   # pdl1 tm / prop m2      # 0.06585484, p-value = 0.5045 (YTMA_471)
  #tb[i,2] = length(idxs_m2) / length(idxs)
  
  #tb[i,1] = mean(expr[idxs_tm,21]) # pdl1 tm / pdl1 m2 :  0.6298228, p-value < 2.2e-16 (all)
  #tb[i,2] = mean(expr[idxs_m2,21])   
  
  #tb[i,1] = mean(expr[idxs,21])   # pdl1 all-cts / prop m2       # 0.2578501, p-value = 0.007332 (YTMA_471)
  #tb[i,2] = length(idxs_m2) / length(idxs)
  
  #tb[i,1] = mean(expr[idxs,21])   # pdl1 all-cts / prop tm       # -0.269804 , p-value = 0.004948 (YTMA_471)
  #tb[i,2] = length(idxs_tm) / length(idxs)  
  
  #tb[i,1] = sum(expr[idxs,43]) # pck / pdl1 :  0.1337014 , p-value = 0.1698 (YTMA_471)
  #tb[i,2] = sum(expr[idxs,21])  
  
  tb[i,1] = mean(expr[idxs,25]) # cd163 / pdl1 :  0.3625109  , p-value = 0.0001245 (YTMA_471)
  tb[i,2] = mean(expr[idxs,21])    
}
rownames(tb) = pats
colnames(tb) = c('tm_pld1','m2_prop')
tb
res = cor.test(tb[,1],tb[,2])
res

# rank celltypes by PDL1 expression
cellTypes  = unique(clinInfo$CT_final)
nCT = length(cellTypes)
PDL1_expr = rep(0,nCT)
PDL1_expr_tot = rep(0,nCT)
PDL1_idx = which(colnames(expr)=="PD-L1")
PD1_idx = which(colnames(expr)=="PD-1")
GzmB_idx = which(colnames(expr)=="Granzyme B")
ICOS_idx = which(colnames(expr)=="ICOS")
Ki67_idx = which(colnames(expr)=="Ki67")
CD163_idx = which(colnames(expr)=="CD163")
LAG3_idx = which(colnames(expr)=="LAG3")

for (i in 1:nCT){
  idx = which(clinInfo$CT_final==cellTypes[i])
  PDL1_expr[i] = mean(expr[idx,PDL1_idx])
  PDL1_expr_tot[i] = sum(expr[idx,PDL1_idx])
}
names(PDL1_expr) = cellTypes
names(PDL1_expr_tot) = cellTypes

PDL1_expr = sort(PDL1_expr,decreasing=TRUE)
print(PDL1_expr)

PDL1_expr_tot = sort(PDL1_expr_tot,decreasing=TRUE)
print(PDL1_expr_tot)

### umap plot

plot.Umap <- 
  function(x, labels,
           main="A UMAP visualization of the NSCLC-ITx cohort",
           colors=c("coral2", "cyan"),
           pad=0.1, cex=0.65, pch=19, add=FALSE, legend.suffix="",
           cex.main=1, cex.legend=1) {
    
    layout = x
    if (is(x, "umap")) {
      layout = x$layout
    } 
    
    xylim = range(layout)
    xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
    if (!add) {
      par(mar=c(0.2,0.7,1.2,0.7), ps=10)
      plot(xylim, xylim, type="n", axes=F, frame=F)
      rect(xylim[1], xylim[1], xylim[2], xylim[2], border="white", lwd=0.25)  
    }
    points(layout[,1], layout[,2], col=colors[as.integer(labels)],
           cex=cex, pch=pch)
    mtext(side=3, main, cex=cex.main)
    
    labels.u = unique(labels)
    legend.pos = "topright"
    legend.text = as.character(labels.u)
    if (add) {
      legend.pos = "bottomright"
      legend.text = paste(as.character(labels.u), legend.suffix)
    }
    legend(legend.pos, legend=legend.text,
           col=colors[as.integer(labels.u)],
           bty="n", pch=pch, cex=cex.legend)
  }

# select data

sel_idx = which(clinInfo$TMA=="YTMA_471")

data_only <- expr[sel_idx,]
data_only <- scale(data_only)
data_category <- clinInfo$Response[sel_idx]
data_cellType <- clinInfo$CT_final[sel_idx]

umap.defaults # The default configuration object is called umap.defaults. 
# This is a list encoding default values for all the parameters used within the algorithm.
custom.config = umap.defaults
custom.config$min_dist = 0.1 # change the default minimum distance, default is 0.1
custom.config$random_state = 123 # change the default random_state (seed)
custom.config$n_neighbors = 15
custom.config$n_components = 2

um.umap = umap(data_only, custom.config) # use umap to generate two dimentional data. 
# use custom.config to tell what parameters you are using.
um.umap

plot.Umap(um.umap,data_category)

color_group1 <- c("gold","gold1","gold2","gold3","gold4","goldenrod",
                  "goldenrod1","goldenrod2","goldenrod3",
                  "aquamarine","aquamarine1","aquamarine2","aquamarine3",
                  "aquamarine4","chartreuse","chartreuse1","chartreuse2",
                  "chartreuse3","chartreuse4","green","green1") 

color_group2 <- c("green2","green3","green4","greenyellow","darkgreen","cyan",
                  "orangered","orangered1","orangered2","orangered3",
                  "orangered4","red","red1","red2","red3","red4",
                  "firebrick","firebrick1","firebrick2","firebrick3",
                  "firebrick4","blue","blue1","blue2","blue3","blue4","darkblue",
                  "dodgerblue","dodgerblue1","dodgerblue2","dodgerblue3",
                  "dodgerblue4","slateblue","slateblue1","slateblue2",
                  "slateblue3","slateblue4","steelblue","steelblue4",
                  "chartreuse3","chartreuse4","green","green1")

color_group3 <- c("red","blue")

color_group4_bcells <- c("yellow","black","black","black","black","black","black","black","black",
                  "black","black","black","black","black","black","black","black","black",
                  "black","black","black","black","black","black","black","black","black")
color_group4_mphage <- c("black","black","black","black","black","black","black","black","black",
                  "red","yellow","black","black","green","black","black","black","black",
                  "black","black","black","black","black","black","black","black","black")

color_group5_allCellTypes <- c("black","dodgerblue","dodgerblue3","slateblue","slateblue3","cyan","slateblue4","#771155",
                  "green","red","yellow","#AAAA44","blue","aquamarine3","azure2","azure3")


tol21rainbow= c("#771122", "#AA4455", "#DD7788","#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77")

um.umap$layout
um.umap_data_label_merged  <- data.frame(uMap1 = um.umap$layout[, 1], uMap2 = um.umap$layout[, 2], 
                                         PDL1 = expr[sel_idx,PDL1_idx], PD1 = expr[sel_idx,PD1_idx],
                                         GzmB = expr[sel_idx,GzmB_idx], Ki67 = expr[sel_idx,Ki67_idx],
                                         ICOS = expr[sel_idx,ICOS_idx], CD163 = expr[sel_idx,CD163_idx],
                                         LAG3 = expr[sel_idx,LAG3_idx])
um.umap_data_label_merged$PDL1[um.umap_data_label_merged$PDL1>20] = 20 # removing outlier cells
um.umap_data_label_merged$PD1[um.umap_data_label_merged$PD1>5] = 10
um.umap_data_label_merged$GzmB[um.umap_data_label_merged$GzmB>5] = 10
um.umap_data_label_merged$Ki67[um.umap_data_label_merged$Ki67>20] = 20
um.umap_data_label_merged$ICOS[um.umap_data_label_merged$ICOS>5] = 10
um.umap_data_label_merged$CD163[um.umap_data_label_merged$CD163>100] = 100
um.umap_data_label_merged$LAG3[um.umap_data_label_merged$LAG3>2.5] = 10

um.umap_data_label_merged$uMap2[um.umap_data_label_merged$uMap2<=-12] = -12
um.umap_data_label_merged$uMap1[um.umap_data_label_merged$uMap1<=-15] = -15
um.umap_data_label_merged$uMap1[um.umap_data_label_merged$uMap1>=15] = 15
ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = as.factor(data_cellType))) +
  geom_point(size = 0.05) +
  theme_bw() +
  scale_color_manual(values = color_group5_allCellTypes) + theme_pubr()+
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 3))

ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = as.factor(data_category))) +
  geom_point(size = 0.05) +
  theme_bw() +
  scale_color_manual(values = color_group3) + theme_pubr()+
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 3))

# ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = PDL1)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("PD-L1 Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
# ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = PD1)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("PD-1 Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = PDL1)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("PD-L1 Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlGn")))(50))
ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = PD1)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("PD-1 Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlGn")))(50))
ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = GzmB)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("GzmB Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlGn")))(50))
ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = Ki67)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("Ki67 Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlGn")))(50))
ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = ICOS)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("ICOS Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlGn")))(50))
ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = CD163)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("CD163 Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlGn")))(50))
ggplot(um.umap_data_label_merged,  aes(x = uMap1, y = uMap2, color = LAG3)) + geom_point(size = 0.01) + theme_classic() + scale_color_gradientn("LAG3 Expression", colours = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlGn")))(50))


