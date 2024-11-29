library(tidyverse)
library(readr)
library(mclust)
library(clusterCrit)
library(readxl)

library(dplyr)
library(tidyr)
library(cluster)
library(caret)

library(TraMineR)
library(seqhandbook)
library(dtw)
library(cluster)
library(datasets)

library(ComplexHeatmap)
library(cowplot)
library(latex2exp)


# Set work space
setwd("F:/phd/1-Longitudinal/Multiview/Experimental_journal")

#Input clustering algorithm
file.sources <- list.files(path = "Utils",pattern="*.R$", full.names=TRUE,ignore.case=TRUE, recursive = TRUE)
file.sources <- file.sources[!grepl("evaluation.R", file.sources, ignore.case = TRUE)]
sapply(file.sources,source,.GlobalEnv)

#Input Dataset
source('Code/datasets.R')


# Data
data(iris)
sepal.v <- iris %>% select(Sepal.Length, Sepal.Width)
petal.v <- iris %>% select(Petal.Length, Petal.Width)
species <- iris %>% select(Species) 
y <- species$Species %>% as.factor() %>% as.numeric()


sepal.v_ <- as.matrix(dist(sepal.v, method = "euclidean"))
petal.v_ <- as.matrix(dist(petal.v, method = "euclidean"))
x <- list(sepal.v_, petal.v_)


# clustering
# mecmdd_lp <- mecmdd.rwl(Xlist=x, c=3, type='full',
#                         alpha=2,
#                         beta=1.5,
#                         delta=10, epsi=1e-3, disp=TRUE,
#                         gamma=0.5, eta=1.5, weight='prod')
# mecmdd_gp <- mecmdd.rwg(Xlist=x, c=3, type='full',
#                         alpha=2,
#                         beta=1.5,
#                         delta=10, epsi=1e-3, disp=TRUE,
#                         gamma=0.5, eta=1.5, weight='prod')
# 
# mecmdd_ls <- mecmdd.rwl(Xlist=x, c=3, type='full',
#                         alpha=2,
#                         beta=1.5,
#                         delta=10, epsi=1e-3, disp=TRUE,
#                         gamma=0, eta=1, weight='sum')
# mecmdd_gs <- mecmdd.rwg(Xlist=x, c=3, type='full',
#                         alpha=2,
#                         beta=1.5,
#                         delta=10, epsi=1e-3, disp=TRUE,
#                         gamma=0, eta=1, weight='sum')
# 
# 
# df_l <- data.frame(ls=c(mecmdd_ls$param$lambda[,1],mecmdd_ls$param$lambda[,2]),
#                    lp=c(mecmdd_lp$param$lambda[,1],mecmdd_lp$param$lambda[,2]))
# df_g <- data.frame(gs=c(mecmdd_gs$param$lambda),
#                    gp=c(mecmdd_gp$param$lambda))
# write.csv(df_l, "Iris/df_l_weight.csv")
# write.csv(df_g, "Iris/df_g_weight.csv")


# save(mecmdd_ls, mecmdd_lp, mecmdd_gs, mecmdd_gp,
#      file = paste0("Iris/iris_model_weigt.Rdata"))







load("Iris/iris_model.Rdata")

# performance
Ptrue <- pairwise_mass(create_hard_credpart(y))
ri <- matrix(0, nrow = 4, ncol = 1, dimnames = list(c("ls", "lp", "gs", "gp")))
colnames(ri) <- 'RI'
ri['ls',] <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(mecmdd_ls),type="c")
ri['lp',] <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(mecmdd_lp),type="c")
ri['gs',] <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(mecmdd_gs),type="c")
ri['gp',] <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(mecmdd_gp),type="c")


ns <- matrix(c(mecmdd_ls$N,mecmdd_lp$N,mecmdd_gs$N,mecmdd_gp$N), 
             nrow = 4, ncol = 1, dimnames = list(c("ls", "lp", "gs", "gp")))
colnames(ns) <- 'nonspecificity'



asw <- matrix(0, nrow = 4, ncol = 1, dimnames = list(c("ls", "lp", "gs", "gp")))
colnames(asw) <- 'asw'
labels <- as.integer(as.vector(apply(mecmdd_ls$mass,1,which.max)))
labels <- true_lab(labels,y)
asw['ls',] <- weighted_silhouette(labels, x, mecmdd_ls$param$lambda, local=TRUE)$asw

labels <- as.integer(as.vector(apply(mecmdd_lp$mass,1,which.max)))
labels <- true_lab(labels,y)
asw['lp',] <- weighted_silhouette(labels, x, mecmdd_lp$param$lambda, local=TRUE)$asw

labels <- as.integer(as.vector(apply(mecmdd_gs$mass,1,which.max)))
labels <- true_lab(labels,y)
asw['gs',] <- weighted_silhouette(labels, x, mecmdd_gs$param$lambda, local=FALSE)$asw

labels <- as.integer(as.vector(apply(mecmdd_gp$mass,1,which.max)))
labels <- true_lab(labels,y)
asw['gp',] <- weighted_silhouette(labels, x, mecmdd_gp$param$lambda, local=FALSE)$asw




### ARI
load("Iris/iris_model.Rdata")

y_ <- ifelse(y==1,5,y)

ls <- as.integer(as.vector(apply(mecmdd_ls$mass,1,which.max)))
gs <- as.integer(as.vector(apply(mecmdd_gs$mass,1,which.max)))
gs <- ifelse(gs==3,5, ifelse(gs==5,3,gs))
lp <- as.integer(as.vector(apply(mecmdd_lp$mass,1,which.max)))
gp  <- as.integer(as.vector(apply(mecmdd_gp$mass,1,which.max)))
gp <- ifelse(gp==3,5, ifelse(gp==5,3,gp))

####
tables <- list(table(y_, ls), table(y_, gs), table(y_, lp), table(y_, gp))

# nam <- c(TeX('$\\omega_1$'), TeX('$\\omega_2$'), TeX('$\\omega_1,\\omega_2$'),
#          TeX('$\\omega_3$'), TeX('$\\omega_1,\\omega_3$'),
#          TeX('$\\omega_2,\\omega_3$'), TeX('$\\Omega'))
nam <- c('ω₁', 'ω₂', 'ω₁,ω₂', 'ω₃', 'ω₁,ω₃', 'ω₂,ω₃', 'Ω')
for (i in seq_along(tables)) {
  colnames(tables[[i]]) <- nam[as.numeric(colnames(tables[[i]]))-1]
  rownames(tables[[i]]) <- nam[as.numeric(rownames(tables[[i]]))-1]
}
table1 <- tables[[1]]
table2 <- tables[[2]]
table3 <- tables[[3]]
table4 <- tables[[4]]

col_fun <- circlize::colorRamp2(c(0, 1, 10, 20,50),
                                c("white","#00feef", "#19ceeb","white", "white"))

make_cell_fun <- function(mat) {
  function(j, i, x, y, width, height, fill) {
    grid::grid.text(sprintf("%.0f", mat[i, j]), x, y, gp = gpar(fontsize = 12,
                                                                fontfamily = "serif", fontface = "bold"))
    grid::grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "grey", fill = NA))
  }}


p1 <- Heatmap(table1, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
              show_column_names = TRUE, column_names_rot = 45, show_heatmap_legend=FALSE,
              cell_fun = make_cell_fun(table1), row_names_side = "left",
              row_title = "True Clusters", column_title="MECMdd-RWL-S",
              column_title_side ="bottom", column_names_side="top",
              column_names_gp = gpar(fontsize = 16, fontfamily = "serif", fontface = "bold", col = "black"),
              row_names_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "black"),
              row_title_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "grey40"),
              column_title_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "grey40"),
              border_gp = gpar(col = "white"))

p2 <- Heatmap(table2, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 45,
              cell_fun = make_cell_fun(table2), show_heatmap_legend=FALSE,
              column_title="MECMdd-RWG-S", column_title_side ="bottom",
              column_names_side="top", border_gp = gpar(col = "white"),
              column_names_gp = gpar(fontsize = 16, fontfamily = "serif", fontface = "bold", col = "black"),
              row_names_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "grey40"),
              column_title_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "grey40"))

p3 <- Heatmap(table3, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 45,
              cell_fun = make_cell_fun(table3), show_heatmap_legend=FALSE,
              column_title="MECMdd-RWL-P", column_title_side ="bottom",
              column_names_side="top", border_gp = gpar(col = "white"),
              column_names_gp = gpar(fontsize = 16, fontfamily = "serif", fontface = "bold", col = "black"),
              row_names_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "grey40"),
              column_title_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "grey40"))

p4 <- Heatmap(table4, col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = FALSE, show_column_names = TRUE, column_names_rot = 45,
              cell_fun = make_cell_fun(table4), column_title="MECMdd-RWG-P",
              column_title_side ="bottom", column_names_side="top", show_heatmap_legend=T,
              column_names_gp = gpar(fontsize = 16, fontfamily = "serif", fontface = "bold", col = "black"),
              row_names_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "grey40"),
              column_title_gp = gpar(fontsize = 14, fontfamily = "serif", fontface = "bold", col = "grey40"),
              heatmap_legend_param = list(title = "Confusion",
                                          labels_gp = gpar(fontsize = 10, fontfamily = "serif", fontface = "bold", col = "grey40"),
                                          title_gp = gpar(fontsize = 10, fontfamily = "serif", fontface = "bold", col = "grey40")),
              border_gp = gpar(col = "white"))

List <-  HeatmapList(p1 + p2 + p3 + p4)
ComplexHeatmap::draw(List, gap = unit(8, "mm"))




# -------------------------------------------------------------------
# Analyse parametres
beta_init <- c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0) 
delta_init <- c(10) #,'Q'
gamma_init <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
eta_init <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
alpha_init <- 2
Ptrue <- evclust::pairwise_mass(evclust::create_hard_credpart(y))
param <- data.frame(algorithme=as.character(), beta = as.numeric(), 
                    gamma = as.numeric(), eta=as.numeric(), cri = as.numeric())


for(b in beta_init){
  for(g  in gamma_init){
    for(e in eta_init){
      tryCatch({
        ls <- mecmdd.rwl(Xlist=x, c=3, type='simple',
                         alpha=2,
                         beta=b,
                         delta=5, epsi=1e-3, disp=FALSE,
                         gamma=g, eta=e, weight='sum')
        ls_ri <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(ls),type="c")
        param_ <- c('ls',b, g, e,ls_ri)
      }, error = function(e) {
        param_ <- c('ls',b, g, e, NA)
      })
      param <- rbind.data.frame(param, param_)
      
      tryCatch({
        lp <- mecmdd.rwl(Xlist=x, c=3, type='simple',
                         alpha=2,
                         beta=b,
                         delta=5, epsi=1e-3, disp=FALSE,
                         gamma=g, eta=e, weight='prod')
        lp_ri <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(lp),type="c")
        param_ <- c('lp',b, g, e,lp_ri)
      }, error = function(e) {
        param_ <- c('lp',b, g, e, NA)
      })
      param <- rbind.data.frame(param, param_)
      
      tryCatch({
        gs <- mecmdd.rwg(Xlist=x, c=3, type='simple',
                         alpha=2,
                         beta=b,
                         delta=5, epsi=1e-3, disp=FALSE,
                         gamma=g, eta=e, weight='sum')
        gs_ri <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(gs),type="c")
        param_ <- c('gs',b, g, e,gs_ri)
      }, error = function(e) {
        param_ <- c('gs',b, g, e, NA)
      })
      param <- rbind.data.frame(param, param_)
      
      tryCatch({
        gp <- mecmdd.rwg(Xlist=x, c=3, type='simple',
                         alpha=2,
                         beta=b,
                         delta=5, epsi=1e-3, disp=FALSE,
                         gamma=g, eta=e, weight='prod')
        gp_ri <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(gp),type="c")
        param_ <- c('gp',b, g, e,gs_ri)
      }, error = function(e) {
        param_ <- c('gp',b, g, e, NA)
      })
      param <- rbind.data.frame(param, param_)

      cat('beta=',b,' gamma=',g, ' eta=', e, '\n')
    }
  }
}
colnames(param) <- c('algorithme', 'beta', 'gamma', 'eta', 'cri')
write.csv(param, "Iris/cri_parameters.csv")





# -------------------------------------------------------------------
# Delta Parameters
delta_init <- c(0.9,'Q')

Ptrue <- evclust::pairwise_mass(evclust::create_hard_credpart(y))
param <- data.frame(iter =as.numeric(), algorithme=as.character(), delta = as.character(), 
                    cri = as.numeric(), asw = as.numeric())

iter <- 10
for(i in 1:iter){
  for(d in delta_init){
    if(d=='Q'){
      del <- NULL
    } else{
      del <- as.numeric(d)
    }
    tryCatch({
      ls <- mecmdd.rwl(Xlist=x, c=3, type='full',
                       alpha=2,
                       beta=1.5,
                       delta=del, epsi=1e-3, disp=FALSE,
                       gamma=0, eta=1, weight='sum')
      ls_ri <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(ls),type="c")
      labels <- as.integer(as.vector(apply(ls$mass,1,which.max)))
      labels <- true_lab(labels,y)
      ls_asw <- weighted_silhouette(labels, x, ls$param$lambda, local=TRUE)$asw
      param_ <- c(i,'ls',d,ls_ri,ls_asw)
    }, error = function(e) {
      param_ <- c(i,'ls',d, NA, NA)
    })
    param <- rbind(param, param_)
    
    tryCatch({
      lp <- mecmdd.rwl(Xlist=x, c=3, type='full',
                       alpha=2,
                       beta=1.5,
                       delta=del, epsi=1e-3, disp=FALSE,
                       gamma=0, eta=1, weight='prod')
      lp_ri <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(lp),type="c")
      labels <- as.integer(as.vector(apply(lp$mass,1,which.max)))
      labels <- true_lab(labels,y)
      lp_asw <- weighted_silhouette(labels, x, lp$param$lambda, local=TRUE)$asw
      param_ <- c(i,'lp',d,lp_ri,lp_asw)
    }, error = function(e) {
      param_ <- c(i,'lp',d, NA, NA)
    })
    param <- rbind(param, param_)
    
    tryCatch({
      gs <- mecmdd.rwg(Xlist=x, c=3, type='full',
                       alpha=2,
                       beta=1.5,
                       delta=del, epsi=1e-3, disp=FALSE,
                       gamma=0, eta=1, weight='sum')
      gs_ri <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(gs),type="c")
      labels <- as.integer(as.vector(apply(gs$mass,1,which.max)))
      labels <- true_lab(labels,y)
      gs_asw <- weighted_silhouette(labels, x, gs$param$lambda, local=FALSE)$asw
      param_ <- c(i,'gs',d,gs_ri,gs_asw)
    }, error = function(e) {
      param_ <- c(i,'gs',d, NA,NA)
    })
    param <- rbind(param, param_)
    
    tryCatch({
      gp <- mecmdd.rwg(Xlist=x, c=3, type='full',
                       alpha=2,
                       beta=1.5,
                       delta=del, epsi=1e-3, disp=FALSE,
                       gamma=0, eta=1, weight='prod')
      gp_ri <- evclust::credal_RI(Ptrue,evclust::pairwise_mass(gp),type="c")
      labels <- as.integer(as.vector(apply(gp$mass,1,which.max)))
      labels <- true_lab(labels,y)
      gp_asw <- weighted_silhouette(labels, x, gp$param$lambda, local=FALSE)$asw
      param_ <- c(i,'gp',d,gp_ri,gp_asw)
    }, error = function(e) {
      param_ <- c(i,'gp',d, NA,NA)
    })
    param <- rbind(param, param_)
    
    cat(i,' delta=',del, '\n')
  }
}
colnames(param) <- c('iter', 'algorithme', 'delta', 'cri', 'asw' )
write.csv(param, "Iris/cri_parameters_deltat_2.csv")
param_<- read_csv("Iris/cri_parameters_deltat_2.csv")
paramM <- aggregate(cbind(cri, asw) ~ algorithme + delta, data = param_, 
                    FUN = function(x) round(mean(x, na.rm = FALSE), 3))



## Courbe de convergence
ls <- mecmdd.rwl(Xlist=x, c=3, type='full',
                       alpha=2,
                       beta=1.5,
                       delta=10, epsi=1e-4, disp=TRUE,
                       gamma=0.5, eta=1, weight='sum')

lp <- mecmdd.rwl(Xlist=x, c=3, type='full',
                        alpha=2,
                        beta=1.5,
                        delta=10, epsi=1e-4, disp=TRUE,
                        gamma=0.5, eta=1, weight='prod')

gs <- mecmdd.rwg(Xlist=x, c=3, type='full',
                        alpha=2,
                        beta=1.5,
                        delta=10, epsi=1e-4, disp=TRUE,
                        gamma=0.5, eta=1, weight='sum')

gp <- mecmdd.rwg(Xlist=x, c=3, type='full',
                        alpha=2,
                        beta=1.5,
                        delta=10, epsi=1e-4, disp=TRUE,
                        gamma=0.5, eta=1, weight='prod')
