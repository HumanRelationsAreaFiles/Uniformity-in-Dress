library(rethinking)
library(tidyverse)
library(phytools)
library(geosphere)
library(igraph)
library(qgraph)
library(corrplot)
library(corpcor)
library(ggridges)
library(plgp)
library(viridis)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(cubelyr)

# Read in study data
d <- read.csv("data.csv")

# Bring in phylogenetic supertree
tree <- ape::read.tree("SCCS_supertree.tre")

setdiff(tree$tip.label, d$socname) # checking for discrepancies between phylo tree names and dataframe names, should return "character(0)" if all is well

d$id <- 1:nrow(d) # Give a new integer identifier to each society

##############################################################
### Extract clothing and adornment codes
d_CA <- d[,c(2:12, 25:28, ncol(d))]

# We need to recode so that the lowest level of each code = 1
d_CA_long <- d_CA %>%
  # Convert to long form
  pivot_longer(-c(socname,id)) %>%
  # Create new index within each outcome
  group_by(name) %>%
  arrange(value, .by_group=T) %>%
  mutate(y = match(value, unique(value[!is.na(value)])))


# Give each outcome an integer code
d_CA_long$resp <- match(d_CA_long$name, unique(d_CA_long$name))

# Get number of unique ordinal categories for each response
K_resp <- d_CA_long %>% 
  group_by(resp) %>%
  summarise(K_resp = max(y, na.rm=T), name=unique(name))

# Reader-friendly labels
K_resp$label <- c(bquote(AO), bquote(A[P]*S[M]), bquote(A[P]*S[F]), bquote(A[P]*R[M]), bquote(A[P]*R[F]), bquote(A[NP]*S[M]), bquote(A[NP]*S[F]), bquote(A[NP]*R[M]), bquote(A[NP]*R[F]), bquote(CS[M]), bquote(CS[F]), bquote(CR[M]), bquote(CR[F]), bquote(CO))

#################################
### Extract TL codes ############
d_TL <- select(d, c(General_TL5_Final, Soc_TL5_Final, Gender_TL5_Final, Mar_TL5_Final, Sex_TL5_Final, FM_TL5_Final, socname, id))

d_TL_long <- d_TL %>%
  pivot_longer(-c(socname,id)) %>%
  group_by(name) %>%
  mutate(TL = round(value), TL_diff = value - round(value))

# give an integer identifier for each TL measure
d_TL_long$TL_j <- match(d_TL_long$name, unique(d_TL_long$name))

#################################
# Index NAs as -99 for Stan model
d_CA_long[is.na(d_CA_long)] <- -99
d[is.na(d)] <- -99

## Create phylogenetic distance matrix
dist_phy <- cophenetic.phylo(tree) # pairwise distance matrix using branch lengths
dist_phy <- dist_phy / max(dist_phy) # scaling matrix to [0,1]

## Create temporal distance matrix, adding very small jitter on date
set.seed(5)
dist_EP <- as.matrix( dist(jitter(d$EP), method='euclidean'), nrow=nrow(d), ncol=nrow(d) )
# Scaling time marix
dist_EP <- dist_EP / max(dist_EP)

## Create geographic distance matrix, little jitter again to keep matrix positive definite
long_lat <- c("Revised.longitude", "Revised.latitude")
dist_geo <- distm( jitter(as.matrix(d[,long_lat])) ) # pairwise geographic distance matrix; great circle distance
dist_geo <- dist_geo / max(dist_geo)
colnames(dist_geo) <- colnames(dist_phy)
rownames(dist_geo) <- colnames(dist_geo)

#############################################################
##### Part 1: Assess C&A dimensionality #####################
#############################################################

### Organize data into a list
data_list <- list(
  N = max(d_CA_long$id),
  N_obs = nrow(d_CA_long),
  J = max(d_CA_long$resp),
  y = d_CA_long$y,
  K = K_resp$K_resp,
  resp = d_CA_long$resp,
  id = d_CA_long$id,
  cor_phy = 1 - dist_phy,
  cor_EP = 1 - dist_EP,
  cor_geo = 1- dist_geo
)

# Fit model with RStan, uncomment to run for first time
# m_EGA <- stan_model("stan_code/m_EGA.stan")
#fit_EGA <- sampling( m_EGA, iter=1000, chains=8, cores=8, data=data_list, control=list(adapt_delta=0.9), init="0" )
#save(fit_EGA, file="fit_EGA.RData")
load("fit_EGA.RData")

post <- extract.samples(fit_EGA) # extract posterior samples

##### Plot correlation and regularized partial correlation network, along with EGA clustering
cor_med <- cor( apply(post$eta_v, 2:3, median) )
cor_med <- round(cor_med, digits=3) # greater consistency across MCMC runs
rownames(cor_med) <- K_resp$label
colnames(cor_med) <- rownames(cor_med)

## export corr mat
cor_df <- as.data.frame(cor_med)
cor_df$var <- names(cor_df)
cor_df <- cor_df[, c(names(cor_df)[length(cor_df)], names(cor_df)[-length(cor_df)]) ]

write_csv(cor_df, "fig_2_left_correlations.csv")

## Sparse GGM with GLASSO algorithm
ggm <- as.matrix(EBICglasso(cor_med,n=data_list$N),gamma=0.5, threshold=F, refit=TRUE)
ggm.g <- graph_from_adjacency_matrix(abs(ggm), weighted = T, mode="undirected", diag=F)

## Community detection
cw <- cluster_walktrap(ggm.g)
cl <- cluster_louvain(ggm.g) # multilevel louvain has higher modularity and better performance in simulation studies than walktrap, therefore use this community detection algorithm

K_resp$factor <- c(
  2, # A0
  1, # APSM
  1, # APSF
  4, # APRM
  1, # APRF
  2, # ANPSM
  2, # ANPSF
  3, # ANPRM
  3, # ANPRF
  3, # CSM
  4, # CSF
  4, # CRM
  4, # CRF
  4 # C0
)

#### Plot correlation and partial correlation networks ##
svg( "cor_pcor.svg", width=12, height=6, pointsize=12 )

par(mfrow=c(1,2))

g_cor <- qgraph(cor_med,
       graph="cor",
       layout="spring",
       posCol = "#70a494",
       negCol = "#de8a5a",
       labels = K_resp$label,
       DoNotPlot=T,
       borders = T,
       fade = T ,
       title = ""
)

plot(g_cor)

lwd <- round(quantile(g_cor$graphAttributes$Edges$width)[2:5],1)
col <- g_cor$graphAttributes$Edges$color[match( lwd, round(g_cor$graphAttributes$Edges$width,1))]

max_cor <- cor_med; diag(max_cor) <- 0; max_cor <- max(max_cor)

legend("topleft",legend=round(c(lwd/lwd[4])*max_cor,2),col=col,lwd=lwd,bty="n")


g_pcor <- qgraph(cor_med,
       graph="EBICglasso",
       sampleSize=data_list$N,
       color = c("#b14241", "#2D9486", "#E0972D", "skyblue"),
       posCol = "black", # plotting absolute value of pcor to avoid misleading interpretations from spurious edges!
       negCol = "black",
       labels = K_resp$label,
       groups = as.character(K_resp$factor),
       borders = T,
       tuning = 0.5,
       DoNotPlot=T,
       threshold = F,
       legend = F,
       fade = T ,
       title = "")

plot(g_pcor)

lwd <- round(quantile(g_pcor$graphAttributes$Edges$width, probs=seq(from=0,to=1, length.out = 6))[3:6],1 )
col <- g_pcor$graphAttributes$Edges$color[match( lwd, round(g_pcor$graphAttributes$Edges$width,1))]

max_cor <- ggm; diag(max_cor) <- 0; max_cor <- max(max_cor)

legend("topleft",legend=round(c(lwd/lwd[4])*max_cor,2),col=col,lwd=lwd,bty="n")

dev.off()

## Export partial correlations
ggm <- round( as.matrix(EBICglasso(cor_med,n=data_list$N),gamma=0.5, threshold=F, refit=TRUE), 3 )

pcor_df <- as.data.frame(ggm)
pcor_df$var <- names(pcor_df)
pcor_df <- pcor_df[, c(names(pcor_df)[length(pcor_df)], names(pcor_df)[-length(pcor_df)]) ]

write_csv(pcor_df, "fig_2_right_partial_correlations.csv")
