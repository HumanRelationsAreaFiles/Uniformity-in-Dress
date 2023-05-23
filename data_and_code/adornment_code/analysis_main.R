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
### Extract RS codes ############
d_RS <- select(d, c(Chronic_Scarcity, Natural_Hazards, Famine, socname, id))

d_RS_long <- d_RS %>%
  pivot_longer(-c(socname,id))

# give an integer identifier for each RS measure
d_RS_long$RS_j <- match(d_RS_long$name, unique(d_RS_long$name))

#################################
### Recode political codes #####
records <- ifelse(d$v149 == 5, 1, 0)
strat <- ifelse(d$v158 > 1, 1, 0)

## quick checks for multicolinearity
library(car)
vif( lm(d$C9 ~ records + strat + d$v157 + d$v235) )
vif( lm(d$A11 ~ records + strat + d$v157 + d$v235) )
vif( lm(d$C2a ~ records + strat + d$v157 + d$v235) )
vif( lm(d$A2a ~ records + strat + d$v157 + d$v235) )

#################################
# Index NAs as -99 for Stan model
d_CA_long[is.na(d_CA_long)] <- -99
d_TL_long[is.na(d_TL_long)] <- -99
d_RS_long[is.na(d_RS_long)] <- -99
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

###### Factor model ############################
# Indicate which factor each response belongs to
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

# Indicate which factor loadings to fix (1 per cluster)
K_resp$fix <- ifelse(K_resp$label %in% c("CO", "A[P] * S[M]", "AO", "A[NP] * R[M]"), 1, 0)

#############################################################
##### Part 2: Assess associations with C&A ##################
#############################################################

### Organize data into a list
data_list <- list(
  N = max(d_CA_long$id),
  N_obs = nrow(d_CA_long),
  J = max(d_CA_long$resp),
  N_f = 4,
  y = d_CA_long$y,
  K = K_resp$K_resp,
  factor = K_resp$factor,
  fix = K_resp$fix,
  resp = d_CA_long$resp,
  id = d_CA_long$id,
  
  y_RS = d_RS_long$value,
  id_RS = d_RS_long$id,
  RS_j = d_RS_long$RS_j,
  N_RS = nrow(d_RS_long),
  
  y_TL = d_TL_long$TL,
  TL_diff = d_TL_long$TL_diff,
  id_TL = d_TL_long$id,
  TL_j = d_TL_long$TL_j,
  N_TL = nrow(d_TL_long),
  
  strat = strat,
  records = records,
  comm_size = d$v235 - 1, # subtracting 1 so that 0 is the minimum
  pol_int = d$v157 - 1, # same
  
  cor_phy = 1 - dist_phy,
  cor_EP = 1 - dist_EP,
  cor_geo = 1- dist_geo
)

### Uncomment and run each of these models and save the results 

##### Assess correlations with Tight-Loose ############
#m_TL <- stan_model("stan_code/factor_GLMM_TL.stan")
#fit_TL <- sampling( m_TL, iter=1000, chains=8, cores=8, data=data_list, control=list(adapt_delta=0.9), init="0" )
#save(fit_TL, file = "fit_TL.RData")

##### Assess correlations with resource stress ############
#m_RS <- stan_model("stan_code/factor_GLMM_RS.stan")
#fit_RS <- sampling( m_RS, iter=1000, chains=8, cores=8, data=data_list, control=list(adapt_delta=0.9), init="0" )
#save(fit_RS, file = "fit_RS.RData")

##### Assess correlations with strat ############
#m_strat <- stan_model("stan_code/factor_GLMM_strat.stan")
#fit_strat <- sampling( m_strat, iter=1000, chains=8, cores=8, data=data_list, control=list(adapt_delta=0.9), init="0" )
#save(fit_strat, file = "fit_strat.RData")

##### Assess correlations with pol int ############
#m_pli <- stan_model("stan_code/factor_GLMM_pli.stan")
#fit_pli <- sampling( m_pli, iter=1000, chains=8, cores=8, data=data_list, control=list(adapt_delta=0.9), init="0" )
#save(fit_pli, file = "fit_pli.RData")

##### Assess correlations with records ############
#m_records <- stan_model("stan_code/factor_GLMM_records.stan")
#fit_records <- sampling( m_records, iter=1000, chains=8, cores=8, data=data_list, control=list(adapt_delta=0.9), init="0" )
#save(fit_records, file = "fit_records.RData")

##### Assess correlations with comm size ############
#m_comm <- stan_model("stan_code/factor_GLMM_comm.stan")
#fit_comm <- sampling( m_comm, iter=1000, chains=8, cores=8, data=data_list, control=list(adapt_delta=0.9), init="0" )
#save(fit_comm, file = "fit_comm.RData")

##### Multivariate model: all predictors ############
#m_mv <- stan_model("stan_code/factor_GLMM_mv.stan")
#fit_mv <- sampling( m_mv, iter=1000, chains=8, cores=8, data=data_list, control=list(adapt_delta=0.9), init="0" )
#save(fit_mv, file = "fit_mv.RData")
#####################################################

### What are the correlations between C&A and Tight-Loose?
load("fit_TL.RData")
post <- extract.samples( fit_TL )
n_samps <- length(post$lp__)

eta_labels <- c(bquote(eta[1]), bquote(eta[2]), bquote(eta[3]), bquote(eta[4]), bquote(eta[TL]))

Rho_eta <- apply(post$Rho_eta, 2:3, median)

svg("Rho_eta.svg", width = 6, height=6)

par(cex=1.5)

qgraph(Rho_eta,
       graph="cor",
       layout="circle",
       posCol = "#70a494",
       negCol = "#de8a5a",
       color = c("#b14241", "#2D9486", "#E0972D", "skyblue", "slategray"),
       labels = eta_labels,
       edge.label.color = "#70a494",
       edge.labels=T,
       borders = T,
       fade = T ,
       title = "")

dev.off()

Rho_eta_long <-  as.data.frame(post$Rho_eta[,1:4,5]) 
colnames(Rho_eta_long) <- c(expression(eta[1]), expression(eta[2]), expression(eta[3]), expression(eta[4]))
Rho_eta_long$samp <- 1:n_samps

Rho_eta_summary <- Rho_eta_long %>% pivot_longer(-samp) %>% 
  group_by(name) %>% 
  summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value, prob=0.9)[2])

svg("eta_TL_cor.svg", width = 4, height=6.5)

ggplot(Rho_eta_summary, aes(x=med, y=fct_rev(name))) + 
  geom_point(size=4) +
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0, lwd=2) +
  scale_x_continuous(expand = c(0.00, 0), limits = c(-1,1))  + 
  scale_y_discrete(labels = rev(c(expression(eta[1]), expression(eta[2]), expression(eta[3]), expression(eta[4])))) +
  geom_vline(xintercept = 0, lty="dashed") +
  xlab("Correlation with Tight-Loose") +
  ylab("") +
  theme_bw(base_size=16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(), legend.position = "none", axis.text.y = element_text(colour=rev(c("#b14241", "#2D9486", "#E0972D", "skyblue")), size=18))

dev.off()

### What is the variance captured for each latent variable by phy, geo, EP, res
S_median <- as.data.frame(t(apply(post$S, 2:3, median)))
colnames(S_median) <- c(expression(eta[1]), expression(eta[2]), expression(eta[3]), expression(eta[4]), expression(eta[TL]))
S_median$source <- c("Phylogenetic", "Geographic", "Ethnographic\n Present", "Residual")

S_median_long <- S_median %>% pivot_longer(-source)

svg("S_median.svg", width = 6, height=3)

ggplot(S_median_long, aes(x=name, y=value, fill=fct_rev(source))) +
  geom_bar(position="fill", stat="identity", alpha=0.9) +
  scale_fill_viridis(discrete = T, option = "E") +
  scale_x_discrete(labels = c(expression(eta[1]), expression(eta[2]), expression(eta[3]), expression(eta[4]), expression(eta[TL]))) +
  scale_y_continuous(labels=scales::percent) +
  xlab("") +
  ylab("% Variance")  +
  theme_bw(base_size=15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(), legend.position = "top",axis.text.x = element_text(colour=(c("#b14241", "#2D9486", "#E0972D", "skyblue", "slategray")), size=18))

dev.off()

#######################################################################################

fac_names <- c(expression(paste(eta[1], "(z-score)")), expression(paste(eta[2], "(z-score)")), expression(paste(eta[3], "(z-score)")), expression(paste(eta[4], "(z-score)")), expression(paste(eta[TL], "(z-score)")))
fac_cols <- c("#b14241", "#2D9486", "#E0972D", "skyblue", "slategray")

load("fit_mv.RData")

#### Correlations between resource stress and eta ############################
# Univariate #####
load("fit_RS.RData")
post <- extract.samples( fit_RS )
n_samps <- length(post$lp__)

sigma_RS <- post$sigma_eta_RS
b_RS <- post$b_RS

for (f in 1:5) b_RS[,f] <- (b_RS[,f] * sigma_RS*2) / (apply(post$eta_v[,,f], 1, sd))

b_RS_long <- as.data.frame(b_RS) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_RS_uv_summary <- b_RS_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Multivariate #####
post <- extract.samples( fit_mv )
n_samps <- length(post$lp__)

sigma_RS <- post$sigma_eta_RS
b_RS <- post$b_RS

for (f in 1:5) b_RS[,f] <- (b_RS[,f] * sigma_RS*2) / (apply(post$eta_v[,,f], 1, sd))

b_RS_long <- as.data.frame(b_RS) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_RS_mv_summary <- b_RS_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Put univariate and multivaraite estimates together
b_RS_all <- bind_rows(b_RS_uv_summary, b_RS_mv_summary)
b_RS_all$model <- c( rep("Univariate", nrow(b_RS_uv_summary)), rep("Multivariate", nrow(b_RS_mv_summary)))

#### Correlations between strat and eta ############################
# Univariate #####
load("fit_strat.RData")
post <- extract.samples( fit_strat )
n_samps <- length(post$lp__)

b_strat <- post$b_strat

for (f in 1:5) b_strat[,f] <- (b_strat[,f]) / (apply(post$eta_v[,,f], 1, sd))

b_strat_long <- as.data.frame(b_strat) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_strat_uv_summary <- b_strat_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Multivariate #####
post <- extract.samples( fit_mv )
n_samps <- length(post$lp__)

b_strat <- post$b_strat

for (f in 1:5) b_strat[,f] <- (b_strat[,f]) / (apply(post$eta_v[,,f], 1, sd))

b_strat_long <- as.data.frame(b_strat) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_strat_mv_summary <- b_strat_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Put univariate and multivaraite estimates together
b_strat_all <- bind_rows(b_strat_uv_summary, b_strat_mv_summary)
b_strat_all$model <- c( rep("Univariate", nrow(b_strat_uv_summary)), rep("Multivariate", nrow(b_strat_mv_summary)))
################################################################

#### Correlations between records and eta ############################
# Univariate #####
load("fit_records.RData")
post <- extract.samples( fit_records )
n_samps <- length(post$lp__)

b_records <- post$b_records

for (f in 1:5) b_records[,f] <- (b_records[,f]) / (apply(post$eta_v[,,f], 1, sd))

b_records_long <- as.data.frame(b_records) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_records_uv_summary <- b_records_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Multivariate #####
post <- extract.samples( fit_mv )
n_samps <- length(post$lp__)

b_records <- post$b_records

for (f in 1:5) b_records[,f] <- (b_records[,f]) / (apply(post$eta_v[,,f], 1, sd))

b_records_long <- as.data.frame(b_records) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_records_mv_summary <- b_records_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Put univariate and multivaraite estimates together
b_records_all <- bind_rows(b_records_uv_summary, b_records_mv_summary)
b_records_all$model <- c( rep("Univariate", nrow(b_records_uv_summary)), rep("Multivariate", nrow(b_records_mv_summary)))
################################################################

#### Correlations between pol int and eta ############################
# Univariate #####
load("fit_pli.RData")
post <- extract.samples( fit_pli )
n_samps <- length(post$lp__)

b_int <- post$b_int

for (f in 1:5) b_int[,f] <- (b_int[,f]) / (apply(post$eta_v[,,f], 1, sd))

b_pli_long <- as.data.frame(b_int) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_pli_uv_summary <- b_pli_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Multivariate #####
post <- extract.samples( fit_mv )
n_samps <- length(post$lp__)

b_int <- post$b_int

for (f in 1:5) b_int[,f] <- (b_int[,f]) / (apply(post$eta_v[,,f], 1, sd))

b_pli_long <- as.data.frame(b_int) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_pli_mv_summary <- b_pli_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Put univariate and multivaraite estimates together
b_pli_all <- bind_rows(b_pli_uv_summary, b_pli_mv_summary)
b_pli_all$model <- c( rep("Univariate", nrow(b_pli_uv_summary)), rep("Multivariate", nrow(b_pli_mv_summary)))
################################################################

#### Correlations between pol int and eta ############################
# Univariate #####
load("fit_comm.RData")
post <- extract.samples( fit_comm )
n_samps <- length(post$lp__)

b_comm <- post$b_comm

for (f in 1:5) b_comm[,f] <- (b_comm[,f]) / (apply(post$eta_v[,,f], 1, sd))

b_comm_long <- as.data.frame(b_comm) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_comm_uv_summary <- b_comm_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Multivariate #####
post <- extract.samples( fit_mv )
n_samps <- length(post$lp__)

b_comm <- post$b_comm

for (f in 1:5) b_comm[,f] <- (b_comm[,f]) / (apply(post$eta_v[,,f], 1, sd))

b_comm_long <- as.data.frame(b_comm) %>% mutate(samp = 1:n_samps) %>% pivot_longer(-samp)
b_comm_mv_summary <- b_comm_long %>% group_by(name) %>% summarise(med = median(value), lower=HPDI(value, prob=0.9)[1], upper=HPDI(value,prob=0.9)[2])

# Put univariate and multivaraite estimates together
b_comm_all <- bind_rows(b_comm_uv_summary, b_comm_mv_summary)
b_comm_all$model <- c( rep("Univariate", nrow(b_comm_uv_summary)), rep("Multivariate", nrow(b_comm_mv_summary)))
################################################################
#### Put ALL of the univariate and multivariate estimates together

b_est_all <- bind_rows(b_RS_all, b_strat_all, b_records_all, b_pli_all, b_comm_all)
b_est_all$resp <- rep( c("Resource Stress", "Social Stratification", "Written Records", "Political Integration", "Mean Community Size"), each = 10)


png("Fig4.png", width=8, height=6, units='in', res=600)

ggplot(b_est_all, aes(x=med, y=fct_rev(name), color=model)) + 
  facet_wrap(~fct_relevel(resp, "Resource Stress", "Social Stratification", "Political Integration", "Written Records","Mean Community Size"), scales="free_x") +
  geom_point(position = position_dodge(0.5),size=2) +
  scale_y_discrete(labels = rev(c(expression(eta[1]), expression(eta[2]), expression(eta[3]), expression(eta[4]), expression(eta[TL])))) +
  geom_vline(xintercept=0, lty="dashed") +
  geom_errorbarh(aes(xmin=lower,xmax=upper), height=0, position = position_dodge(0.5),lwd=1.2) + 
  scale_color_manual(values=c("coral","black")) +
  theme_bw(base_size=15) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),strip.background = element_blank(),legend.title = element_blank(), axis.text.y = element_text(colour=rev(c("#b14241", "#2D9486", "#E0972D", "skyblue", "slategray")), size=18), panel.spacing = unit(1.5, "lines"), legend.position = c(0.78,0.35)) +
  xlab("Standardized Mean Difference (+2 SD)") + 
    ylab("")

dev.off()


b_est_all$name <- fct_recode(b_est_all$name, eta1 = "V1", eta2 = "V2", eta3 = "V3", eta4 = "V4", eta_TL = "V5")


write_csv(b_est_all, "SMD_estimates.csv")

######################################################################
#### Exploratory plots: more standardization in "mid-range" societies?
load("fit_TL.RData")

post <- extract.samples(fit_TL)

eta <- as.data.frame(apply(post$eta_v, 2:3, median))

# Standardize variables
for (j in 1:ncol(eta))
  eta[,j] = as.numeric(scale(eta[,j]))

names(eta) <- c("eta1","eta2","eta3","eta4","etaTL")
eta$id <- d$id

eta_long <- eta %>% pivot_longer(-id)
eta_long <- left_join(eta_long, select(d, c(id, v157, v158, v235, v151)))


eta_long$name <- factor(eta_long$name, labels = c(
  "eta1" = expression(eta[1]),
  "eta2" = expression(eta[2]),
  "eta3" = expression(eta[3]),
  "eta4" = expression(eta[4]),
  "etaTL" = expression(eta[TL])
)
)

#### Political Integration Figure
png("Fig_polint.png", width=9.5, height=6, units='in', res=600)

ggplot(eta_long, aes(x=as.factor(v157), y=value, color=as.factor(v157), fill=as.factor(v157))) + 
  facet_wrap(~name, scales="free",labeller =label_parsed) +
  geom_boxplot(alpha=0.2) +
  geom_jitter(alpha=0.5,width=0.05) +
  scale_color_viridis(discrete=TRUE, option="D", labels=c("None", "Autonomous local communities", "1 level above community", "2 levels above community", "3 levels above community")) +
  scale_fill_viridis(discrete=TRUE, option="D",labels=c("None", "Autonomous local communities", "1 level above community", "2 levels above community", "3 levels above community")) +
  scale_x_discrete(breaks=c(1,5), labels=c("Less \nIntegration", "More \nIntegration")) +
  ggtitle("") +
  ylab("z-score") +
  xlab("") +
  theme_minimal(base_size=12) +
  theme(legend.position = c(0.82,0.23),panel.spacing = unit(2, "lines"), legend.title = element_blank(),strip.text.x = element_text(size = 15)) +
  ggtitle("Political Integration")

dev.off()


#### Strat Figure
png("Fig_strat.png", width=9.5, height=6, units='in', res=600)

ggplot(eta_long, aes(x=as.factor(v158), y=value, color=as.factor(v158), fill=as.factor(v158))) + 
  facet_wrap(~name, scales="free",labeller =label_parsed) +
  geom_boxplot(alpha=0.2) +
  geom_jitter(alpha=0.5,width=0.05) +
  scale_color_viridis(discrete=TRUE, option="D", labels=c("Egalitarian", "Status/Wealth Inequality", "2 classes/castes, no slavery", "2 classes/castes, with slavery", "3+ classes/castes")) +
  scale_fill_viridis(discrete=TRUE, option="D",labels=c("Egalitarian", "Status/Wealth Inequality", "2 classes/castes, no slavery", "2 classes/castes, with slavery", "3+ classes/castes")) +
  scale_x_discrete(breaks=c(1,5), labels=c("Less \nStratification", "More \nStratification")) +
  ggtitle("") +
  ylab("z-score") +
  xlab("") +
  theme_minimal(base_size=12) +
  theme(legend.position = c(0.82,0.23),panel.spacing = unit(2, "lines"), legend.title = element_blank(),strip.text.x = element_text(size = 15)) +
  ggtitle("Social Stratification")

dev.off()


#### Agriculture Figure
png("Fig_agr.png", width=9.5, height=6, units='in', res=600)

ggplot(eta_long, aes(x=as.factor(v151), y=value, color=as.factor(v151), fill=as.factor(v151))) + 
  facet_wrap(~name, scales="free",labeller =label_parsed) +
  geom_boxplot(alpha=0.2) +
  geom_jitter(alpha=0.5,width=0.05) +
  scale_color_viridis(discrete=TRUE, option="D", labels=c("None","<10% of food supply", ">10% of food supply", "Primary, not intensive", "Primary, intensive")) +
  scale_fill_viridis(discrete=TRUE, option="D",labels=c("None","<10% of food supply", ">10% of food supply", "Primary, not intensive", "Primary, intensive")) +
  scale_x_discrete(breaks=c(1,5), labels=c("Less \nIntensification", "More \nIntensification")) +
  ggtitle("") +
  ylab("z-score") +
  xlab("") +
  theme_minimal(base_size=12) +
  theme(legend.position = c(0.82,0.23),panel.spacing = unit(2, "lines"), legend.title = element_blank(),strip.text.x = element_text(size = 15)) +
  ggtitle("Agriculture")

dev.off()


#### Community size Figure
png("Fig_comm.png", width=9.5, height=6, units='in', res=600)

ggplot(filter(eta_long, v235 > 0), aes(x=as.factor(v235), y=value, color=as.factor(v235), fill=as.factor(v235))) + 
  facet_wrap(~name, scales="free",labeller =label_parsed) +
  geom_boxplot(alpha=0.2) +
  geom_jitter(alpha=0.5,width=0.05) +
  scale_color_viridis(discrete=TRUE, option="D", labels=c("<50","50-99", "100-199", "200-399", "400-1000", ">1000", "One or more towns of 5000-50,0000", "One or more cities of 50,000+")) +
  scale_fill_viridis(discrete=TRUE, option="D",labels=c("<50","50-99", "100-199", "200-399", "400-1000", ">1000", "One or more towns of 5000-50,0000", "One or more cities of 50,000+")) +
  scale_x_discrete(breaks=c(1,7), labels=c("Smaller \nCommunities", "Larger \nCommunities")) +
  ggtitle("") +
  ylab("z-score") +
  xlab("") +
  theme_minimal(base_size=12) +
  theme(legend.position = c(0.82,0.23),panel.spacing = unit(2, "lines"), legend.title = element_blank(),strip.text.x = element_text(size = 15)) +
  ggtitle("Mean size of local community")

dev.off()

###############################################################
#### Correlation between RS latent and observed ###############
load("fit_RS.RData")

post <- extract.samples(fit_RS)

cor_chrs <- sqrt( apply(post$eta_RS_v * 1, 1, var) / (apply(post$eta_RS_v * 1, 1, var) + pi^2/3)  )
cor_haz <- sqrt( apply(post$eta_RS_v * post$lambda_RS[,1], 1, var) / (apply(post$eta_RS_v * post$lambda_RS[,1], 1, var) + pi^2/3)  ) * sign(post$lambda_RS[,1])
cor_fam <- sqrt( apply(post$eta_RS_v * post$lambda_RS[,2], 1, var) / (apply(post$eta_RS_v * post$lambda_RS[,2], 1, var) + pi^2/3)  ) * sign(post$lambda_RS[,2])

cor_RS <- cbind(cor_chrs, cor_haz, cor_fam)
round(apply(cor_RS, 2, median), 4)

round(apply(cor_RS, 2, HPDI, prob=0.9), 4)

#### Correlation between TL latent and observed ###############
load("fit_TL.RData")

post <- extract.samples(fit_TL)

cor_general <- sqrt( apply(post$eta_v[,,5] * 1, 1, var) / (apply(post$eta_v[,,5] * 1, 1, var) + pi^2/3)  )

cor_soc <- sqrt( apply(post$eta_v[,,5] * post$lambda_TL[,1], 1, var) / (apply(post$eta_v[,,5] * post$lambda_TL[,1], 1, var) + pi^2/3)  ) * sign(post$lambda_TL[,1])

cor_gender <- sqrt( apply(post$eta_v[,,5] * post$lambda_TL[,2], 1, var) / (apply(post$eta_v[,,5] * post$lambda_TL[,2], 1, var) + pi^2/3)  ) * sign(post$lambda_TL[,2])

cor_marr <- sqrt( apply(post$eta_v[,,5] * post$lambda_TL[,3], 1, var) / (apply(post$eta_v[,,5] * post$lambda_TL[,3], 1, var) + pi^2/3)  ) * sign(post$lambda_TL[,3])

cor_sex <- sqrt( apply(post$eta_v[,,5] * post$lambda_TL[,4], 1, var) / (apply(post$eta_v[,,5] * post$lambda_TL[,4], 1, var) + pi^2/3)  ) * sign(post$lambda_TL[,4])

cor_FM <- sqrt( apply(post$eta_v[,,5] * post$lambda_TL[,5], 1, var) / (apply(post$eta_v[,,5] * post$lambda_TL[,5], 1, var) + pi^2/3)  ) * sign(post$lambda_TL[,5])


cor_TL <- cbind(cor_general, cor_soc, cor_gender, cor_marr, cor_sex, cor_FM)
round(apply(cor_TL, 2, median), 4)

round(apply(cor_TL, 2, HPDI, prob=0.9), 4)

############################################################
#### Changes in variables across ethnographic present ######
load("fit_TL.RData")
post <- extract.samples(fit_TL)

n_samps <- length(post$lp__)

post_EP <- array( post$EP_v, dim = dim(post$EP_v), dimnames = list(samp = 1:n_samps, obs = 1:data_list$N, var = c(expression(eta[1]), expression(eta[2]), expression(eta[3]), expression(eta[4]), expression(eta[TL])) ) )

post_eta_sd <- as.data.frame(apply(post$eta_v, c(1,3), sd))
names(post_eta_sd) <- c(expression(eta[1]), expression(eta[2]), expression(eta[3]), expression(eta[4]), expression(eta[TL]))
post_eta_sd$samp <- 1:n_samps

eta_sd_long <- pivot_longer(post_eta_sd, -samp, values_to = "sd")

EP_long <- post_EP %>% 
  cubelyr::as.tbl_cube(met_name = "est") %>% 
  as_tibble %>% 
  left_join(eta_sd_long) %>% 
  mutate(est_z = est/sd)

d_EP_obs <- data.frame(
  obs = 1:data_list$N,
  EP = d$EP)

EP_summary <- EP_long %>% 
  group_by(obs, var) %>% 
  summarise(
    med = median(est_z),
    lower = HPDI(est_z, prob=0.9)[1],
    upper = HPDI(est_z, prob=0.9)[2]
  ) %>% 
  left_join(d_EP_obs)

head(EP_summary)

png("Fig_EP.png", width=9.5, height=6, units='in', res=600)

ggplot(filter(EP_summary, EP > 1800), aes(x = EP, y = med)) +
  facet_wrap(~var,labeller =label_parsed) +
  geom_point(alpha=0.5) +
  geom_errorbar(aes(x = EP, ymin=lower, ymax=upper), alpha=0.5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  ylab("EP parameter z-score") +
  theme_minimal(base_size=12)

dev.off()


















