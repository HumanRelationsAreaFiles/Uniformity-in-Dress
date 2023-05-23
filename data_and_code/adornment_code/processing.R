library(tidyverse)
library(readxl)
library(corrplot)
library(phytools)

d <- read_xlsx("data-raw/CA_June_2020.xlsx")

##### Data dictionary #####################

### Cluster 1: standardization and constraint
# RESTL-C2a - "To what extent is clothing standardized among adult males?
## 1 = not standardized
## 2 = slightly standardized
## 3 = moderately standardized
## 4 = very standardized
## 4.5 = inferred very standardized

# RESTL-C2b - "To what extent is clothing standardized among adult females?
## same scale as C2a

# C3 - Rules about clothing

# RESTL-C9 - "Overall, with specific regard to clothing, how free (vs. constrained) are people in this culture to act as they please?"
## 1 = Not at all constrained
## 2 = Somewhat constrained
## 3 = Moderately constrained
## 4 = Very constrained

# RESTL-A11 - "Overall, with regard to adornment, how constrained are individuals to adorn themselves a certain way?"
## 1 = Not at all constrained
## 2 = Somewhat constrained
## 3 = Moderately constrained
## 4 = Very constrained

# RESTL-A2a - "To what extent is permanent adornment standardized among adult males?"
## 1 = not standardized
## 2 = slightly standardized
## 3 = moderately standardized
## 4 = very standardized
## 4.5 = inferred very standardized

# RESTL-A2b - "To what extent is permanent adornment standardized among adult females?"
## same scale as A2a

# RESTL-A3 - Permament adornment rules

# RESTL-A6a - "To what extent is nonpermanent adornment standardized among adult males?"
## 1 = not standardized
## 2 = slightly standardized
## 3 = moderately standardized
## 4 = very standardized
## 4.5 = inferred very standardized

# RESTL-A7 - nonpermanent adornment rules

###################################

### Cluster 2: Rules ##############
# RESTL-C3a - "To what extent do rules dictate adult male clothing?”
## Dichtomize into "absent" (0) or "present" (>0)
# RESTL-C3b - "To what extent do rules dictate adult female clothing?”
## Same as C3a
###################################

### Cluster 3: Gender #############
# RESTL-C1 - "To what extent (excluding special purpose clothing, and other than accommodations for body build differences) do males and females dress differently?"
## 0 = little to no difference, 3 = great differentiation
# RESTL-A1c - "If there are permanent markers, to what extent do males and females have different permanent markers?”
## same scale as C1
###################################

#### Now, extract relevant variables
CA_vars <- c("sccsn", "RESTL-C2a", "RESTL-C2b", "RESTL-C9", "RESTL-A2a", "RESTL-A2b", "RESTL-A6a", "RESTL-A6b", "RESTL-A11", "RESTL-C3a", "RESTL-C3b", "RESTL-C1", "RESTL-A1c")

other_vars <- c("sccsn","socname","EP","Famine","Natural_Hazards","Chronic_Scarcity")

# Subset  CA variables
d_CA <- d[,CA_vars] %>%
  # Change codes to numeric
  mutate_if( is.character, as.numeric ) %>%
  # "floor" rounds down codes, so that we get rid of the 0.5s in an appropriate way
  mutate_all( floor )

# Clean up the variable names, getting rid of "RESTL-"
names(d_CA)[-1] <- substr( names(d_CA)[-1], 7, nchar(names(d_CA)[-1]) ) 
d_CA <- as.data.frame(d_CA)

# Calculate how many NAs per society
num_NA <- apply(d_CA, 1, function(x) sum(is.na(x)))

# Indicate rows where ALL of the CA codes are NA, because those societies not relevant to this analysis
d_CA$no_CA <- ifelse( num_NA < (ncol(d_CA)-1), 0, 1)

# Subset other variables
d2 <- d[,other_vars]
d2[,4:6] <- d2[,4:6] %>%
  mutate_if( is.character, as.numeric ) %>%
  mutate_all( ceiling ) # resource stress codes are rounded up in the case of rare .5, as in previous publication on subdiv
 
# Bring together
d_both <- left_join(d_CA, d2)

## Bring in some new vars
d_SCCS <- read_csv("data-raw/SCCS-var1-2000.csv")
d_SCCS <- d_SCCS[,c("socname","v149", "v157", "v158", "v235", "v151")]

d_all <- left_join(d_both, d_SCCS)
d_all$strat <- ifelse(d_all$v158 > 1, 1, 0)

# Drop CA variables we won't use in this analysis
d_all <- d_all %>% select(-c(C1, A1c))
##################################################
##### Bring in phylogeny #########################
sccs_tree <- read.nexus("data-raw/Time-calibrated SCCS supertree.nex")
setdiff(sccs_tree$tip.label, d_all$socname) # checking for discrepancies between phylo tree names and dataframe names

# "MISSING" indicates that genetic data missing for the target pop, position based on linguistic data
d_all[c(1,2,14,25,27,29,30,33,39,46,48,53,59,63,66,82,93,97,102,106,124,133,134,135,136,137,138,139,140,143,144,146,158,164,167,174,176,178), "socname"] <-  c("Nama_Hottentot", "Kung_Bushmen", "Nkundo_Mongo", "Pastoral_Fulani", "Massa_Masa", "Fur_Darfur", "Otoro_Nuba", "MISSING_Kafa_Kaffa", "Kenuzi_Nubians","Rwala_Bedouin", "Gheg_Albanians", "Yurak_Samoyed", "Punjabi_West", "Uttar_Pradesh", "Khalka_Mongols", "Negri_Sembilan", "MISSING_Kimam", "New_Ireland", "Mbau_Fijians", "Western_Samoans", "Copper_Eskimo", "MISSING_Twana", "MISSING_Yurok", "MISSING_Pomo_Eastern", "Yokuts_Lake", "Paiute_North", "MISSING_Klamath","MISSING_Kutenai", "Gros_Ventre", "MISSING_Omaha", "MISSING_Huron", "MISSING_Natchez", "Cuna_Tule", "Carib_Barama", "Cubeo_Tucano", "MISSING_Nambicuara", "MISSING_Timbira", "MISSING_Botocudo")

setdiff(sccs_tree$tip.label, d_all$socname) # should be 0 now

# Drop non CA cases, won't be useful for this analysis
d_all <- d_all[d_all$no_CA == 0,]

# Prune tree to cases relevant for this analysis
sccs_tree <- drop.tip(sccs_tree, subset(sccs_tree$tip.label, !(sccs_tree$tip.label %in% d_all$socname)))

# Putting the phylogeny and the data frame in the same order, this is very important!
d_all <- d_all[match(sccs_tree$tip.label, d_all$socname),]

# Now we can remove the "MISSING" from the labels
d_all$socname <- ifelse(substr(d_all$socname, 1, nchar("MISSING")) == "MISSING", substr(d_all$socname, nchar("MISSING")+2, nchar(d_all$socname)), d_all$socname)

# Do the same for phylogeny
sccs_tree$tip.label <- ifelse(substr(sccs_tree$tip.label, 1, nchar("MISSING")) == "MISSING", substr(sccs_tree$tip.label, nchar("MISSING")+2, nchar(sccs_tree$tip.label)), sccs_tree$tip.label)

# Check one more time!
setdiff(sccs_tree$tip.label, d_all$socname) # should be 0
#########################################################
##### Bring in location data from DPLACE ################
d_loc <- read.csv("data-raw/location.csv", stringsAsFactors = F)

d_loc <- d_loc %>%
  filter(Source == "Standard cross-cultural sample") %>%
  select(Society.id, Revised.latitude, Revised.longitude) %>% 
  mutate(sccsn = substr(Society.id, 5, nchar(Society.id))) %>%
  mutate(sccsn = as.numeric(sccsn)) %>%
  select(-Society.id)

d_all <- left_join(d_all, d_loc)

# Manually adding lat/lon for one society where data missing from DPLACE
d_all$Revised.latitude[d_all$socname == "Iban"] <- 2
d_all$Revised.longitude[d_all$socname == "Iban"] <- 110

########################################################
##### Adding in extra CA codes and TL codes ############
d_CA_2 <- read_xlsx("data-raw/CA rules2.xlsx")
d_TL <- read_xlsx("data-raw/TL Codes.xlsx")

names(d_CA_2)[3:ncol(d_CA_2)] <- substr( names(d_CA_2)[3:ncol(d_CA_2)], 7, nchar(names(d_CA_2)[3:ncol(d_CA_2)]))

d_all <- left_join(d_all, select(d_CA_2, -c(C3a, C3b, socname)))
d_all <- left_join(d_all, select(d_TL, c(sccsn, General_TL5_Final, Soc_TL5_Final, Gender_TL5_Final, Mar_TL5_Final, Sex_TL5_Final, FM_TL5_Final) ))

# Export processed files for analysis
write.csv(d_all %>% select(-no_CA), "data.csv", row.names = F)
write.tree(sccs_tree, "SCCS_supertree.tre")
