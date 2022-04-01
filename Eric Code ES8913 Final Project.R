

# Packages ----------------------------------------------------------------

library(readr)
library(data.table)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(tidyr)

# Data Manipulations ------------------------------------------------------

#data<-fread(paste0(getwd(),"/contaminant_database.csv"))
data <- fread("C:/Users/ericf/Desktop/Ryerson Masters/ES8913/contaminant_database.csv")

##data manipulations
data<-data %>% rowwise() %>%
  mutate(Total_Ph =
           sum(DMP,
               DEP,
               DBP,
               BBP,
               DEHP,
               DNOP, na.rm = TRUE))

#total PBDE in ng/g
data<-data %>% rowwise() %>%
  mutate(PBDE =
           sum(BDE17_Kim,
               BDE28_Kim,
               BDE47_Kim,
               BDE49_Kim,
               BDE66_Kim,
               BDE85_Kim,
               BDE99_Kim,
               BDE100_Kim,
               BDE138_Kim,
               BDE153_Kim,
               BDE154_BB153_Kim,
               BDE183_Kim,
               BDE190_Kim,
               BDE209_Kim,
               BDE100_Rob,
               BDE119_Rob,
               BDE138_Rob,
               BDE15_Rob,
               BDE153_Rob,
               BDE154_Rob,
               BDE17_Rob,
               BDE181_Rob,
               BDE183_Rob,
               BDE203_Rob,
               BDE205_Rob,
               BDE206_Rob,
               BDE207_Rob,
               BDE209_Rob,
               BDE28_Rob,
               BDE3_Rob,
               BDE47_Rob,
               BDE49_Rob,
               BDE66_Rob,
               BDE7_Rob,
               BDE71_Rob,
               BDE77_Rob,
               BDE85_BDE155_Rob,
               BDE99_Rob, na.rm = TRUE))

#total HBCDD in ng/g
data<-data %>% rowwise() %>%
  mutate(Total_HBCDD =
           sum(HBCD_Kim,
               alpha_HBCDD_Rob, na.rm = TRUE))

#total DP in ng/g
data<-data %>% rowwise() %>%
  mutate(Total_DP =
           sum(syn_DP_Kim,
               anti_DP_Kim,
               anti_DDCCO_Rob,
               syn_DDCCO_Rob, na.rm = TRUE))


#total NBFRs in ng/g
data<-data %>% rowwise() %>%
  mutate(Total_NBFR =
           sum(b_TBECH_BDE15_Kim,
               BTBPE_Kim,
               HBB_Kim,
               a_TBECH_Kim,
               alpha_beta_DBE_DBCH_Rob,
               BEHTEBP_Rob,
               BTBPE_Rob,
               DBDPE_Rob,
               DBHCTD_Rob,
               EHTBB_Rob,
               HBB_Rob,
               OBTMPI_Rob,
               PBBAcr_Rob,
               PBEB_Rob,
               PBPAE_Rob,
               PBPdbpe_Rob,
               PBT_Rob,
               TBCT_Rob,
               TBPAE_Rob,
               TBPDBPE_Rob,
               TBX_Rob, na.rm = TRUE))

#Total UV stabilizers [ng/g]
data<- data %>% rowwise() %>%
  mutate(Total_UV =
           sum(C4,
               C4C4,
               C8,
               C9,
               C4C8,
               C8C8,
               C9C9,
               diAMS,
               UV234,
               UV326,
               UV327,
               UV328,
               UV329,
               UV350,
               C8_mono_DPA.peak_1, na.rm = TRUE))

#Total metals
data<- data %>% rowwise() %>%
  mutate(Metals =   sum(Lead, Chromium, Arsenic, Cadmium, Copper, Manganese, Rubidium, Aluminum, Mercury, Molybdenum, Nickel, Lithium, Strontium, Boron, Cobalt, Bismuth, Silver, na.rm = TRUE))

#Total PFAS [ng/g]
data<- data %>% rowwise() %>%
  mutate(PFAS =
           sum(FBSA_Rob,
               FOSA_Rob,
               N_MeFOSA_Rob,
               N_EtFOSA_Rob,
               PFEtCHxS_Rob,
               PFBS_Rob,
               PFHxS_Rob,
               PFOS_Rob,
               PFDS_Rob,
               PFBA_Rob,
               PFPeA_Rob,
               PFHxA_Rob,
               PFHpA_Rob,
               PFOA_Rob,
               PFNA_Rob,
               PFDA_Rob,
               PFUdA_Rob,
               PFDoA_Rob,
               PFTrDA_Rob,
               PFTeDA_Rob,
               PFHxDA_Rob,
               PFODA_Rob, na.rm = TRUE))

#Total OPEs [ng/g]
data<- data %>%  rowwise() %>%
  mutate(OPE =
           sum(TEP,
               TPrP,
               TNBP,
               TBOEP,
               TEHP,
               TCEP,
               TCIPP,
               TDCIPP,
               BCMP_BCEP,
               T2B4MP,
               T3B4MP,
               T4B3MP,
               TDBPP,
               TTBNPP,
               TPHP,
               EHDPP,
               TMPP, na.rm = TRUE))

#Total OCPs [ng/g]
data<- data %>%   rowwise() %>%
  mutate(Total_OCP =
           sum(cis_Chlordane,
               trans_Chlordane,
               p.p._DDT,
               p.p._DDD,
               p.p._DDE,
               Dieldrin,
               Heptachlor_epoxide,
               alpha_Hexachlorocyclohexane,
               beta_Hexachlorocyclohehexane,
               gamma_Hexachlorocyclohexane,
               Mirex,
               cis_Nonachlor,
               trans_Nonachlor,
               Octachlorostyrene,
               Oxychlordane,
               Pentachlorobenzene,
               Photomirex,
               x1.2.3.4_Tetrachlorobenzene,
               x1.2.4.5_Tetrachlorobenzene_1.2.3.5_Tetrachlorobenzene, na.rm = TRUE))


# making a new df with only the summarized conaminants

sum_data <- select(data, 1:8, 870:879)


# k-means for all ------------------------------------------------------------


# scaling the data
  # may not be necessary, have to check the manuscript and see what units everything is in

scale_data <- sum_data
scale_data[9:18] <- scale(scale_data[9:18])


# k-means test
  # total contaminants (PDBE, metals, PFAS, OPEs)

kmeans_test <- kmeans(x = scale_data[c(10,15:17)],
                     centers = 5,
                     iter.max = 100,
                     nstart = 25,
                     algorithm = "Hartigan-Wong",
                     trace = FALSE)

kmeans_test$cluster
kmeans_test$centers
kmeans_test$size



#visualizing the clusters

fviz_cluster(kmeans_test, geom = "point", data = scale_data[c(10,15:17)]) + 
  scale_color_manual(values = wesanderson::wes_palette("IsleofDogs1",length(kmeans_test$size))) + 
  scale_fill_manual(values = wesanderson::wes_palette("IsleofDogs1",length(kmeans_test$size))) + 
  labs(title = "Cluster Plot",
       fill = "Cluster",
       shape = "Cluster",
       color = "Cluster") + 
  theme_bw() 


# Elbow method for determining number of clusters
  # 5 clusters is appropriate for the data in its current state

set.seed(69420)
fviz_nbclust(scale_data[c(10,15:17)], kmeans, method = "wss")


# adding a column that shows which row each cluster correlates to

sum_data$cluster <- kmeans_test$cluster
sum_data2 <- sum_data[c(1:8,10,15:17)]


# Pivoting long and facet boxplots ---------------------------------------



# taking the scaled data for only metals, PBDE, OPEs, and PFAS
  # may not be necessary, have to check the manuscript and see what units everything is in

scale_data2 <- scale_data[c(1:8,10,15:17)]


# pivoting the un-scaled data long

long_df <- pivot_longer(sum_data2,
                       cols = c("Metals", "PBDE", "PFAS", "OPE"),
                       names_to = "Contaminant",
                       values_to = "Concentration")

# pivoting the scaled data long

long_df_scale <- pivot_longer(scale_data2,
                        cols = c("Metals", "PBDE", "PFAS", "OPE"),
                        names_to = "Contaminant",
                        values_to = "Concentration")


# facet boxplots template
  # feel free to update the theme/colourings, etc.

ggplot(data = long_df_scale, aes(x = Contaminant, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) +
  facet_wrap( ~ Tissue, scales = "free") + 
  ylab("Concentration") + 
  xlab("Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 10, color = "black"), 
        strip.text.y = element_text(size = 15, color = "black")) +
  theme(text = element_text(size=15), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(color = "black"), 
        legend.position = "bottom")





# k-means for metals ------------------------------------------------------

# attempt to pivot the df wide for metals
  # too many NA values to do any valuable analysis


metals_df <- sum_data[c(1:8,15)]

metals_tissue <- pivot_wider(metals_df,
                             id_cols = c("Bird_ID", "Collection.Location", "Sex", "species"),
                             names_from = "Tissue",
                             values_from = "metals")

metals_location <- pivot_wider(metals_df,
                               id_cols = c("Bird_ID", "Tissue", "Sex", "species"),
                               names_from = "Collection.Location",
                               values_from = "metals")

metals_species <- pivot_wider(metals_df,
                              id_cols = c("Bird_ID", "Collection.Location", "Tissue", "Sex"),
                              names_from = "species",
                              values_from = "metals")

metals_sex <- pivot_wider(metals_df,
                          id_cols = c("Bird_ID", "Collection.Location", "Tissue", "species"),
                          names_from = "Sex",
                          values_from = "metals")



# attempt at doing k-means on one of these pivoted dfs
  # hint: didn't go well

kmeans_metals <- kmeans(x = metals_tissue[8:14],
                      centers = 5,
                      iter.max = 100,
                      nstart = 25,
                      algorithm = "Hartigan-Wong",
                      trace = FALSE)

kmeans_test$cluster
kmeans_test$centers


fviz_cluster(kmeans_metals, geom = "point", data = scale_data[15]) + 
  scale_color_manual(values = wesanderson::wes_palette("IsleofDogs1",length(kmeans_metals$size))) + 
  scale_fill_manual(values = wesanderson::wes_palette("IsleofDogs1",length(kmeans_metals$size))) + 
  labs(title = "Cluster Plot",
       fill = "Cluster",
       shape = "Cluster",
       color = "Cluster") + 
  theme_bw() 


# PCA ---------------------------------------------------------------------


# attempt of doing PCA on the data
  # hint: this also didn't go well

contam_PCA <- PCA(scale_data[c(10,15:17)], scale.unit = TRUE, ncp = 2, graph = TRUE)


fviz_eig(contam_PCA,
         addlabels = TRUE,
         ylim = c(0, 100))



fviz_pca_biplot(contam_PCA, repel = TRUE,
                col.var = "blue",
                col.ind = "red")



# Things to do ------------------------------------------------------------


#k-s test
  #compares distribution
#try other stats test comparing distributions
#compare CVs
#correlation tests between variables
#boxplots for each grouping variable to accompany these tests
#scatterplots comparing different variables, manually colour by a given grouping variable
  #try to manually identify any clusters


