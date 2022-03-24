

# Packages ----------------------------------------------------------------

library(readr)
library(data.table)
library(dplyr)
library(FactoMineR)
library(factoextra)

# Data Manipulations ------------------------------------------------------


data<-fread(paste0(getwd(),"/contaminant_database.csv"))

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
  mutate(Total_PBDE =
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
  mutate(metals =   sum(Lead, Chromium, Arsenic, Cadmium, Copper, Manganese, Rubidium, Aluminum, Mercury, Molybdenum, Nickel, Lithium, Strontium, Boron, Cobalt, Bismuth, Silver, na.rm = TRUE))

#Total PFAS [ng/g]
data<- data %>% rowwise() %>%
  mutate(Total_PFAS =
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
  mutate(Total_OPE =
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

data3 <- select(data, 1:8, 871, 876:878)


# PCA ---------------------------------------------------------------------


contam_PCA <- PCA(data3[9:11], scale.unit = TRUE, ncp = 2, graph = TRUE)


fviz_eig(contam_PCA,
         addlabels = TRUE,
         ylim = c(0, 100))



fviz_pca_biplot(contam_PCA, repel = TRUE,
                col.var = "blue",
                col.ind = "red")



# k-means Clustering ------------------------------------------------------


library(ClusterR)
library(cluster)


reduced_data <- data2[15:16]


k_means_tissue <- kmeans(reduced_data, centers = 8, nstart = 20)
k_means_tissue

k_means_tissue$cluster



