

# Packages ----------------------------------------------------------------

library(readr)
library(data.table)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(tidyr)
#library(Rmisc)
#library(dgof)

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

# use this df

sum_data2 <- sum_data[c(1:8,10,15:17)]


sum_data2$Tissue <- as.factor(sum_data2$Tissue)
sum_data2$species <- as.factor(sum_data2$species)
sum_data2$Sex <- as.factor(sum_data2$Sex)
sum_data2$Collection.Location <- as.factor(sum_data2$Collection.Location)


# scaling the data
  # may not be necessary, have to check the manuscript and see what units everything is in

scale_data <- sum_data
scale_data[9:18] <- scale(scale_data[9:18])

# taking the scaled data for only metals, PBDE, OPEs, and PFAS
  # may not be necessary, have to check the manuscript and see what units everything is in

scale_data2 <- scale_data[c(1:8,10,15:17)]



# replacing non-detects (0's) with randomly generated value that is 
# approximately 0.5*LOD

set.seed(12345)

sum_data2$PBDE[sum_data2$PBDE == 0] <- runif(sum(sum_data2$PBDE == 0),
                                             min = 0.0145000,
                                             max=0.0155000)

sum_data2$Metals[sum_data2$Metals == 0] <- runif(sum(sum_data2$Metals == 0),
                                             min = 0.0145000,
                                             max=0.0155000)

sum_data2$PFAS[sum_data2$PFAS == 0] <- runif(sum(sum_data2$PFAS == 0),
                                             min = 0.0145000,
                                             max=0.0155000)

sum_data2$OPE[sum_data2$OPE == 0] <- runif(sum(sum_data2$OPE == 0),
                                             min = 0.0145000,
                                             max=0.0155000)




# Pivoting long and facet boxplots ---------------------------------------


# pivoting the un-scaled data long

long_df <- pivot_longer(sum_data2,
                       cols = c("Metals", "PBDE", "PFAS", "OPE"),
                       names_to = "Contaminant",
                       values_to = "Concentration")


# subsetting by contaminant class

long_df_metals <- filter(long_df, Contaminant == "Metals")

long_df_OPE <- filter(long_df, Contaminant == "OPE")

long_df_PBDE <- filter(long_df, Contaminant == "PBDE")

long_df_PFAS <- filter(long_df, Contaminant == "PFAS")





# pivoting the long_df to a wide one for statistical analysis
  # pivotting by tissue

tissue_df <- pivot_wider(long_df,
                             id_cols = c("Bird_ID",
                                         "Element",
                                         "USOX",
                                         "species",
                                         "Sex",
                                         "Collection.Date",
                                         "Collection.Location",
                                         "Contaminant"),
                             names_from = "Tissue",
                             values_from = "Concentration")


tissue_metals <- subset(tissue_df, Contaminant == "Metals")
tissue_PBDE <- subset(tissue_df, Contaminant == "PBDE")
tissue_PFAS <- subset(tissue_df, Contaminant == "PFAS")
tissue_OPE <- subset(tissue_df, Contaminant == "OPE")

# pivotting by location

location_df <- pivot_wider(long_df,
                         id_cols = c("Bird_ID",
                                     "Element",
                                     "USOX",
                                     "species",
                                     "Sex",
                                     "Collection.Date",
                                     "Tissue",
                                     "Contaminant"),
                         names_from = "Collection.Location",
                         values_from = "Concentration")


location_metals <- subset(location_df, Contaminant == "Metals")
location_PBDE <- subset(location_df, Contaminant == "PBDE")
location_PFAS <- subset(location_df, Contaminant == "PFAS")
location_OPE <- subset(location_df, Contaminant == "OPE")


# pivotting by species

species_df <- pivot_wider(long_df,
                           id_cols = c("Bird_ID",
                                       "Element",
                                       "USOX",
                                       "Collection.Location",
                                       "Sex",
                                       "Collection.Date",
                                       "Tissue",
                                       "Contaminant"),
                           names_from = "species",
                           values_from = "Concentration")


species_metals <- subset(species_df, Contaminant == "Metals")
species_PBDE <- subset(species_df, Contaminant == "PBDE")
species_PFAS <- subset(species_df, Contaminant == "PFAS")
species_OPE <- subset(species_df, Contaminant == "OPE")


# pivotting by sex

sex_df <- pivot_wider(long_df,
                          id_cols = c("Bird_ID",
                                      "Element",
                                      "USOX",
                                      "Collection.Location",
                                      "species",
                                      "Collection.Date",
                                      "Tissue",
                                      "Contaminant"),
                          names_from = "Sex",
                          values_from = "Concentration")


sex_metals <- subset(sex_df, Contaminant == "Metals")
sex_PBDE <- subset(sex_df, Contaminant == "PBDE")
sex_PFAS <- subset(sex_df, Contaminant == "PFAS")
sex_OPE <- subset(sex_df, Contaminant == "OPE")




# facet boxplots template
# feel free to update the theme/colourings, etc.

ggplot(data = long_df, aes(x = Contaminant, y = Concentration)) + 
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



# Median, Min, Max --------------------------------------------------------

basic_stats <- function(vect){
  
  values1 <- as.vector(vect)
  values1 <- na.omit(values1)
  
  med <- median(values1)
  maximum <- max(values1)
  minimum <- min(values1)
  sample_mean <- mean(values1)
  
  return(c(minimum, med, maximum, sample_mean))
  
}


# basic stats for tissue variable

basic_stats(tissue_metals$preen_oil)
basic_stats(tissue_metals$liver)
basic_stats(tissue_metals$egg)
basic_stats(tissue_metals$muscle)
basic_stats(tissue_metals$brain)
basic_stats(tissue_metals$blood)
basic_stats(tissue_metals$fat)

basic_stats(tissue_OPE$preen_oil)
basic_stats(tissue_OPE$liver)
basic_stats(tissue_OPE$egg)
basic_stats(tissue_OPE$muscle)
basic_stats(tissue_OPE$brain)
basic_stats(tissue_OPE$blood)
basic_stats(tissue_OPE$fat)

basic_stats(tissue_PBDE$preen_oil)
basic_stats(tissue_PBDE$liver)
basic_stats(tissue_PBDE$egg)
basic_stats(tissue_PBDE$muscle)
basic_stats(tissue_PBDE$brain)
basic_stats(tissue_PBDE$blood)
basic_stats(tissue_PBDE$fat)

basic_stats(tissue_PFAS$preen_oil)
basic_stats(tissue_PFAS$liver)
basic_stats(tissue_PFAS$egg)
basic_stats(tissue_PFAS$muscle)
basic_stats(tissue_PFAS$brain)
basic_stats(tissue_PFAS$blood)
basic_stats(tissue_PFAS$fat)

# basic stats for species variable

basic_stats(species_metals$NOFU)
basic_stats(species_metals$BLKI)
basic_stats(species_OPE$NOFU)
basic_stats(species_OPE$BLKI)
basic_stats(species_PBDE$NOFU)
basic_stats(species_PBDE$BLKI)
basic_stats(species_PFAS$NOFU)
basic_stats(species_PFAS$BLKI)


# basic stats for location variable

basic_stats(location_metals$`Labrador Sea`)
basic_stats(location_metals$`Prince Leopold Island, NU`)
basic_stats(location_OPE$`Labrador Sea`)
basic_stats(location_OPE$`Prince Leopold Island, NU`)
basic_stats(location_PBDE$`Labrador Sea`)
basic_stats(location_PBDE$`Prince Leopold Island, NU`)
basic_stats(location_PFAS$`Labrador Sea`)
basic_stats(location_PFAS$`Prince Leopold Island, NU`)


# basic stats for sex variable

basic_stats(sex_metals$F)
basic_stats(sex_metals$M)
basic_stats(sex_OPE$F)
basic_stats(sex_OPE$M)
basic_stats(sex_PBDE$F)
basic_stats(sex_PBDE$M)
basic_stats(sex_PFAS$F)
basic_stats(sex_PFAS$M)


# ks-test tissue -------------------------------------------------------------

# KS-test for tissues

ks_tissue_matrix <- matrix(data = NA, nrow = 7, ncol = 7)
colnames(ks_tissue_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")
rownames(ks_tissue_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")


ks_tissue <- function(df){
  
  for (i in rownames(ks_tissue_matrix)) {
    
    for (j in colnames(ks_tissue_matrix)) {
      
      subset_df1 <- filter(df, Tissue == i)
      subset_df2 <- filter(df, Tissue == j)
      
      ks_tissue_matrix[i,j] <- ks.test(x = subset_df1$Concentration,
                                       y = subset_df2$Concentration)$p.value
    }
  }
  
  return(ks_tissue_matrix)
  
}

metals_ks_matrix <- ks_tissue(long_df_metals)
write.csv(metals_ks_matrix, "metals_ks_matrix.csv")

OPE_ks_matrix <- ks_tissue(long_df_OPE)
write.csv(OPE_ks_matrix, "OPE_ks_matrix.csv")

PBDE_ks_matrix <- ks_tissue(long_df_PBDE)
write.csv(PBDE_ks_matrix, "PBDE_ks_matrix.csv")

PFAS_ks_matrix <- ks_tissue(long_df_PFAS)
write.csv(PFAS_ks_matrix, "PFAS_ks_matrix.csv")






# ks-test location --------------------------------------------------------


ks_location_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(ks_location_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")
rownames(ks_location_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")


ks_location <- function(df){
  
  for (i in rownames(ks_location_matrix)) {
    
    for (j in colnames(ks_location_matrix)) {
      
      subset_df1 <- filter(df, Collection.Location == i)
      subset_df2 <- filter(df, Collection.Location == j)
      
      ks_location_matrix[i,j] <- ks.test(x = subset_df1$Concentration,
                                         y = subset_df2$Concentration)$p.value
    }
  }
  
  return(ks_location_matrix)
  
}


metals_ks_location <- ks_location(long_df_metals)
write.csv(metals_ks_location, "metals_ks_location.csv")

OPE_ks_location <- ks_location(long_df_OPE)
write.csv(OPE_ks_location, "OPE_ks_location.csv")

PBDE_ks_location <- ks_location(long_df_PBDE)
write.csv(PBDE_ks_location, "PBDE_ks_location.csv")

PFAS_ks_location <- ks_location(long_df_PFAS)
write.csv(PFAS_ks_location, "PFAS_ks_location.csv")


# ks-test species ---------------------------------------------------------


ks_species_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(ks_species_matrix) <- c("BLKI", "NOFU")
rownames(ks_species_matrix) <- c("BLKI", "NOFU")


ks_species <- function(df){
  
  for (i in rownames(ks_species_matrix)) {
    
    for (j in colnames(ks_species_matrix)) {
      
      subset_df1 <- filter(df, species == i)
      subset_df2 <- filter(df, species == j)
      
      ks_species_matrix[i,j] <- ks.test(x = subset_df1$Concentration,
                                        y = subset_df2$Concentration)$p.value
    }
  }
  
  return(ks_species_matrix)
  
}


# saving each matrix as a CSV

metals_ks_species <- ks_species(long_df_metals)
write.csv(metals_ks_species, "metals_ks_species.csv")

OPE_ks_species <- ks_species(long_df_OPE)
write.csv(OPE_ks_species, "OPE_ks_species.csv")

PBDE_ks_species <- ks_species(long_df_PBDE)
write.csv(PBDE_ks_species, "PBDE_ks_species.csv")

PFAS_ks_species <- ks_species(long_df_PFAS)
write.csv(PFAS_ks_species, "PFAS_ks_species.csv")


# ks-test sex -------------------------------------------------------------


ks_sex_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(ks_sex_matrix) <- c("F", "M")
rownames(ks_sex_matrix) <- c("F", "M")


ks_sex <- function(df){
  
  for (i in rownames(ks_sex_matrix)) {
    
    for (j in colnames(ks_sex_matrix)) {
      
      subset_df1 <- filter(df, Sex == i)
      subset_df2 <- filter(df, Sex == j)
      
      ks_sex_matrix[i,j] <- ks.test(x = subset_df1$Concentration,
                                    y = subset_df2$Concentration)$p.value
    }
  }
  
  return(ks_sex_matrix)
  
}


metals_ks_sex <- ks_sex(long_df_metals)
write.csv(metals_ks_sex, "metals_ks_sex.csv")

OPE_ks_sex <- ks_sex(long_df_OPE)
write.csv(OPE_ks_sex, "OPE_ks_sex.csv")

PBDE_ks_sex <- ks_sex(long_df_PBDE)
write.csv(PBDE_ks_sex, "PBDE_ks_sex.csv")

PFAS_ks_sex <- ks_sex(long_df_PFAS)
write.csv(PFAS_ks_sex, "PFAS_ks_sex.csv")


# CVs ---------------------------------------------------------------------

# CVs for metals in tissue

(sd(tissue_metals$blood, na.rm = TRUE)/mean(tissue_metals$blood, na.rm = TRUE))*100
(sd(tissue_metals$brain, na.rm = TRUE)/mean(tissue_metals$brain, na.rm = TRUE))*100
(sd(tissue_metals$egg, na.rm = TRUE)/mean(tissue_metals$egg, na.rm = TRUE))*100
(sd(tissue_metals$fat, na.rm = TRUE)/mean(tissue_metals$fat, na.rm = TRUE))*100
(sd(tissue_metals$liver, na.rm = TRUE)/mean(tissue_metals$liver, na.rm = TRUE))*100
(sd(tissue_metals$muscle, na.rm = TRUE)/mean(tissue_metals$muscle, na.rm = TRUE))*100
(sd(tissue_metals$preen_oil, na.rm = TRUE)/mean(tissue_metals$preen_oil, na.rm = TRUE))*100


# CVs for OPEs in tissue

(sd(tissue_OPE$blood, na.rm = TRUE)/mean(tissue_OPE$blood, na.rm = TRUE))*100
(sd(tissue_OPE$brain, na.rm = TRUE)/mean(tissue_OPE$brain, na.rm = TRUE))*100
(sd(tissue_OPE$egg, na.rm = TRUE)/mean(tissue_OPE$egg, na.rm = TRUE))*100
(sd(tissue_OPE$fat, na.rm = TRUE)/mean(tissue_OPE$fat, na.rm = TRUE))*100
(sd(tissue_OPE$liver, na.rm = TRUE)/mean(tissue_OPE$liver, na.rm = TRUE))*100
(sd(tissue_OPE$muscle, na.rm = TRUE)/mean(tissue_OPE$muscle, na.rm = TRUE))*100
(sd(tissue_OPE$preen_oil, na.rm = TRUE)/mean(tissue_OPE$preen_oil, na.rm = TRUE))*100


# CVs for PBDEs in tissue

(sd(tissue_PBDE$blood, na.rm = TRUE)/mean(tissue_PBDE$blood, na.rm = TRUE))*100
(sd(tissue_PBDE$brain, na.rm = TRUE)/mean(tissue_PBDE$brain, na.rm = TRUE))*100
(sd(tissue_PBDE$egg, na.rm = TRUE)/mean(tissue_PBDE$egg, na.rm = TRUE))*100
(sd(tissue_PBDE$fat, na.rm = TRUE)/mean(tissue_PBDE$fat, na.rm = TRUE))*100
(sd(tissue_PBDE$liver, na.rm = TRUE)/mean(tissue_PBDE$liver, na.rm = TRUE))*100
(sd(tissue_PBDE$muscle, na.rm = TRUE)/mean(tissue_PBDE$muscle, na.rm = TRUE))*100
(sd(tissue_PBDE$preen_oil, na.rm = TRUE)/mean(tissue_PBDE$preen_oil, na.rm = TRUE))*100


# CVs for PFAS in tissue

(sd(tissue_PFAS$blood, na.rm = TRUE)/mean(tissue_PFAS$blood, na.rm = TRUE))*100
(sd(tissue_PFAS$brain, na.rm = TRUE)/mean(tissue_PFAS$brain, na.rm = TRUE))*100
(sd(tissue_PFAS$egg, na.rm = TRUE)/mean(tissue_PFAS$egg, na.rm = TRUE))*100
(sd(tissue_PFAS$fat, na.rm = TRUE)/mean(tissue_PFAS$fat, na.rm = TRUE))*100
(sd(tissue_PFAS$liver, na.rm = TRUE)/mean(tissue_PFAS$liver, na.rm = TRUE))*100
(sd(tissue_PFAS$muscle, na.rm = TRUE)/mean(tissue_PFAS$muscle, na.rm = TRUE))*100
(sd(tissue_PFAS$preen_oil, na.rm = TRUE)/mean(tissue_PFAS$preen_oil, na.rm = TRUE))*100






# CVs for metals by location

(sd(location_metals$`Labrador Sea`, na.rm = TRUE)/mean(location_metals$`Labrador Sea`, na.rm = TRUE))*100
(sd(location_metals$`Prince Leopold Island, NU`, na.rm = TRUE)/mean(location_metals$`Prince Leopold Island, NU`, na.rm = TRUE))*100


# CVs for OPEs by location

(sd(location_OPE$`Labrador Sea`, na.rm = TRUE)/mean(location_OPE$`Labrador Sea`, na.rm = TRUE))*100
(sd(location_OPE$`Prince Leopold Island, NU`, na.rm = TRUE)/mean(location_OPE$`Prince Leopold Island, NU`, na.rm = TRUE))*100


# CVs for PBDE by location

(sd(location_PBDE$`Labrador Sea`, na.rm = TRUE)/mean(location_PBDE$`Labrador Sea`, na.rm = TRUE))*100
(sd(location_PBDE$`Prince Leopold Island, NU`, na.rm = TRUE)/mean(location_PBDE$`Prince Leopold Island, NU`, na.rm = TRUE))*100


# CVs for PFAS by location

(sd(location_PFAS$`Labrador Sea`, na.rm = TRUE)/mean(location_PFAS$`Labrador Sea`, na.rm = TRUE))*100
(sd(location_PFAS$`Prince Leopold Island, NU`, na.rm = TRUE)/mean(location_PFAS$`Prince Leopold Island, NU`, na.rm = TRUE))*100





# CVs for metals by species

(sd(species_metals$BLKI, na.rm = TRUE)/mean(species_metals$BLKI, na.rm = TRUE))*100
(sd(species_metals$NOFU, na.rm = TRUE)/mean(species_metals$NOFU, na.rm = TRUE))*100


# CVs for OPEs by species

(sd(species_OPE$BLKI, na.rm = TRUE)/mean(species_OPE$BLKI, na.rm = TRUE))*100
(sd(species_OPE$NOFU, na.rm = TRUE)/mean(species_OPE$NOFU, na.rm = TRUE))*100


# CVs for PBDE by species

(sd(species_PBDE$BLKI, na.rm = TRUE)/mean(species_PBDE$BLKI, na.rm = TRUE))*100
(sd(species_PBDE$NOFU, na.rm = TRUE)/mean(species_PBDE$NOFU, na.rm = TRUE))*100

# CVs for PFAS by species

(sd(species_PFAS$BLKI, na.rm = TRUE)/mean(species_PFAS$BLKI, na.rm = TRUE))*100
(sd(species_PFAS$NOFU, na.rm = TRUE)/mean(species_PFAS$NOFU, na.rm = TRUE))*100




# CVs for metals by sex

(sd(sex_metals$F, na.rm = TRUE)/mean(sex_metals$F, na.rm = TRUE))*100
(sd(sex_metals$M, na.rm = TRUE)/mean(sex_metals$M, na.rm = TRUE))*100


# CVs for OPEs by sex

(sd(sex_OPE$F, na.rm = TRUE)/mean(sex_OPE$F, na.rm = TRUE))*100
(sd(sex_OPE$M, na.rm = TRUE)/mean(sex_OPE$M, na.rm = TRUE))*100


# CVs for PBDE by sex

(sd(sex_PBDE$F, na.rm = TRUE)/mean(sex_PBDE$F, na.rm = TRUE))*100
(sd(sex_PBDE$M, na.rm = TRUE)/mean(sex_PBDE$M, na.rm = TRUE))*100

# CVs for PFAS by sex

(sd(sex_PFAS$F, na.rm = TRUE)/mean(sex_PFAS$F, na.rm = TRUE))*100
(sd(sex_PFAS$M, na.rm = TRUE)/mean(sex_PFAS$M, na.rm = TRUE))*100


# PCA ---------------------------------------------------------------------


# attempt of doing PCA on the data
  # hint: this also didn't go well

contam_PCA <- PCA(sum_data2[c(9:12)], scale.unit = TRUE, ncp = 2, graph = TRUE)


fviz_eig(contam_PCA,
         addlabels = TRUE,
         ylim = c(0, 30))



fviz_pca_biplot(contam_PCA, repel = TRUE,
                col.var = "blue",
                col.ind = "red")


# k-means for all ------------------------------------------------------------


# correlation matrix
# shows there are no strong correlations between any of the contaminants

cor(sum_data2[9:12],
    method = "spearman")


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

#set.seed(69420)
fviz_nbclust(scale_data[c(10,15:17)], kmeans, method = "wss")


# adding a column that shows which row each cluster correlates to

sum_data$cluster <- kmeans_test$cluster


