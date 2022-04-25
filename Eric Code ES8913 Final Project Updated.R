

# Packages ----------------------------------------------------------------


set.seed(12345)

library(readr)
library(data.table)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
source(paste0(getwd(), "./utils.R"))

# Data Import and Manipulation -------------------------------------------------------

data <- read.csv(paste0(getwd(), "./contaminant_database.csv"))

# Change name of Sex variable

for (i in 1:dim(data)[1]) {
  if (data$Sex[i] == "F") {
    data$Sex[i] <- "Female"
  }
  else {
    data$Sex[i] <- "Male"
  }
}

# Adding contaminants into summary columns -----------------------------------

#total Ph
data <- data %>% rowwise() %>%
  mutate(Total_Ph =
           sum(DMP,
               DEP,
               DBP,
               BBP,
               DEHP,
               DNOP, na.rm = TRUE))

#total PBDE in ng/g ***
data <- data %>% rowwise() %>%
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
data <- data %>% 
  rowwise() %>%
  mutate(Total_HBCDD =
           sum(HBCD_Kim,
               alpha_HBCDD_Rob, na.rm = TRUE))

#total DP in ng/g
data <- data %>% 
  rowwise() %>%
  mutate(Total_DP =
           sum(syn_DP_Kim,
               anti_DP_Kim,
               anti_DDCCO_Rob,
               syn_DDCCO_Rob, na.rm = TRUE))


#total NBFRs in ng/g
data <- data %>% 
  rowwise() %>%
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
data <- data %>% 
  rowwise() %>%
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

#Total metals ***
data <- data %>% 
  rowwise() %>%
  mutate(Metals = sum(Lead, 
                        Chromium, 
                        Arsenic,
                        Cadmium, 
                        Copper, 
                        Manganese, 
                        Rubidium, 
                        Aluminum, 
                        Mercury, 
                        Molybdenum, 
                        Nickel, 
                        Lithium, 
                        Strontium, 
                        Boron, 
                        Cobalt,
                        Bismuth, 
                        Silver, na.rm = TRUE))

#Total PFAS [ng/g] ***
data <- data %>% 
  rowwise() %>%
  mutate(PFAS =
           sum(FBSA._Rob,
               FOSA._Rob,
               N_MeFOSA._Rob,
               N_EtFOSA._Rob,
               PFEtCHxS._Rob,
               PFBS._Rob,
               PFHxS._Rob,
               PFOS._Rob,
               PFDS._Rob,
               PFBA._Rob,
               PFPeA._Rob,
               PFHxA._Rob,
               PFHpA._Rob,
               PFOA._Rob,
               PFNA._Rob,
               PFDA._Rob,
               PFUdA._Rob,
               PFDoA._Rob,
               PFTrDA._Rob,
               PFTeDA._Rob,
               PFHxDA._Rob,
               PFODA._Rob, na.rm = TRUE))

#Total OPEs [ng/g] ***
data <- data %>%  
  rowwise() %>%
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
data <- data %>%   
  rowwise() %>%
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
               X1.2.3.4_Tetrachlorobenzene,
               X1.2.4.5_Tetrachlorobenzene_1.2.3.5_Tetrachlorobenzene, na.rm = TRUE))


# Extracting 4 best groups of contaminant: metals, PBDE, PFAS, OPE --------


subset_data <- select(data, Tissue, species, Sex, Collection.Location, 
                      PBDE, Metals, PFAS, OPE)


# Making each grouping variable (tissue, species, sex, location) a factor


subset_data$Tissue <- as.factor(subset_data$Tissue)
subset_data$species <- as.factor(subset_data$species)
subset_data$Sex <- as.factor(subset_data$Sex)
subset_data$Collection.Location <- as.factor(subset_data$Collection.Location)

# Replace zero with NA values for plotting scatter plot -----------------------------------------------------------
subset_data_NA <- zero_it_out(subset_data, 
                              col = c("PBDE", "Metals", "PFAS", "OPE"),
                              options = "by_NA")

# Non-detect (zero) value replacement -------------------------------------


# All zero values were replaced with a random value around the half-LOD (0.015)


# replacing zeros in the Metals group

subset_data$Metals[subset_data$Metals == 0]<-runif(sum(subset_data$Metals == 0),
                                                   min = 0.0145000,
                                                   max=0.0155000)


# replacing zeros in the OPEs group

subset_data$OPE[subset_data$OPE == 0]<-runif(sum(subset_data$OPE == 0),
                                             min = 0.0145000,
                                             max=0.0155000)


# replacing zeros in the PBDEs group

subset_data$PBDE[subset_data$PBDE == 0]<-runif(sum(subset_data$PBDE == 0),
                                               min = 0.0145000,
                                               max=0.0155000)


# replacing zeros in the PFAS group

subset_data$PFAS[subset_data$PFAS == 0]<-runif(sum(subset_data$PFAS == 0),
                                               min = 0.0145000,
                                               max=0.0155000)


# checking for duplicate values in each of the variables with 0's replaced
# duplicate returns a TRUE/FALSE for each duplicate, so by taking the sum
# any value greater than 0 will mean there are still ties in the data

sum(duplicated(subset_data$Metals))
sum(duplicated(subset_data$OPE))
sum(duplicated(subset_data$PBDE))
sum(duplicated(subset_data$PFAS))


# Pivoting the Data Long -------------------------------------------------


# pivoting the data long by contaminant type

long_df <- pivot_longer(subset_data,
                        cols = c("Metals", "PBDE", "PFAS", "OPE"),
                        names_to = "Contaminant",
                        values_to = "Concentration")


# sub-setting by contaminant class

long_df_metals <- filter(long_df, Contaminant == "Metals")

long_df_OPE <- filter(long_df, Contaminant == "OPE")

long_df_PBDE <- filter(long_df, Contaminant == "PBDE")

long_df_PFAS <- filter(long_df, Contaminant == "PFAS")


# PCA ---------------------------------------------------------------------


# Attempt of doing PCA on the four contaminant groups


contaminant_PCA <- PCA(subset_data[c(5:8)],
                       scale.unit = TRUE,
                       ncp = 2,
                       graph = TRUE)


# visualizing the Scree plot of the PCA

fviz_eig(contaminant_PCA,
         addlabels = TRUE,
         ylim = c(0, 30))


# correlation matrix to identify any linear relationships
# between the contaminants

cor(subset_data[5:8],
    method = "spearman")


# k-means clustering ---------------------------------------------------------


# Elbow method for determining number of clusters
# 5 clusters is appropriate for the data


fviz_nbclust(scale(subset_data[5:8]), kmeans, method = "wss")


# k-means clustering on the total contaminant levels (PDBE, metals, PFAS, OPEs)


kmeans_contaminants <- kmeans(x = scale(subset_data[5:8]),
                      centers = 5,
                      iter.max = 100,
                      nstart = 25,
                      algorithm = "Hartigan-Wong",
                      trace = FALSE)


# identifying the cluster centres

kmeans_contaminants$centers


# identifying the cluster sizes

kmeans_contaminants$size


#visualizing the clusters

fviz_cluster(kmeans_contaminants,
             geom = "point",
             data = scale(subset_data[5:8])) + 
  labs(title = "K-means Cluster Plot",
       fill = "Cluster",
       shape = "Cluster",
       color = "Cluster") + 
  theme_bw()


# Median, Min, Max, Mean, CV -------------------------------------------------


# calculating the minimum, maximum, median, mean, and CVs for each
# contaminant, grouped by tissue type

tissue_stats <- rbind((subset_data %>% 
                          group_by(Tissue) %>%
                          summarise("Minimum (ng/g ww)" = min(Metals),
                                    "Median (ng/g ww)" = median(Metals),
                                    "Maximum (ng/g ww)" = max(Metals),
                                    "Mean (ng/g ww)" = mean(Metals),
                                    "CV (%)" = (sd(Metals)/mean(Metals))*100)),
                       
                       (subset_data %>% 
                          group_by(Tissue) %>%
                          summarise("Minimum (ng/g ww)" = min(OPE),
                                    "Median (ng/g ww)" = median(OPE),
                                    "Maximum (ng/g ww)" = max(OPE),
                                    "Mean (ng/g ww)" = mean(OPE),
                                    "CV (%)" = (sd(OPE)/mean(OPE))*100)),          
                       
                       (subset_data %>% 
                          group_by(Tissue) %>%
                          summarise("Minimum (ng/g ww)" = min(PBDE),
                                    "Median (ng/g ww)" = median(PBDE),
                                    "Maximum (ng/g ww)" = max(PBDE),
                                    "Mean (ng/g ww)" = mean(PBDE),
                                    "CV (%)" = (sd(PBDE)/mean(PBDE))*100)),
                       
                       (subset_data %>% 
                          group_by(Tissue) %>%
                          summarise("Minimum (ng/g ww)" = min(PFAS),
                                    "Median (ng/g ww)" = median(PFAS),
                                    "Maximum (ng/g ww)" = max(PFAS),
                                    "Mean (ng/g ww)" = mean(PFAS),
                                    "CV (%)" = (sd(PFAS)/mean(PFAS))*100)))

tissue_stats <- add_column(tissue_stats,
                            Contaminant = c(rep("Metals", 7),
                                            rep("OPE", 7),
                                            rep("PBDE", 7),
                                            rep("PFAS", 7)),
                            .before = 1)


# calculating the minimum, maximum, median, mean, and CVs for each
# contaminant, grouped by bird species

species_stats <- rbind((subset_data %>% 
                           group_by(species) %>%
                           summarise("Minimum (ng/g ww)" = min(Metals),
                                     "Median (ng/g ww)" = median(Metals),
                                     "Maximum (ng/g ww)" = max(Metals),
                                     "Mean (ng/g ww)" = mean(Metals),
                                     "CV (%)" = (sd(Metals)/mean(Metals))*100)),
                        
                        (subset_data %>% 
                           group_by(species) %>%
                           summarise("Minimum (ng/g ww)" = min(OPE),
                                     "Median (ng/g ww)" = median(OPE),
                                     "Maximum (ng/g ww)" = max(OPE),
                                     "Mean (ng/g ww)" = mean(OPE),
                                     "CV (%)" = (sd(OPE)/mean(OPE))*100)),          
                        
                        (subset_data %>% 
                           group_by(species) %>%
                           summarise("Minimum (ng/g ww)" = min(PBDE),
                                     "Median (ng/g ww)" = median(PBDE),
                                     "Maximum (ng/g ww)" = max(PBDE),
                                     "Mean (ng/g ww)" = mean(PBDE),
                                     "CV (%)" = (sd(PBDE)/mean(PBDE))*100)),
                        
                        (subset_data %>% 
                           group_by(species) %>%
                           summarise("Minimum (ng/g ww)" = min(PFAS),
                                     "Median (ng/g ww)" = median(PFAS),
                                     "Maximum (ng/g ww)" = max(PFAS),
                                     "Mean (ng/g ww)" = mean(PFAS),
                                     "CV (%)" = (sd(PFAS)/mean(PFAS))*100)))

species_stats <- add_column(species_stats,
                             Contaminant = c("Metals", "Metals",
                                             "OPE", "OPE",
                                             "PBDE", "PBDE",
                                             "PFAS", "PFAS"),
                             .before = 1)


# calculating the minimum, maximum, median, mean, and CVs for each
# contaminant, grouped by location

location_stats <- rbind((subset_data %>% 
  group_by(Collection.Location) %>%
  summarise("Minimum (ng/g ww)" = min(Metals),
            "Median (ng/g ww)" = median(Metals),
            "Maximum (ng/g ww)" = max(Metals),
            "Mean (ng/g ww)" = mean(Metals),
            "CV (%)" = (sd(Metals)/mean(Metals))*100)),

(subset_data %>% 
  group_by(Collection.Location) %>%
  summarise("Minimum (ng/g ww)" = min(OPE),
            "Median (ng/g ww)" = median(OPE),
            "Maximum (ng/g ww)" = max(OPE),
            "Mean (ng/g ww)" = mean(OPE),
            "CV (%)" = (sd(OPE)/mean(OPE))*100)),          

(subset_data %>% 
  group_by(Collection.Location) %>%
  summarise("Minimum (ng/g ww)" = min(PBDE),
            "Median (ng/g ww)" = median(PBDE),
            "Maximum (ng/g ww)" = max(PBDE),
            "Mean (ng/g ww)" = mean(PBDE),
            "CV (%)" = (sd(PBDE)/mean(PBDE))*100)),

(subset_data %>% 
  group_by(Collection.Location) %>%
  summarise("Minimum (ng/g ww)" = min(PFAS),
            "Median (ng/g ww)" = median(PFAS),
            "Maximum (ng/g ww)" = max(PFAS),
            "Mean (ng/g ww)" = mean(PFAS),
            "CV (%)" = (sd(PFAS)/mean(PFAS))*100)))

location_stats <- add_column(location_stats,
                             Contaminant = c("Metals", "Metals",
                                             "OPE", "OPE",
                                             "PBDE", "PBDE",
                                             "PFAS", "PFAS"),
                             .before = 1)



# calculating the minimum, maximum, median, mean, and CVs for each
# contaminant, grouped by bird sex

sex_stats <- rbind((subset_data %>% 
                      group_by(Sex) %>%
  summarise("Minimum (ng/g ww)" = min(Metals),
            "Median (ng/g ww)" = median(Metals),
            "Maximum (ng/g ww)" = max(Metals),
            "Mean (ng/g ww)" = mean(Metals),
            "CV (%)" = (sd(Metals)/mean(Metals))*100)),

(subset_data %>% 
  group_by(Sex) %>%
  summarise("Minimum (ng/g ww)" = min(OPE),
            "Median (ng/g ww)" = median(OPE),
            "Maximum (ng/g ww)" = max(OPE),
            "Mean (ng/g ww)" = mean(OPE),
            "CV (%)" = (sd(OPE)/mean(OPE))*100)),           

(subset_data %>% 
  group_by(Sex) %>%
  summarise("Minimum (ng/g ww)" = min(PBDE),
            "Median (ng/g ww)" = median(PBDE),
            "Maximum (ng/g ww)" = max(PBDE),
            "Mean (ng/g ww)" = mean(PBDE),
            "CV (%)" = (sd(PBDE)/mean(PBDE))*100)),

(subset_data %>% 
  group_by(Sex) %>%
  summarise("Minimum (ng/g ww)" = min(PFAS),
            "Median (ng/g ww)" = median(PFAS),
            "Maximum (ng/g ww)" = max(PFAS),
            "Mean (ng/g ww)" = mean(PFAS),
            "CV (%)" = (sd(PFAS)/mean(PFAS))*100)))

sex_stats <- add_column(sex_stats,
                        Contaminant = c("Metals", "Metals",
                                        "OPE", "OPE",
                                        "PBDE", "PBDE",
                                        "PFAS", "PFAS"),
                        .before = 1)


# ks-test tissue -------------------------------------------------------------

# KS-test for tissues


# Creating an empty matrix to store the ks-test p-values for tissue variables

ks_tissue_matrix <- matrix(data = NA, nrow = 7, ncol = 7)

colnames(ks_tissue_matrix) <- c("blood",
                                "brain",
                                "egg",
                                "fat",
                                "liver",
                                "muscle",
                                "preen_oil")

rownames(ks_tissue_matrix) <- c("blood",
                                "brain",
                                "egg",
                                "fat",
                                "liver",
                                "muscle",
                                "preen_oil")


# creating a function, ks_tissue, to run the ks-test comparing two samples
# (grouped by tissue type), and store the p-value in the ks_tissue_matrix

# function takes one of the long contaminant data frames as an argument (df)

ks_tissue <- function(df){
  
  # iterates through every tissue name in the rows
  
  for (i in rownames(ks_tissue_matrix)) {
    
    # iterates through every tissue name in the columns
    
    for (j in colnames(ks_tissue_matrix)) {
      
      # creates two data frames for the two tissue types being compared
      
      subset_df1 <- filter(df, Tissue == i)
      subset_df2 <- filter(df, Tissue == j)
      
      # runs the ks-test for the contaminant levels in the two tissue types
      
      ks_tissue_matrix[i,j] <- ks.test(x = subset_df1$Concentration,
                                       y = subset_df2$Concentration)$p.value
    }
  }
  
  # returns the matrix filled with all the ks-test p-values
  
  return(ks_tissue_matrix)
  
}

tissue_KS <- as.data.frame(rbind(ks_tissue(long_df_metals),
                   ks_tissue(long_df_OPE),
                   ks_tissue(long_df_PBDE),
                   ks_tissue(long_df_PFAS)))

tissue_KS <- add_column(tissue_KS,
                           Contaminant = c(rep("Metals", 7),
                                           rep("OPE", 7),
                                           rep("PBDE", 7),
                                           rep("PFAS", 7)),
                           .before = 1)



# ks-test location --------------------------------------------------------


# Creating an empty matrix to store the ks-test p-values for location variables

ks_location_matrix <- matrix(data = NA, nrow = 2, ncol = 2)

colnames(ks_location_matrix) <- c("Labrador Sea",
                                  "Prince Leopold Island, NU")

rownames(ks_location_matrix) <- c("Labrador Sea",
                                  "Prince Leopold Island, NU")


# creating a function, ks_location, to run the ks-test comparing two samples
# (grouped by location), and store the p-value in the ks_location_matrix

# function takes one of the long contaminant data frames as an argument (df)

ks_location <- function(df){
  
  # iterates through every location name in the rows
  
  for (i in rownames(ks_location_matrix)) {
    
    # iterates through every location name in the columns
    
    for (j in colnames(ks_location_matrix)) {
      
      # creates two data frames for the two locations being compared
      
      subset_df1 <- filter(df, Collection.Location == i)
      subset_df2 <- filter(df, Collection.Location == j)
      
      # runs the ks-test for the contaminant levels in the two locations 
      
      ks_location_matrix[i,j] <- ks.test(x = subset_df1$Concentration,
                                         y = subset_df2$Concentration)$p.value
    }
  }
  
  # returns the matrix filled with all the ks-test p-values
  
  return(ks_location_matrix)
  
}


location_KS <- as.data.frame(rbind(ks_location(long_df_metals),
                                   ks_location(long_df_OPE),
                                   ks_location(long_df_PBDE),
                                   ks_location(long_df_PFAS)))


location_KS <- add_column(location_KS,
                         Contaminant = c("Metals", "Metals",
                                         "OPE", "OPE",
                                         "PBDE", "PBDE",
                                         "PFAS", "PFAS"),
                         .before = 1)


# ks-test species ---------------------------------------------------------


# Creating an empty matrix to store the ks-test p-values for species variables

ks_species_matrix <- matrix(data = NA, nrow = 2, ncol = 2)

colnames(ks_species_matrix) <- c("BLKI",
                                 "NOFU")

rownames(ks_species_matrix) <- c("BLKI",
                                 "NOFU")


# creating a function, ks_species, to run the ks-test comparing two samples
# (grouped by species), and store the p-value in the ks_species_matrix

# function takes one of the long contaminant data frames as an argument (df)

ks_species <- function(df){
  
  # iterates through every species name in the rows
  
  for (i in rownames(ks_species_matrix)) {
    
    # iterates through every species name in the columns
    
    for (j in colnames(ks_species_matrix)) {
      
      # creates two data frames for the two species being compared
      
      subset_df1 <- filter(df, species == i)
      subset_df2 <- filter(df, species == j)
      
      # runs the ks-test for the contaminant levels in the two species 
      
      ks_species_matrix[i,j] <- ks.test(x = subset_df1$Concentration,
                                        y = subset_df2$Concentration)$p.value
    }
  }
  
  # returns the matrix filled with all the ks-test p-values
  
  return(ks_species_matrix)
  
}

species_KS <- as.data.frame(rbind(ks_species(long_df_metals),
                                  ks_species(long_df_OPE),
                                  ks_species(long_df_PBDE),
                                  ks_species(long_df_PFAS)))


species_KS <- add_column(species_KS,
                            Contaminant = c("Metals", "Metals",
                                            "OPE", "OPE",
                                            "PBDE", "PBDE",
                                            "PFAS", "PFAS"),
                            .before = 1)

# ks-test sex -------------------------------------------------------------


# Creating an empty matrix to store the ks-test p-values for sex variables

ks_sex_matrix <- matrix(data = NA, nrow = 2, ncol = 2)

colnames(ks_sex_matrix) <- c("Female",
                             "Male")

rownames(ks_sex_matrix) <- c("Female",
                             "Male")


# creating a function, ks_sex, to run the ks-test comparing two samples
# (grouped by sex), and store the p-value in the ks_sex_matrix

# function takes one of the long contaminant data frames as an argument (df)

ks_sex <- function(df){
  
  # iterates through every sex in the rows
  
  for (i in rownames(ks_sex_matrix)) {
    
    # iterates through every sex in the columns
    
    for (j in colnames(ks_sex_matrix)) {
      
      # creates two data frames for the two sexes being compared
      
      subset_df1 <- filter(df, Sex == i)
      subset_df2 <- filter(df, Sex == j)
      
      # runs the ks-test for the contaminant levels in the two sexes 
      
      ks_sex_matrix[i,j] <- ks.test(x = subset_df1$Concentration,
                                    y = subset_df2$Concentration)$p.value
    }
  }
  
  # returns the matrix filled with all the ks-test p-values
  
  return(ks_sex_matrix)
  
}

metals_ks_sex <- ks_sex(long_df_metals)
OPE_ks_sex <- ks_sex(long_df_OPE)
PBDE_ks_sex <- ks_sex(long_df_PBDE)
PFAS_ks_sex <- ks_sex(long_df_PFAS)



sex_KS <- as.data.frame(rbind(metals_ks_sex,
                              OPE_ks_sex,
                              PBDE_ks_sex,
                              PFAS_ks_sex))


sex_KS <- add_column(sex_KS,
                          Contaminant = c("Metals", "Metals",
                                          "OPE", "OPE",
                                          "PBDE", "PBDE",
                                          "PFAS", "PFAS"),
                          .before = 1)

# Mann-Whitney test -------------------------------------------------------

# MW Tissue -------------------------------------------------------
mw_tissue_matrix <- matrix(data = NA, nrow = 7, ncol = 7)

colnames(mw_tissue_matrix) <- c("blood",
                                "brain",
                                "egg",
                                "fat",
                                "liver",
                                "muscle",
                                "preen_oil")

rownames(mw_tissue_matrix) <- c("blood",
                                "brain",
                                "egg",
                                "fat",
                                "liver",
                                "muscle",
                                "preen_oil")

mw_tissue <- function(df){
  
  # iterates through every tissue name in the rows
  
  for (i in rownames(mw_tissue_matrix)) {
    
    # iterates through every tissue name in the columns
    
    for (j in colnames(mw_tissue_matrix)) {
      
      # creates two data frames for the two tissue types being compared
      
      subset_df1 <- filter(df, Tissue == i)
      subset_df2 <- filter(df, Tissue == j)
      
      # runs the ks-test for the contaminant levels in the two tissue types
      
      mw_tissue_matrix[i,j] <- wilcox.test(x = subset_df1$Concentration,
                                           y = subset_df2$Concentration)$p.value
    }
  }
  
  # returns the matrix filled with all the ks-test p-values
  
  return(mw_tissue_matrix)
  
}

tissue_MW <- as.data.frame(rbind(mw_tissue(long_df_metals),
                                 mw_tissue(long_df_OPE),
                                 mw_tissue(long_df_PBDE),
                                 mw_tissue(long_df_PFAS)))

tissue_MW <- add_column(tissue_MW,
                        Contaminant = c(rep("Metals", 7),
                                        rep("OPE", 7),
                                        rep("PBDE", 7),
                                        rep("PFAS", 7)),
                        .before = 1)


# MW Species --------------------------------------------------------

mw_species_matrix <- matrix(data = NA, nrow = 2, ncol = 2)

colnames(mw_species_matrix) <- c("BLKI",
                                 "NOFU")

rownames(mw_species_matrix) <- c("BLKI",
                                 "NOFU")

mw_species <- function(df){
  
  # iterates through every species name in the rows
  
  for (i in rownames(mw_species_matrix)) {
    
    # iterates through every species name in the columns
    
    for (j in colnames(mw_species_matrix)) {
      
      # creates two data frames for the two species being compared
      
      subset_df1 <- filter(df, species == i)
      subset_df2 <- filter(df, species == j)
      
      # runs the ks-test for the contaminant levels in the two species 
      
      mw_species_matrix[i,j] <- wilcox.test(x = subset_df1$Concentration,
                                            y = subset_df2$Concentration)$p.value
    }
  }
  
  # returns the matrix filled with all the ks-test p-values
  
  return(mw_species_matrix)
  
}

species_MW <- as.data.frame(rbind(mw_species(long_df_metals),
                                  mw_species(long_df_OPE),
                                  mw_species(long_df_PBDE),
                                  mw_species(long_df_PFAS)))


species_MW <- add_column(species_MW,
                         Contaminant = c("Metals", "Metals",
                                         "OPE", "OPE",
                                         "PBDE", "PBDE",
                                         "PFAS", "PFAS"),
                         .before = 1)



# MW Sex------------------------------------------------------------

mw_sex_matrix <- matrix(data = NA, nrow = 2, ncol = 2)

colnames(mw_sex_matrix) <- c("Male",
                             "Female")

rownames(mw_sex_matrix) <- c("Male",
                             "Female")


# creating a function, ks_sex, to run the ks-test comparing two samples
# (grouped by sex), and store the p-value in the ks_sex_matrix

# function takes one of the long contaminant data frames as an argument (df)

mw_sex <- function(df){
  
  # iterates through every sex in the rows
  
  for (i in rownames(mw_sex_matrix)) {
    
    # iterates through every sex in the columns
    
    for (j in colnames(mw_sex_matrix)) {
      
      # creates two data frames for the two sexes being compared
      
      subset_df1 <- filter(df, Sex == i)
      subset_df2 <- filter(df, Sex == j)
      
      # runs the ks-test for the contaminant levels in the two sexes 
      
      mw_sex_matrix[i,j] <- wilcox.test(x = subset_df1$Concentration,
                                    y = subset_df2$Concentration)$p.value
    }
  }
  
  # returns the matrix filled with all the ks-test p-values
  
  return(mw_sex_matrix)
  
}

metals_mw_sex <- mw_sex(long_df_metals)
OPE_mw_sex <- mw_sex(long_df_OPE)
PBDE_mw_sex <- mw_sex(long_df_PBDE)
PFAS_mw_sex <- mw_sex(long_df_PFAS)



sex_MW <- as.data.frame(rbind(metals_mw_sex,
                              OPE_mw_sex,
                              PBDE_mw_sex,
                              PFAS_mw_sex))

sex_MW <- add_column(sex_MW,
                     Contaminant = c("Metals", "Metals",
                                     "OPE", "OPE",
                                     "PBDE", "PBDE",
                                     "PFAS", "PFAS"),
                     .before = 1)

# MW Location-------------------------------------------------------
# Creating an empty matrix to store the ks-test p-values for location variables

mw_location_matrix <- matrix(data = NA, nrow = 2, ncol = 2)

colnames(mw_location_matrix) <- c("Labrador Sea",
                                  "Prince Leopold Island, NU")

rownames(mw_location_matrix) <- c("Labrador Sea",
                                  "Prince Leopold Island, NU")

mw_location <- function(df){
  
  # iterates through every location name in the rows
  
  for (i in rownames(mw_location_matrix)) {
    
    # iterates through every location name in the columns
    
    for (j in colnames(mw_location_matrix)) {
      
      # creates two data frames for the two locations being compared
      
      subset_df1 <- filter(df, Collection.Location == i)
      subset_df2 <- filter(df, Collection.Location == j)
      
      # runs the ks-test for the contaminant levels in the two locations 
      
      mw_location_matrix[i,j] <- wilcox.test(x = subset_df1$Concentration,
                                         y = subset_df2$Concentration)$p.value
    }
  }
  
  # returns the matrix filled with all the ks-test p-values
  
  return(mw_location_matrix)
  
}


location_MW <- as.data.frame(rbind(mw_location(long_df_metals),
                                   mw_location(long_df_OPE),
                                   mw_location(long_df_PBDE),
                                   mw_location(long_df_PFAS)))


location_MW <- add_column(location_MW,
                          Contaminant = c("Metals", "Metals",
                                          "OPE", "OPE",
                                          "PBDE", "PBDE",
                                          "PFAS", "PFAS"),
                          .before = 1)

# Heat Maps ---------------------------------------------------------------


library(pheatmap)

tissue_nondetect <- c((14/14)*100,(31/31)*100,(0/11)*100,(6/6)*100,(10/31)*100,(27/27)*100,(10/10)*100,
                      (14/14)*100,(22/31)*100,(11/11)*100,(0/6)*100,(31/31)*100,(16/27)*100,(10/10)*100,
                      (14/14)*100,(29/31)*100,(0/11)*100,(0/6)*100,(20/31)*100,(18/27)*100,(10/10)*100,
                      (14/14)*100,(31/31)*100,(11/11)*100,(6/6)*100,(21/31)*100,(27/27)*100,(10/10)*100)

tissue_col_name <- c("Blood", "Brain", "Egg", "Fat", "Liver", "Muscle", "Preen Oil")
row_names_heat <- c("Metals", "OPEs", "PBDEs", "PFAS")



tissue_nondetect_matrix <- matrix(data = tissue_nondetect,
                                  nrow = 4,
                                  byrow = TRUE,
                                  dimnames = list(row_names_heat, tissue_col_name))


pheatmap(tissue_nondetect_matrix,
         display_numbers = T,
         color = colorRampPalette(c('white','red'))(100),
         breaks = seq(0,100,1),
         cluster_rows = F,
         cluster_cols = F,
         fontsize_number = 15,
         number_color = "black",
         number_format = "%g",
         fontsize = 15,
         angle_col = 0,
         annotation_legend = TRUE)

species_nondetect <- c(82,64,
                       67,100,
                       69,74,
                       88,100)

species_colnames <- c("Northern Fulmar", "Black-Legged Kittiwake")

species_nondetect_matrix <- matrix(data = species_nondetect,
                                  nrow = 4,
                                  byrow = TRUE,
                                  dimnames = list(row_names_heat, species_colnames))


colour_gradient <- colorRampPalette(c("white", "red"))(100)


pheatmap(species_nondetect_matrix,
         display_numbers = T,
         color = colour_gradient,
         breaks = seq(0,100,1),
         cluster_rows = F,
         cluster_cols = F,
         fontsize_number = 15,
         number_color = "black",
         number_format = "%i",
         fontsize = 15,
         angle_col = 0,
         annotation_legend = TRUE)

location_nondetect <- c(100,64,
                        48,94,
                       74,69,
                       86,95)

location_colnames <- c("Labrador Sea", "Prince Leopold Island")

location_nondetect_matrix <- matrix(data = location_nondetect,
                                   nrow = 4,
                                   byrow = TRUE,
                                   dimnames = list(row_names_heat, location_colnames))


colour_gradient <- colorRampPalette(c("white", "red"))(100)

pheatmap(location_nondetect_matrix,
         display_numbers = T,
         color = colour_gradient,
         breaks = seq(0,100,1),
         cluster_rows = F,
         cluster_cols = F,
         fontsize_number = 15,
         number_color = "black",
         number_format = "%g",
         fontsize = 15,
         angle_col = 0,
         annotation_legend = TRUE)

sex_nondetect <- c(72,77,
                   100,71,
                   92,63,
                   89,94)

sex_colnames <- c("Male", "Female")

sex_nondetect_matrix <- matrix(data = sex_nondetect,
                                    nrow = 4,
                                    byrow = TRUE,
                                    dimnames = list(row_names_heat,
                                                    sex_colnames))


colour_gradient <- colorRampPalette(c("white", "red"))(100)

pheatmap(sex_nondetect_matrix,
         display_numbers = T,
         color = colour_gradient,
         breaks = seq(0,100,1),
         cluster_rows = F,
         cluster_cols = F,
         fontsize_number = 15,
         number_color = "black",
         number_format = "%i",
         fontsize = 15,
         angle_col = 0,
         annotation_legend = TRUE)

# Box plots Tissue with non-detects -------------------------------------------------
tissue_boxplot <- ggplot(data = long_df, 
                         aes(x = Tissue, y = Concentration)) + 
  facet_wrap(~Contaminant, scales = "free_y") +
  geom_boxplot(aes(fill = Contaminant)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Tissue Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 20,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 20,
                                    vjust = 1),
        axis.text.x = element_text(size = 15,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden",
        plot.margin = margin(t = 0.7, r = 0.7, b = 0.7, l = 0.7, "cm"))

# Box plots Location with non-detects---------------------------------------------------------

location_boxplot <- ggplot(data = long_df, aes(x = Collection.Location, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) +
  facet_wrap( ~ Contaminant, scales = "free") + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Location") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 20,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 20,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden",
        plot.margin = margin(t = 0.7, r = 0.7, b = 0.7, l = 0.7, "cm"))

# Box plots Sex with non-detects-----------------------------------------------------

sex_boxplot <- ggplot(data = long_df, aes(x = Sex, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) +
  facet_wrap( ~ Contaminant, scales = "free") + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Sex") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 20,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 20,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden",
        plot.margin = margin(t = 0.7, r = 0.7, b = 0.7, l = 0.7, "cm"))

# Box plots Species with non-detects-------------------------------------------------

spe_boxplot <- ggplot(data = long_df, aes(x = species, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) +
  facet_wrap( ~ Contaminant, scales = "free") + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Species") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 20,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 20,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden",
        plot.margin = margin(t = 0.7, r = 0.7, b = 0.7, l = 0.7, "cm"))


# Scatter plot Tissue with vs. without non-detects------------------------------------------------------------------

metal_PBDE <- ggplot(subset_data, aes(x = Metals, 
                                      y = PBDE)) +
  geom_point(aes(color = Tissue)) +
  labs(title = "metal vs. PBDE with non-detects")

metal_PBDE_clean_scatter <- ggplot(subset_data_NA, aes(x = Metals, 
                                                       y = PBDE)) +
  geom_point(aes(color = Tissue)) +
  labs(title = "metal vs. PBDE without non-detects")

scatter_plot_tissue <- ggarrange(metal_PBDE, metal_PBDE_clean_scatter, 
                                 ncol=2, 
                                 common.legend = TRUE,
                                 legend = "right")

# Scatter plot Species with vs without non-detects-----------------------------------------------------------------
metal_PBDE <- ggplot(subset_data, aes(x = Metals, 
                                      y = PBDE)) +
  geom_point(aes(color = species)) +
  labs(title = "metal vs. PBDE with non-detects")

metal_PBDE_clean <- ggplot(subset_data_NA, aes(x = Metals, 
                                               y = PBDE)) +
  geom_point(aes(color = species)) +
  labs(title = "metal vs. PBDE without non-detects")

scatter_plot_species <- ggarrange(metal_PBDE, metal_PBDE_clean, 
                                  ncol=2, 
                                  common.legend = TRUE,
                                  legend = "right")

# Scatter plot Location with vs. without non-detects---------------------------------------------------------------------

OPE_PBDE <- ggplot(subset_data, aes(x = OPE,
                                         y = PBDE)) +
  geom_point(aes(color = Collection.Location)) +
  labs(title = "Total_OPE vs. PBDE with non-detects")

OPE_PBDE_clean <- ggplot(subset_data_NA, aes(x = OPE,
                                                    y = PBDE)) +
  geom_point(aes(color = Collection.Location)) +
  labs(title = "Total_OPE vs. PBDE without non-detects")

scatter_plot_location <- ggarrange(OPE_PBDE, OPE_PBDE_clean, 
                                   ncol=2, 
                                   common.legend = TRUE,
                                   legend = "right")

# Scatter plot Sex with vs without non-detects----------------------------------------------------------------
PBDE_PFAS <- ggplot(subset_data, aes(x = PBDE,
                                          y = PFAS)) +
  geom_point(aes(color = Sex)) +
  labs(title = "PBDE vs. PFAS with non-detects")

PBDE_PFAS_clean <- ggplot(subset_data_NA, aes(x = PBDE,
                                                     y = PFAS)) +
  geom_point(aes(color = Sex)) +
  labs(title = "PBDE vs. PFAS without non-detects")

scatter_plot_sex <- ggarrange(PBDE_PFAS, PBDE_PFAS_clean, 
                              ncol=2, 
                              common.legend = TRUE,
                              legend = "right")
