# Loading packages --------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(data.table)
library(viridis)
library(wesanderson)
library(scales)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(vegan)
library(rmarkdown) 
library(knitr)
library(FactoMineR)
library(factoextra)
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

# Adding sum columns ------------------------------------------------------

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
  mutate(metals =   sum(Lead, 
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
  mutate(Total_PFAS =
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
                           Total_PBDE, metals, Total_PFAS, Total_OPE)

# Replace zero with NA values for plotting scatter plot
subset_data_NA <- zero_it_out(subset_data, 
                                col = c("Total_PBDE", "metals", "Total_PFAS", "Total_OPE"),
                                options = "by_NA")

# Option 1. Replace zero values with a value that is 0.5*LOD  -------------
subset_data_LOD <- zero_it_out(subset_data, 
                               col = c("Total_PBDE", "metals", "Total_PFAS", "Total_OPE"),
                               options = "LOD")

# Option 2. REMOVING all zero values --------------------------------------
sum_4_conta_data_metals_clean <- sum_4_conta_data %>%
  filter(metals > 0) 
sum_4_conta_data_PBDE_clean <- sum_4_conta_data %>%
  filter(Total_PBDE > 0) 
sum_4_conta_data_PFAS_clean <- sum_4_conta_data %>%
  filter(Total_PFAS > 0) 
sum_4_conta_data_OPE_clean <- sum_4_conta_data %>%
  filter(Total_OPE > 0) 

  
# Scale data to reduce outlier effect -------------------------------------
sum_4_conta_scale_data <- copy(sum_4_conta_data)
sum_4_conta_scale_data[5:8] <- scale(sum_4_conta_scale_data[5:8])

# Pivot longer --------------------------------------------------------------
long_df <- pivot_longer(sum_4_conta_data_LOD,
                        cols = c("metals", "Total_PBDE", 
                                 "Total_PFAS", "Total_OPE"),
                        names_to = "Contaminant",
                        values_to = "Concentration")

# With non-detects, aka. LOD data
long_df_metals <- filter(long_df, Contaminant == "metals")
long_df_PBDE <- filter(long_df, Contaminant == "Total_PBDE")
long_df_PFAS <- filter(long_df, Contaminant == "Total_PFAS")
long_df_OPE <- filter(long_df, Contaminant == "Total_OPE")

# Without non-detects data
long_df_metals_clean<- pivot_longer(sum_4_conta_data_metals_clean,
                                    cols = c("metals", "Total_PBDE", 
                                             "Total_PFAS", "Total_OPE"),
                                    names_to = "Contaminant",
                                    values_to = "Concentration")

long_df_PBDE_clean <- pivot_longer(sum_4_conta_data_PBDE_clean,
                                   cols = c("metals", "Total_PBDE", 
                                            "Total_PFAS", "Total_OPE"),
                                   names_to = "Contaminant",
                                   values_to = "Concentration")

long_df_PFAS_clean <- pivot_longer(sum_4_conta_data_PFAS_clean,
                                   cols = c("metals", "Total_PBDE", 
                                            "Total_PFAS", "Total_OPE"),
                                   names_to = "Contaminant",
                                   values_to = "Concentration")

long_df_OPE_clean <- pivot_longer(sum_4_conta_data_OPE_clean,
                                  cols = c("metals", "Total_PBDE", 
                                           "Total_PFAS", "Total_OPE"),
                                  names_to = "Contaminant",
                                  values_to = "Concentration")


long_scale_df <- pivot_longer(sum_4_conta_scale_data,
                        cols = c("metals", "Total_PBDE", 
                                 "Total_PFAS", "Total_OPE"),
                        names_to = "Contaminant",
                        values_to = "Concentration")


# Visualize Data Distribution ---------------------------------------------

# Facet box plots Tissue with non-detects -------------------------------------------------
tissue_boxplot_metals <- ggplot(data = filter(long_df, Contaminant == "metals"), 
                                aes(x = Tissue, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 25, color = "black"),
        axis.title.x = element_text(size = 25,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 25,
                                    vjust = 1),
        axis.text.x = element_text(size = 25,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 25,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=25), 
        legend.text = element_text(size=25),
        legend.position = "hidden")

tissue_boxplot_PBDE <- ggplot(data = filter(long_df, Contaminant == "Total_PBDE"), 
                              aes(x = Tissue, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 25, color = "black"),
        axis.title.x = element_text(size = 25,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 25,
                                    vjust = 1),
        axis.text.x = element_text(size = 25,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 25,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=25), 
        legend.text = element_text(size=25),
        legend.position = "hidden")

tissue_boxplot_PFAS <- ggplot(data = filter(long_df, Contaminant == "Total_PFAS"), 
                              aes(x = Tissue, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) +
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 25, color = "black"),
        axis.title.x = element_text(size = 25,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 25,
                                    vjust = 1),
        axis.text.x = element_text(size = 25,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 25,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=25), 
        legend.text = element_text(size=25),
        legend.position = "hidden")

tissue_boxplot_OPE <- ggplot(data = filter(long_df, Contaminant == "Total_OPE"), 
                                aes(x = Tissue, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 25, color = "black"),
        axis.title.x = element_text(size = 25,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 25,
                                    vjust = 1),
        axis.text.x = element_text(size = 25,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 25,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=25), 
        legend.text = element_text(size=25),
        legend.position = "hidden")

# Exporting plot 
ggsave(paste0(getwd(), "/tissue_boxplot_metals.png"), 
       tissue_boxplot_metals,
       dpi = 480, 
       height = 10)
ggsave(paste0(getwd(), "/tissue_boxplot_PBDE.png"), 
       tissue_boxplot_PBDE,
       dpi = 480, 
       height = 10)
ggsave(paste0(getwd(), "/tissue_boxplot_PFAS.png"), 
       tissue_boxplot_PFAS,
       dpi = 480, 
       height = 10)
ggsave(paste0(getwd(), "/tissue_boxplot_OPE.png"), 
       tissue_boxplot_OPE,
       dpi = 480, 
       height = 10)



# Facet box plots Tissue without non-detects --------------------------------
tissue_boxplot_metals_clean <- ggplot(data = filter(long_df_metals_clean, 
                                                    Contaminant == "metals"), 
                                aes(x = Tissue, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 25, color = "black"),
        axis.title.x = element_text(size = 25,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 25,
                                    vjust = 1),
        axis.text.x = element_text(size = 25,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 25,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=25), 
        legend.text = element_text(size=25),
        legend.position = "hidden")

tissue_boxplot_PBDE_clean <- ggplot(data = filter(long_df_PBDE_clean, 
                                                  Contaminant == "Total_PBDE"), 
                              aes(x = Tissue, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 25, color = "black"),
        axis.title.x = element_text(size = 25,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 25,
                                    vjust = 1),
        axis.text.x = element_text(size = 25,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 25,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=25), 
        legend.text = element_text(size=25),
        legend.position = "hidden")

tissue_boxplot_PFAS_clean <- ggplot(data = filter(long_df_PFAS_clean, 
                                                  Contaminant == "Total_PFAS"), 
                              aes(x = Tissue, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 25, color = "black"),
        axis.title.x = element_text(size = 25,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 25,
                                    vjust = 1),
        axis.text.x = element_text(size = 25,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 25,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=25), 
        legend.text = element_text(size=25),
        legend.position = "hidden")

tissue_boxplot_OPE_clean <- ggplot(data = filter(long_df_OPE_clean, 
                                                 Contaminant == "Total_OPE"), 
                             aes(x = Tissue, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 25, color = "black"),
        axis.title.x = element_text(size = 25,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 25,
                                    vjust = 1),
        axis.text.x = element_text(size = 25,
                                   angle = 45, 
                                   hjust = 1), 
        axis.text.y = element_text(size = 25,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=25), 
        legend.text = element_text(size=25),
        legend.position = "hidden")

# Exporting plot 
ggsave(paste0(getwd(), "/tissue_boxplot_metals_clean.png"), 
       tissue_boxplot_metals_clean,
       dpi = 480, 
       height = 10)
ggsave(paste0(getwd(), "/tissue_boxplot_PBDE_clean.png"), 
       tissue_boxplot_PBDE_clean,
       dpi = 480, 
       height = 10)
ggsave(paste0(getwd(), "/tissue_boxplot_PFAS_clean.png"), 
       tissue_boxplot_PFAS_clean,
       dpi = 480, 
       height = 10)
ggsave(paste0(getwd(), "/tissue_boxplot_OPE_clean.png"), 
       tissue_boxplot_OPE_clean,
       dpi = 480, 
       height = 10)

# Facet box plots Location with non-detects---------------------------------------------------------

location_boxplot <- ggplot(data = long_df, aes(x = Collection.Location, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) +
  facet_wrap( ~ Contaminant, scales = "free") + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Location") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

# Exporting plot 
ggsave(paste0(getwd(), "/location_boxplot_2.png"), 
       location_boxplot,
       dpi = 480, 
       height = 10, 
       width = 10)


# Facet box plots Location without non-detects ----------------------------
location_boxplot_metal_clean <- ggplot(data = filter(long_df_metals_clean, 
                                                            Contaminant == "metals"), 
                                       aes(x = Collection.Location, 
                                               y = Concentration)) + 
  geom_boxplot(aes(fill = Collection.Location)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Location",
       title = "Metals") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

location_boxplot_PBDE_clean <- ggplot(data = filter(long_df_PBDE_clean, 
                                                      Contaminant == "Total_PBDE"), 
                                       aes(x = Collection.Location, 
                                           y = Concentration)) + 
  geom_boxplot(aes(fill = Collection.Location)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Location",
       title = "PBDE") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

location_boxplot_PFAS_clean <- ggplot(data = filter(long_df_PFAS_clean, 
                                                            Contaminant == "Total_PFAS"), 
                                       aes(x = Collection.Location, 
                                           y = Concentration)) + 
  geom_boxplot(aes(fill = Collection.Location)) +
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Location",
       title = "PFAS") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

location_boxplot_OPE_clean <- ggplot(data = filter(long_df_OPE_clean, 
                                                            Contaminant == "Total_OPE"), 
                                       aes(x = Collection.Location, 
                                           y = Concentration)) + 
  geom_boxplot(aes(fill = Collection.Location)) +
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Location",
       title = "OPE") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

location_boxplot_clean <- ggarrange(location_boxplot_metal_clean,
                                    location_boxplot_OPE_clean,
                                    location_boxplot_PBDE_clean,
                                    location_boxplot_PFAS_clean,
                                    ncol = 2, 
                                    nrow = 2)

ggsave(paste0(getwd(), "/location_boxplot_clean.png"), 
       location_boxplot_clean,
       dpi = 480, 
       height = 10, 
       width = 10)

# Facet box plots Sex with non-detects-----------------------------------------------------

sex_boxplot <- ggplot(data = long_df, aes(x = Sex, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) +
  facet_wrap( ~ Contaminant, scales = "free") + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Contaminant Type") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

# Exporting plot 
ggsave(paste0(getwd(), "/sex_boxplot_2.png"), 
       sex_boxplot,
       dpi = 320,
       height = 10)


# Facet box plots Sex without non-detects ---------------------------------
sex_boxplot_metal_clean <- ggplot(data = filter(long_df_metals_clean, 
                                                     Contaminant == "metals"), 
                                       aes(x = Sex, 
                                           y = Concentration)) + 
  geom_boxplot(aes(fill = Sex)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Sex",
       title = "Metals") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

sex_boxplot_PBDE_clean <- ggplot(data = filter(long_df_PBDE_clean, 
                                                    Contaminant == "Total_PBDE"), 
                                      aes(x = Sex, 
                                          y = Concentration)) + 
  geom_boxplot(aes(fill = Sex)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Sex",
       title = "PBDE") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

sex_boxplot_PFAS_clean <- ggplot(data = filter(long_df_PFAS_clean, 
                                                    Contaminant == "Total_PFAS"), 
                                      aes(x = Sex, 
                                          y = Concentration)) + 
  geom_boxplot(aes(fill = Sex)) +
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Sex",
       title = "PFAS") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

sex_boxplot_OPE_clean <- ggplot(data = filter(long_df_OPE_clean, 
                                                   Contaminant == "Total_OPE"), 
                                     aes(x = Sex, 
                                         y = Concentration)) + 
  geom_boxplot(aes(fill = Sex)) +
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "Sex",
       title = "OPE") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

sex_boxplot_clean <- ggarrange(sex_boxplot_metal_clean,
                               sex_boxplot_OPE_clean,
                                    sex_boxplot_PBDE_clean,
                                    sex_boxplot_PFAS_clean,
                                    ncol = 2, 
                                    nrow = 2)

ggsave(paste0(getwd(), "/sex_boxplot_clean.png"), 
       sex_boxplot_clean,
       dpi = 480, 
       height = 10, 
       width = 10)

# Facet box plots Species with non-detects-------------------------------------------------

spe_boxplot <- ggplot(data = long_df, aes(x = species, y = Concentration)) + 
  geom_boxplot(aes(fill = Contaminant)) +
  facet_wrap( ~ Contaminant, scales = "free") + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "species") +
  theme_classic() + 
  theme(strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

# Exporting plot 
ggsave(paste0(getwd(), "/spe_boxplot_2.png"), 
       spe_boxplot,
       dpi = 480,
       height = 10)


# Facet box plots Species without non-detects -----------------------------
species_boxplot_metal_clean <- ggplot(data = filter(long_df_metals_clean, 
                                                Contaminant == "metals"), 
                                  aes(x = species, 
                                      y = Concentration)) + 
  geom_boxplot(aes(fill = species)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "species",
       title = "Metals") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

species_boxplot_PBDE_clean <- ggplot(data = filter(long_df_PBDE_clean, 
                                               Contaminant == "Total_PBDE"), 
                                 aes(x = species, 
                                     y = Concentration)) + 
  geom_boxplot(aes(fill = species)) + 
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "species",
       title = "PBDE") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

species_boxplot_PFAS_clean <- ggplot(data = filter(long_df_PFAS_clean, 
                                               Contaminant == "Total_PFAS"), 
                                 aes(x = species, 
                                     y = Concentration)) + 
  geom_boxplot(aes(fill = species)) +
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "species",
       title = "PFAS") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

species_boxplot_OPE_clean <- ggplot(data = filter(long_df_OPE_clean, 
                                              Contaminant == "Total_OPE"), 
                                aes(x = species, 
                                    y = Concentration)) + 
  geom_boxplot(aes(fill = species)) +
  labs(y = bquote("Concentration (ng g"^-1*"ww)"),
       x = "species",
       title = "OPE") +
  theme_classic() + 
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_text(size = 15,
                                    vjust = -0.5),
        axis.title.y = element_text(size = 15,
                                    vjust = 1),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15,
                                   hjust = 1,
                                   vjust = 0.5,
                                   color = "black"), 
        legend.title = element_text(size=15), 
        legend.text = element_text(size=15),
        legend.position = "hidden")

species_boxplot_clean <- ggarrange(species_boxplot_metal_clean,
                               species_boxplot_OPE_clean,
                               species_boxplot_PBDE_clean,
                               species_boxplot_PFAS_clean,
                               ncol = 2, 
                               nrow = 2)

ggsave(paste0(getwd(), "/species_boxplot_clean.png"), 
       species_boxplot_clean,
       dpi = 480, 
       height = 10, 
       width = 10)

# Scatter plot Tissue, Species, Sex, Location ------------------------------------

# Tissue with vs. without non-detects------------------------------------------------------------------
metal_PBDE <- ggplot(sum_4_conta_data_LOD, aes(x = metals, 
                                           y = Total_PBDE)) +
  geom_point(aes(color = Tissue)) +
  labs(title = "metal vs. PBDE with non-detects")

metal_PBDE_clean_scatter_tissue <- ggplot(sum_4_conta_data_copy, aes(x = metals, 
                                                              y = Total_PBDE)) +
  geom_point(aes(color = Tissue)) +
  labs(title = "metal vs. PBDE without non-detects")

scatter_plot_tissue <- ggarrange(metal_PBDE, metal_PBDE_clean_scatter, 
                          ncol=2, 
                          common.legend = TRUE,
                          legend = "right")

ggsave(paste0(getwd(), "/metal_PBDE_clean_scatter_tissue.png"), 
       metal_PBDE_clean_scatter_tissue,
       dpi = 480)


# Species with vs without non-detects-----------------------------------------------------------------
metal_PBDE <- ggplot(sum_4_conta_data, aes(x = metals, 
                                           y = Total_PBDE)) +
  geom_point(aes(color = species)) +
  labs(title = "metal vs. PBDE with non-detects")

metal_PBDE_clean_scatter_species <- ggplot(sum_4_conta_data_copy, aes(x = metals, 
                                                      y = Total_PBDE)) +
  geom_point(aes(color = species)) +
  labs(title = "metal vs. PBDE without non-detects")

scatter_plot_species <- ggarrange(metal_PBDE, metal_PBDE_clean, 
                                 ncol=2, 
                                 common.legend = TRUE,
                                 legend = "right")

ggsave(paste0(getwd(), "/metal_PBDE_clean_scatter_species.png"), 
       metal_PBDE_clean_scatter_species,
       dpi = 480)

# Location with vs. without non-detects---------------------------------------------------------------------

OPE_PBDE <- ggplot(sum_4_conta_data, aes(x = Total_OPE,
                                                    y = Total_PBDE)) +
  geom_point(aes(color = Collection.Location)) +
  labs(title = "Total_OPE vs. PBDE with non-detects")

OPE_PBDE_clean <- ggplot(sum_4_conta_data_copy, aes(x = Total_OPE,
                                                    y = Total_PBDE)) +
  geom_point(aes(color = Collection.Location)) +
  labs(title = "Total_OPE vs. PBDE without non-detects")

scatter_plot_location <- ggarrange(OPE_PBDE, OPE_PBDE_clean, 
                              ncol=2, 
                              common.legend = TRUE,
                              legend = "right")

ggsave(paste0(getwd(), "/scatter_plot_location.png"), 
       scatter_plot_location,
       dpi = 480)

# Sex with vs without non-detects----------------------------------------------------------------
PBDE_PFAS <- ggplot(sum_4_conta_data, aes(x = Total_PBDE,
                                          y = Total_PFAS)) +
  geom_point(aes(color = Sex)) +
  labs(title = "PBDE vs. PFAS with non-detects")

PBDE_PFAS_clean <- ggplot(sum_4_conta_data_copy, aes(x = Total_PBDE,
                                          y = Total_PFAS)) +
  geom_point(aes(color = Sex)) +
  labs(title = "PBDE vs. PFAS without non-detects")



scatter_plot_sex <- ggarrange(PBDE_PFAS, PBDE_PFAS_clean, 
                                   ncol=2, 
                                   common.legend = TRUE,
                                   legend = "right")

ggsave(paste0(getwd(), "/scatter_plot_sex.png"), 
       scatter_plot_sex,
       dpi = 480)

# Bar plot frequency of non-detects  --------------------------------------
PBDE_nondetect <- sum(is.na(sum_4_conta_data_copy$Total_PBDE))/130
metals_nondetect <- sum(is.na(sum_4_conta_data_copy$metals))/130
PFAS_nondetect <- sum(is.na(sum_4_conta_data_copy$Total_PFAS))/130
OPE_nondetect <- sum(is.na(sum_4_conta_data_copy$Total_OPE))/130
nondetect_df <- data.frame(Contaminant = c("metals", "PBDE", "PFAS", "OPE"),
                           Nondetects_quantity = c(metals_nondetect,
                                                   PBDE_nondetect,
                                                   PFAS_nondetect,
                                                   OPE_nondetect))

# Summarize percentage of non-detects by Tissue/Species/Sex/Location ----------------------
sum_4_conta_data_copy %>%
  group_by(Collection.Location) %>%
  summarise(count_non_detects = sum(is.na(Total_OPE))/130*100)
  


# non_detect_barplot <- ggplot(data = nondetect_df,
#                              aes(x = Contaminant,
#                                  y = Nondetects_quantity)) + 
#   geom_bar(aes(fill = Contaminant),
#            stat = "identity")
# 
# ggsave(paste0(getwd(), "/non_detect_barplot.png"), 
#        non_detect_barplot,
#        dpi = 480)

# Mann-Whitney test -------------------------------------------------------

# wilcox_test_tissue <- function(nrow, ncol, ) {
#   # body
#   # output
# }

# Tissue ~ Metals -------------------------------------------------------
tissue_metals_matrix <- matrix(data = NA, nrow = 7, ncol = 7)
colnames(tissue_metals_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")
rownames(tissue_metals_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")


for (i in rownames(tissue_metals_matrix)) {
  for (j in colnames(tissue_metals_matrix)) {
    subset_df_metals_1 <- filter(long_df_metals, Tissue == i)
    subset_df_metals_2 <- filter(long_df_metals, Tissue == j)
    tissue_metals_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                             y = subset_df_metals_2$Concentration)$p.value
  }
}

# Tissue ~ PBDE -----------------------------------------------------------
tissue_PBDE_matrix <- matrix(data = NA, nrow = 7, ncol = 7)
colnames(tissue_PBDE_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")
rownames(tissue_PBDE_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")


for (i in rownames(tissue_PBDE_matrix)) {
  for (j in colnames(tissue_PBDE_matrix)) {
    subset_df_metals_1 <- filter(long_df_PBDE, Tissue == i)
    subset_df_metals_2 <- filter(long_df_PBDE, Tissue == j)
    tissue_PBDE_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                             y = subset_df_metals_2$Concentration)$p.value
  }
}

# Tissue ~ PFAS -----------------------------------------------------------
tissue_PFAS_matrix <- matrix(data = NA, nrow = 7, ncol = 7)
colnames(tissue_PFAS_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")
rownames(tissue_PFAS_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")


for (i in rownames(tissue_PFAS_matrix)) {
  for (j in colnames(tissue_PFAS_matrix)) {
    subset_df_metals_1 <- filter(long_df_PFAS, Tissue == i)
    subset_df_metals_2 <- filter(long_df_PFAS, Tissue == j)
    tissue_PFAS_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                             y = subset_df_metals_2$Concentration)$p.value
  }
}

# Tissue ~ OPE ------------------------------------------------------------
tissue_OPE_matrix <- matrix(data = NA, nrow = 7, ncol = 7)
colnames(tissue_OPE_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")
rownames(tissue_OPE_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")


for (i in rownames(tissue_OPE_matrix)) {
  for (j in colnames(tissue_OPE_matrix)) {
    subset_df_metals_1 <- filter(long_df_OPE, Tissue == i)
    subset_df_metals_2 <- filter(long_df_OPE, Tissue == j)
    tissue_OPE_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                             y = subset_df_metals_2$Concentration)$p.value
  }
}


# Species ~ Metals --------------------------------------------------------
species_metals_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(species_metals_matrix) <- c("BLKI", "NOFU")
rownames(species_metals_matrix) <- c("BLKI", "NOFU")

for (i in rownames(species_metals_matrix)) {
  for (j in colnames(species_metals_matrix)) {
    subset_df_1 <- filter(long_df_metals, species == i)
    subset_df_2 <- filter(long_df_metals, species == j)
    species_metals_matrix[i,j] <- wilcox.test(x = subset_df_1$Concentration,
                                             y = subset_df_2$Concentration)$p.value
  }
}

# Species ~ PBDE ----------------------------------------------------------
species_PBDE_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(species_PBDE_matrix) <- c("BLKI", "NOFU")
rownames(species_PBDE_matrix) <- c("BLKI", "NOFU")

for (i in rownames(species_PBDE_matrix)) {
  for (j in colnames(species_PBDE_matrix)) {
    subset_df_1 <- filter(long_df_PBDE, species == i)
    subset_df_2 <- filter(long_df_PBDE, species == j)
    species_PBDE_matrix[i,j] <- wilcox.test(x = subset_df_1$Concentration,
                                              y = subset_df_2$Concentration)$p.value
  }
}

# Species ~ PFAS ----------------------------------------------------------
species_PFAS_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(species_PFAS_matrix) <- c("BLKI", "NOFU")
rownames(species_PFAS_matrix) <- c("BLKI", "NOFU")

for (i in rownames(species_PFAS_matrix)) {
  for (j in colnames(species_PFAS_matrix)) {
    subset_df_1 <- filter(long_df_PFAS, species == i)
    subset_df_2 <- filter(long_df_PFAS, species == j)
    species_PFAS_matrix[i,j] <- wilcox.test(x = subset_df_1$Concentration,
                                              y = subset_df_2$Concentration)$p.value
  }
}

# Species ~ OPE -----------------------------------------------------------
species_OPE_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(species_OPE_matrix) <- c("BLKI", "NOFU")
rownames(species_OPE_matrix) <- c("BLKI", "NOFU")

for (i in rownames(species_OPE_matrix)) {
  for (j in colnames(species_OPE_matrix)) {
    subset_df_1 <- filter(long_df_OPE, species == i)
    subset_df_2 <- filter(long_df_OPE, species == j)
    species_OPE_matrix[i,j] <- wilcox.test(x = subset_df_1$Concentration,
                                              y = subset_df_2$Concentration)$p.value
  }
}


# Sex ~ Metals ------------------------------------------------------------
sex_metals_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(sex_metals_matrix) <- c("Male", "Female")
rownames(sex_metals_matrix) <- c("Male", "Female")

for (i in rownames(sex_metals_matrix)) {
  for (j in colnames(sex_metals_matrix)) {
    subset_df_metals_1 <- filter(long_df_metals, Sex == i)
    subset_df_metals_2 <- filter(long_df_metals, Sex == j)
    sex_metals_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                           y = subset_df_metals_2$Concentration)$p.value
  }
}


# Sex ~ PBDE --------------------------------------------------------------
sex_PBDE_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(sex_PBDE_matrix) <- c("Male", "Female")
rownames(sex_PBDE_matrix) <- c("Male", "Female")

for (i in rownames(sex_PBDE_matrix)) {
  for (j in colnames(sex_PBDE_matrix)) {
    subset_df_metals_1 <- filter(long_df_PBDE, Sex == i)
    subset_df_metals_2 <- filter(long_df_PBDE, Sex == j)
    sex_PBDE_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                          y = subset_df_metals_2$Concentration)$p.value
  }
}

# Sex ~ PFAS --------------------------------------------------------------
sex_PFAS_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(sex_PFAS_matrix) <- c("Male", "Female")
rownames(sex_PFAS_matrix) <- c("Male", "Female")

for (i in rownames(sex_PFAS_matrix)) {
  for (j in colnames(sex_PFAS_matrix)) {
    subset_df_metals_1 <- filter(long_df_PFAS, Sex == i)
    subset_df_metals_2 <- filter(long_df_PFAS, Sex == j)
    sex_PFAS_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                          y = subset_df_metals_2$Concentration)$p.value
  }
}

# Sex ~ OPE ---------------------------------------------------------------
sex_OPE_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(sex_OPE_matrix) <- c("Male", "Female")
rownames(sex_OPE_matrix) <- c("Male", "Female")

for (i in rownames(sex_OPE_matrix)) {
  for (j in colnames(sex_OPE_matrix)) {
    subset_df_metals_1 <- filter(long_df_OPE, Sex == i)
    subset_df_metals_2 <- filter(long_df_OPE, Sex == j)
    sex_OPE_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                          y = subset_df_metals_2$Concentration)$p.value
  }
}


# Location ~ Metals -------------------------------------------------------
location_metals_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(location_metals_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")
rownames(location_metals_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")

for (i in rownames(location_metals_matrix)) {
  for (j in colnames(location_metals_matrix)) {
    subset_df_metals_1 <- filter(long_df_metals, Collection.Location == i)
    subset_df_metals_2 <- filter(long_df_metals, Collection.Location == j)
    location_metals_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                       y = subset_df_metals_2$Concentration)$p.value
  }
}


# Location ~ PBDE ---------------------------------------------------------
location_PBDE_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(location_PBDE_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")
rownames(location_PBDE_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")

for (i in rownames(location_PBDE_matrix)) {
  for (j in colnames(location_PBDE_matrix)) {
    subset_df_metals_1 <- filter(long_df_PBDE, Collection.Location == i)
    subset_df_metals_2 <- filter(long_df_PBDE, Collection.Location == j)
    location_PBDE_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                               y = subset_df_metals_2$Concentration)$p.value
  }
}

# Location ~ PFAS ---------------------------------------------------------
location_PFAS_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(location_PFAS_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")
rownames(location_PFAS_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")

for (i in rownames(location_PFAS_matrix)) {
  for (j in colnames(location_PFAS_matrix)) {
    subset_df_metals_1 <- filter(long_df_PFAS, Collection.Location == i)
    subset_df_metals_2 <- filter(long_df_PFAS, Collection.Location == j)
    location_PFAS_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                               y = subset_df_metals_2$Concentration)$p.value
  }
}

# Location ~ OPE ----------------------------------------------------------
location_OPE_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
colnames(location_OPE_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")
rownames(location_OPE_matrix) <- c("Labrador Sea", "Prince Leopold Island, NU")

for (i in rownames(location_OPE_matrix)) {
  for (j in colnames(location_OPE_matrix)) {
    subset_df_metals_1 <- filter(long_df_OPE, Collection.Location == i)
    subset_df_metals_2 <- filter(long_df_OPE, Collection.Location == j)
    location_OPE_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                               y = subset_df_metals_2$Concentration)$p.value
  }
}



