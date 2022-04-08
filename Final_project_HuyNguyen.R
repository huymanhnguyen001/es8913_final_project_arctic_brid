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

sum_4_conta_data <- select(data, Tissue, species, Sex, Collection.Location, 
                           Total_PBDE, metals, Total_PFAS, Total_OPE)

## Replace zero values with a value that is 0.5*LOD 
sum_4_conta_data$Total_PBDE[sum_4_conta_data$Total_PBDE == 0] <- runif(sum(sum_4_conta_data$Total_PBDE == 0),
                                             min = 0.0145000,
                                             max=0.0155000)
sum_4_conta_data$metals[sum_4_conta_data$metals == 0] <- runif(sum(sum_4_conta_data$metals == 0),
                                                 min = 0.0145000,
                                                 max=0.0155000)
sum_4_conta_data$Total_PFAS[sum_4_conta_data$Total_PFAS == 0] <- runif(sum(sum_4_conta_data$Total_PFAS == 0),
                                             min = 0.0145000,
                                             max=0.0155000)
sum_4_conta_data$Total_OPE[sum_4_conta_data$Total_OPE == 0] <- runif(sum(sum_4_conta_data$Total_OPE == 0),
                                           min = 0.0145000,
                                           max=0.0155000)

# Scale data to reduce outlier effect
sum_4_conta_scale_data <- copy(sum_4_conta_data)
sum_4_conta_scale_data[5:8] <- scale(sum_4_conta_scale_data[5:8])


# Pivot longer --------------------------------------------------------------
long_df <- pivot_longer(sum_4_conta_data,
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

# Facet box plots Tissue -------------------------------------------------

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
  ylim(0,0.2) +
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


# Facet box plots Location---------------------------------------------------------

location_boxplot <- ggplot(data = long_df, aes(x = Collection.Location, y = Concentration)) + 
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
        legend.position = "bottom")

# Exporting plot 
ggsave(paste0(getwd(), "/location_boxplot_2.png"), 
       location_boxplot,
       dpi = 320, 
       height = 10, 
       width = 12)


# Facet box plots Sex -----------------------------------------------------

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
        legend.position = "bottom")

# Exporting plot 
ggsave(paste0(getwd(), "/sex_boxplot_2.png"), 
       sex_boxplot,
       dpi = 320,
       height = 10)

# Facet box plots Species -------------------------------------------------

spe_boxplot <- ggplot(data = long_df, aes(x = species, y = Concentration)) + 
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
        legend.position = "bottom")

# Exporting plot 
ggsave(paste0(getwd(), "/spe_boxplot_2.png"), 
       spe_boxplot,
       dpi = 320,
       height = 10)

# Scatter plot with ellipses metals, PBDE, PFAS, OPE ------------------------------------

# Ellipse by groups
metal_PFAS <- ggplot(sum_4_conta_data, aes(x = metals, 
                             y = Total_PFAS)) +
  geom_point(aes(color = Tissue))

metal_PBDE <- ggplot(sum_4_conta_data, aes(x = metals, 
                                           y = Total_PBDE)) +
  geom_point(aes(color = Tissue))

scatter_plot <- ggarrange(metal_PBDE, metal_PFAS, 
                          ncol=2, 
                          common.legend = TRUE,
                          legend = "right")

ggsave(paste0(getwd(), "/scatter_plot.png"), 
       scatter_plot,
       dpi = 320)
# Change the type of ellipses: possible values are "t", "norm", "euclid"
# p + stat_ellipse()


# Mann-Whitney test -------------------------------------------------------

wilcox_test_tissue <- function(nrow, ncol, ) {
  # body
  # output
}

# Tissue ~ Metals -------------------------------------------------------
tissue_metals_matrix <- matrix(data = NA, nrow = 7, ncol = 7)
colnames(tissue_metals_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")
rownames(tissue_metals_matrix) <- c("blood", "brain", "egg", "fat", "liver", "muscle", "preen_oil")
long_df_metals <- filter(long_df, Contaminant == "metals")

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
long_df_PBDE <- filter(long_df, Contaminant == "Total_PBDE")

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
long_df_PFAS <- filter(long_df, Contaminant == "Total_PFAS")

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
long_df_OPE <- filter(long_df, Contaminant == "Total_OPE")

for (i in rownames(tissue_OPE_matrix)) {
  for (j in colnames(tissue_OPE_matrix)) {
    subset_df_metals_1 <- filter(long_df_OPE, Tissue == i)
    subset_df_metals_2 <- filter(long_df_OPE, Tissue == j)
    tissue_OPE_matrix[i,j] <- wilcox.test(x = subset_df_metals_1$Concentration,
                                             y = subset_df_metals_2$Concentration)$p.value
  }
}

## Species

## Sex

## Location

