# VirScan_VRC_codes
# Author: Andy Chan
# Date: 7Jan2026

# Libraries ---------------------------------------------------------------------
library(DESeq2)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(viridis)
library(readr)
library(RColorBrewer)
library(gridExtra)
library(gghighlight)
library(ggh4x)
library(readxl)
library(ggpubr)
library(writexl)
library(openxlsx)


# All transformed datasets (VARscore, HFC, TPH)--------------------
VARscore_1_WTMeV_LAMV <- readRDS("VARscore_WTMeV_LAMV_latest_latest.RDS")%>%
  select(-Duplicate)
Merged_hfc_all_WTMeV_LAMV <- readRDS("Merged_hfc_WTMeV_LAMV_latest.RDS")
Total_peptide_hits_all_WTMeV_LAMV<- readRDS("Total_peptide_hits_WTMeV_LAMV_latest.RDS")

VARscore_1_VRC <- readRDS("VARscore_VRC_latest.RDS") %>%
  select(-Duplicate)
Merged_hfc_all_VRC <-readRDS("Merged_hfc_VRC_latest.RDS")
Total_peptide_hits_all_VRC<- readRDS("Total_peptide_hits_VRC_latest.RDS")

# Total_number_of_peptide_tiles_per_pathogen_all_table
Total_number_of_peptide_tiles_per_pathogen_all <- Merged_hfc_all_WTMeV_LAMV %>%
  group_by(taxon_species, pep_aa, pos_start, pos_end, pep_id)%>%
  summarize()

Total_number_of_peptide_tiles_per_pathogen_all_table <- as.data.frame(table(Total_number_of_peptide_tiles_per_pathogen_all$taxon_species))%>%
  filter(!row_number() %in% c(1))%>%
  dplyr:: rename(taxon_species= Var1, 
         Number_of_Peptide_Tiles = Freq)

# Varscore for all viruses
VARscore_1_VRC_vacc_group_DPV_all_viruses <- VARscore_1_VRC %>%
  filter (Treatment %in% c("2rMeVH", "2MMR", "Boost")) %>%
  filter (DPI %in% c("DPV0", "DPV14", "DPV28", "DPV84", "DPV169"))%>%
  mutate(DPI = recode(DPI,
                      "DPV0" = "d0",
                      "DPV14" = "d14_15",
                      "DPV28" = "d28",
                      "DPV84" = "d84_89",
                      "DPV169" = "d112_168_176"))%>%
  mutate (DPI = factor (DPI, levels = c("d0", "d14_15", "d28", "d84_89", "d112_168_176")))%>%
  select(-Duplicate)

VARscore_1_WTMeV_LAMV_all_viruses <- VARscore_1_WTMeV_LAMV %>%
  mutate (Treatment = recode(Treatment,
                             "LAMV" = "1LAMV"))%>%
  filter (Treatment %in% c("WTMeV", "1LAMV"))%>%
  select(-Duplicate)

VARscore_1_Sham_control_all_viruses <- VARscore_1_WTMeV_LAMV %>%
  filter (Treatment %in% c("Sham_control"))%>%
  select(-Duplicate)
  
# Varscore DPV 
VARscore_1_WTMeV_LAMV_selected <- VARscore_1_WTMeV_LAMV %>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))

VARscore_1_VRC_vacc_group_DPV <- VARscore_1_VRC %>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))%>% 
  filter (Treatment %in% c("2rMeVH", "2MMR", "Boost")) %>%
  filter (DPI %in% c("DPV0", "DPV14", "DPV28", "DPV84", "DPV169"))%>%
  mutate(DPI = recode(DPI,
                              "DPV0" = "d0",
                              "DPV14" = "d14_15",
                              "DPV28" = "d28",
                              "DPV84" = "d84_89",
                              "DPV169" = "d112_168_176"))%>%
  mutate (DPI = factor (DPI, levels = c("d0", "d14_15", "d28", "d84_89", "d112_168_176")))

VARscore_1_DPV<- rbind (VARscore_1_WTMeV_LAMV_selected, VARscore_1_VRC_vacc_group_DPV) %>%
  mutate (Treatment = recode(Treatment,
                                  "LAMV" = "1LAMV",
                                  "Boost" = "1MMR/1rMeVH")) %>%
  mutate (Treatment = factor (Treatment, levels = c("WTMeV", "1LAMV", "2MMR", "2rMeVH", "1MMR/1rMeVH", "Sham_control")))

# Varscore DPC 
VARscore_1_WTMeV_selected <- VARscore_1_WTMeV_LAMV %>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))%>%
  filter (Treatment %in% c("WTMeV")) %>% 
  filter (DPI %in% c("d0", "d14_15", "d28", "d84_89"))
  
VARscore_1_VRC_vacc_group_DPC <- VARscore_1_VRC %>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))%>% 
  filter (Treatment %in% c("2rMeVH", "2MMR", "Boost")) %>%
  filter (DPI %in% c("DPC0", "DPC14", "DPC28", "DPC84"))%>%
  mutate(DPI = recode(DPI,
                      "DPC0" = "d0",
                      "DPC14" = "d14_15",
                      "DPC28" = "d28",
                      "DPC84" = "d84_89"))%>%
  mutate (DPI = factor (DPI, levels = c("d0", "d14_15", "d28", "d84_89"))) %>%
  mutate (Treatment = factor (Treatment, levels = c("2MMR", "2rMeVH", "Boost")))

VARscore_1_DPC<- rbind (VARscore_1_WTMeV_selected, VARscore_1_VRC_vacc_group_DPC) %>%
  mutate (Treatment = recode(Treatment,
                             "LAMV" = "1LAMV",
                             "Boost" = "1MMR/1rMeVH")) %>%
  mutate (Treatment = factor (Treatment, levels = c("WTMeV", "1LAMV", "2MMR", "2rMeVH", "1MMR/1rMeVH", "Sham_control")))

# Total_peptide_hits_all_DPV
Total_peptide_hits_all_WTMeV_LAMV_selected <- Total_peptide_hits_all_WTMeV_LAMV%>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))

Total_peptide_hits_all_VRC_vacc_group_DPV <- Total_peptide_hits_all_VRC %>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))%>% 
  filter (Treatment %in% c("2rMeVH", "2MMR", "Boost"))%>%
  filter (DPI %in% c("DPV0", "DPV14", "DPV28", "DPV84", "DPV169"))%>%
  mutate(DPI = recode(DPI,
                      "DPV0" = "d0",
                      "DPV14" = "d14_15",
                      "DPV28" = "d28",
                      "DPV84" = "d84_89",
                      "DPV169" = "d112_168_176"))%>%
  mutate (DPI = factor (DPI, levels = c("d0", "d14_15", "d28", "d84_89", "d112_168_176")))

Total_peptide_hits_all_DPV<- rbind (Total_peptide_hits_all_WTMeV_LAMV_selected, Total_peptide_hits_all_VRC_vacc_group_DPV) %>%
  mutate (Treatment = recode(Treatment,
                             "LAMV" = "1LAMV",
                             "Boost" = "1MMR/1rMeVH")) %>%
  mutate (Treatment = factor (Treatment, levels = c("WTMeV", "1LAMV", "2MMR", "2rMeVH", "1MMR/1rMeVH", "Sham_control")))

# Total_peptide_hits_all_DPC
Total_peptide_hits_all_WTMeV_selected <- Total_peptide_hits_all_WTMeV_LAMV%>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus")) %>%
  filter (Treatment %in% c("WTMeV")) %>% 
  filter (DPI %in% c("d0", "d14_15", "d28", "d84_89"))

Total_peptide_hits_all_VRC_vacc_group_DPC <- Total_peptide_hits_all_VRC %>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))%>% 
  filter (Treatment %in% c("2rMeVH", "2MMR", "Boost")) %>%
  filter (DPI %in% c("DPC0", "DPC14", "DPC28", "DPC84"))%>%
  mutate(DPI = recode(DPI,
                      "DPC0" = "d0",
                      "DPC14" = "d14_15",
                      "DPC28" = "d28",
                      "DPC84" = "d84_89"))%>%
  mutate (DPI = factor (DPI, levels = c("d0", "d14_15", "d28", "d84_89"))) %>%
  mutate (Treatment = factor (Treatment, levels = c("2MMR", "2rMeVH", "Boost")))

Total_peptide_hits_all_DPC<- rbind (Total_peptide_hits_all_WTMeV_selected, Total_peptide_hits_all_VRC_vacc_group_DPC) %>%
  mutate (Treatment = recode(Treatment,
                             "LAMV" = "1LAMV",
                             "Boost" = "1MMR/1rMeVH")) %>%
  mutate (Treatment = factor (Treatment, levels = c("WTMeV", "1LAMV", "2MMR", "2rMeVH", "1MMR/1rMeVH", "Sham_control")))

# PPH DPV
Percentage_peptide_hits_all_DPV <- Total_peptide_hits_all_DPV %>%
  inner_join(Total_number_of_peptide_tiles_per_pathogen_all_table, by = "taxon_species") %>%
  mutate(pph = total_peptide_hits/Number_of_Peptide_Tiles*100) %>%
  mutate (Treatment = recode(Treatment,
                             "LAMV" = "1LAMV",
                             "Boost" = "1MMR/1rMeVH")) %>%
  mutate (Treatment = factor (Treatment, levels = c("WTMeV", "1LAMV", "2MMR", "2rMeVH", "1MMR/1rMeVH", "Sham_control")))

# PPH DPC
Percentage_peptide_hits_all_DPC <- Total_peptide_hits_all_DPC %>%
  inner_join(Total_number_of_peptide_tiles_per_pathogen_all_table, by = "taxon_species") %>%
  mutate(pph = total_peptide_hits/Number_of_Peptide_Tiles*100) %>%
    mutate (Treatment = recode(Treatment,
                                  "LAMV" = "1LAMV",
                                  "Boost" = "1MMR/1rMeVH")) %>%
  mutate (Treatment = factor (Treatment, levels = c("WTMeV", "1LAMV", "2MMR", "2rMeVH", "1MMR/1rMeVH", "Sham_control")))

# Merged_hfc_all_DPV
Merged_hfc_all_WTMeV_LAMV_latest <- Merged_hfc_all_WTMeV_LAMV %>%
  mutate(DPI = recode(DPI,
                    "DPV0" = "d0",
                    "DPV14" = "d14_15",
                    "DPV28" = "d28",
                    "DPV84" = "d84_89",
                    "DPV169" = "d112_168_176"))%>%
  mutate (DPI = factor (DPI, levels = c("d0", "d14_15", "d28", "d84_89", "d112_168_176")))%>%
  mutate (Treatment = recode(Treatment,
                             "LAMV" = "1LAMV"))

Merged_hfc_all_WTMeV_LAMV_selected <- Merged_hfc_all_WTMeV_LAMV%>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))

Merged_hfc_all_VRC_vacc_group_DPV <- Merged_hfc_all_VRC %>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))%>% 
  filter (Treatment %in% c("2rMeVH", "2MMR", "Boost"))%>%
  filter (DPI %in% c("DPV0", "DPV14", "DPV28", "DPV84", "DPV169"))%>%
  mutate(DPI = recode(DPI,
                      "DPV0" = "d0",
                      "DPV14" = "d14_15",
                      "DPV28" = "d28",
                      "DPV84" = "d84_89",
                      "DPV169" = "d112_168_176"))%>%
  mutate (DPI = factor (DPI, levels = c("d0", "d14_15", "d28", "d84_89", "d112_168_176")))

Merged_hfc_all_DPV<- rbind (Merged_hfc_all_WTMeV_LAMV_selected, Merged_hfc_all_VRC_vacc_group_DPV) %>%
  mutate (Treatment = recode(Treatment,
                             "LAMV" = "1LAMV",
                             "Boost" = "1MMR/1rMeVH")) %>%
  mutate (Treatment = factor (Treatment, levels = c("WTMeV", "1LAMV", "2MMR", "2rMeVH", "1MMR/1rMeVH", "Sham_control")))

Merged_hfc_all_DPV_rubella <- Merged_hfc_all_DPV %>% # Rubella gene_symbol missing; need to add them back
  filter(taxon_species == "Rubella virus") %>%
  mutate(gene_symbol = UniProt_acc)%>%
  mutate_at ("gene_symbol", str_replace_all, "^A4KAD6$", "E1") %>%  
  mutate_at ("gene_symbol", str_replace_all, "^A8INB2$", "E1") %>%
  mutate_at ("gene_symbol", str_replace_all, "^A8WCE6$", "E1") %>%
  mutate_at ("gene_symbol", str_replace_all, "^A8WCF3$", "E1") %>%
  mutate_at ("gene_symbol", str_replace_all, "^D5MBS8$", "E1") %>%
  mutate_at ("gene_symbol", str_replace_all, "^D5KJ87$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^P07566$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^P08563$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^P19725$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q107X9$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q6X2U1$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q6X2U3$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q9J6K8$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q8VA10$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q99IE6$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^O40955$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^P13889$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q6X2U2$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q6X2U4$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q86500$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q8BCR0$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q99IE5$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q99IE7$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q9J6K9$", "NSP")

# Merged_hfc_all_DPC
Merged_hfc_all_WTMeV_selected <- Merged_hfc_all_WTMeV_LAMV%>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))  %>%
  filter (Treatment %in% c("WTMeV")) %>% 
  filter (DPI %in% c("d0", "d14_15", "d28", "d84_89"))

Merged_hfc_all_VRC_vacc_group_DPC <- Merged_hfc_all_VRC %>%
  filter (taxon_species %in% c("Measles virus", "Mumps virus", "Rubella virus"))%>% 
  filter (Treatment %in% c("2rMeVH", "2MMR", "Boost")) %>%
  filter (DPI %in% c("DPC0", "DPC14", "DPC28", "DPC84"))%>%
  mutate(DPI = recode(DPI,
                      "DPC0" = "d0",
                      "DPC14" = "d14_15",
                      "DPC28" = "d28",
                      "DPC84" = "d84_89"))%>%
  mutate (DPI = factor (DPI, levels = c("d0", "d14_15", "d28", "d84_89"))) %>%
  mutate (Treatment = factor (Treatment, levels = c("2MMR", "2rMeVH", "Boost")))

Merged_hfc_all_DPC<- rbind (Merged_hfc_all_WTMeV_selected, Merged_hfc_all_VRC_vacc_group_DPC) %>%
  mutate (Treatment = recode(Treatment,
                             "LAMV" = "1LAMV",
                             "Boost" = "1MMR/1rMeVH")) %>%
  mutate (Treatment = factor (Treatment, levels = c("WTMeV", "1LAMV", "2MMR", "2rMeVH", "1MMR/1rMeVH", "Sham_control")))

Merged_hfc_all_DPC_rubella <- Merged_hfc_all_DPC %>% # Rubella gene_symbol missing; need to add them back
  filter(taxon_species == "Rubella virus") %>%
  mutate(gene_symbol = UniProt_acc)%>%
  mutate_at ("gene_symbol", str_replace_all, "^A4KAD6$", "E1") %>%  
  mutate_at ("gene_symbol", str_replace_all, "^A8INB2$", "E1") %>%
  mutate_at ("gene_symbol", str_replace_all, "^A8WCE6$", "E1") %>%
  mutate_at ("gene_symbol", str_replace_all, "^A8WCF3$", "E1") %>%
  mutate_at ("gene_symbol", str_replace_all, "^D5MBS8$", "E1") %>%
  mutate_at ("gene_symbol", str_replace_all, "^D5KJ87$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^P07566$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^P08563$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^P19725$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q107X9$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q6X2U1$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q6X2U3$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q9J6K8$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q8VA10$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q99IE6$", "SP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^O40955$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^P13889$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q6X2U2$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q6X2U4$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q86500$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q8BCR0$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q99IE5$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q99IE7$", "NSP") %>%
  mutate_at ("gene_symbol", str_replace_all, "^Q9J6K9$", "NSP")


# Calculate VARscore and PPH threshold for three viruses ------
VARscore_1_DPV_summary <- VARscore_1_DPV %>%
  filter (Treatment == "Sham_control") %>%
  group_by (taxon_species)%>% 
  summarize (Longitudinal_VarScore_mean = mean(VAR))

Percentage_peptide_hits_all_DPV_summary <- Percentage_peptide_hits_all_DPV %>%
  filter (Treatment == "Sham_control") %>%
  group_by (taxon_species)%>% 
  summarize (Longitudinal_PPH_mean = mean(pph))


# Calculate VARscore and PPH per timepoint and group (means and SD) ------
VARscore_DPV_summary <- VARscore_1_DPV %>%
  group_by (taxon_species, Treatment, DPI)%>% 
  summarize (Longitudinal_VarScore_mean = mean(VAR), SD = sd(VAR))

VARscore_DPC_summary <- VARscore_1_DPC %>%
  group_by (taxon_species, Treatment, DPI)%>% 
  summarize (Longitudinal_VarScore_mean = mean(VAR), SD = sd(VAR))

Percentage_peptide_hits_DPV_summary <- Percentage_peptide_hits_all_DPV %>%
  group_by (taxon_species, Treatment, DPI)%>% 
  summarize (Longitudinal_PPH_mean = mean(pph), SD = sd(pph))

Percentage_peptide_hits_DPC_summary <- Percentage_peptide_hits_all_DPC %>%
  group_by (taxon_species, Treatment, DPI)%>% 
  summarize (Longitudinal_PPH_mean = mean(pph), SD = sd(pph))

# write.xlsx(VARscore_DPV_summary, "VARscore_DPV_summary.xlsx")
# write.xlsx(VARscore_DPC_summary, "VARscore_DPC_summary.xlsx")
# write.xlsx(Percentage_peptide_hits_DPV_summary, "Percentage_peptide_hits_DPV_summary.xlsx")
# write.xlsx(Percentage_peptide_hits_DPC_summary, "Percentage_peptide_hits_DPC_summary.xlsx")


# Fig1 ---------------------
VAR_PPH_graph_function_1 <- function (datatable, viral_species, y_axis_var, x_lab, y_lab, ncol_var, threshold_value) {
  ggplot(subset(datatable, taxon_species  == viral_species), aes(x = DPI, y = {{y_axis_var}}, group = Macaques_ID))+
    geom_line(aes(color = taxon_species), linewidth = 1, colour="lightgrey")+
    labs(x = x_lab,  y= y_lab) +
    facet_wrap(~Treatment, ncol = ncol_var, scales = "fixed") +
    theme_classic(base_size=10)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    stat_summary(aes(group=Treatment), fun=mean, geom="line", colour="red", linewidth = 2)+
    labs(title = paste0(y_lab, " of ", viral_species))+
    geom_hline(yintercept = threshold_value, linetype = "dotted", color = "black", linewidth = 1)
}

plot1 <- VAR_PPH_graph_function_1 (VARscore_1_DPV, "Measles virus", VAR, "Days Post Vaccination (DPV)", "VarScores", 6, -0.03)
plot2 <- VAR_PPH_graph_function_1 (subset(VARscore_1_DPV, Treatment %in% c("WTMeV", "2MMR", "1MMR/1rMeVH")), "Mumps virus", VAR, "Days Post Vaccination (DPV)", "VarScores", 6, 0.013)
plot3 <- VAR_PPH_graph_function_1 (subset(VARscore_1_DPV, Treatment %in% c("WTMeV", "2MMR", "1MMR/1rMeVH")), "Rubella virus", VAR, "Days Post Vaccination (DPV)", "VarScores", 6, 0.288)
grid.arrange(grobs = list (plot1, plot2, plot3), ncol=1)

my_plot_1A <- plot1
ggsave("Fig.1A.tiff", plot = my_plot_1A, width = 7, height = 2.5, dpi = 600, units = "in")

my_plot_1C <- grid.arrange(grobs = list (plot2, plot3), nrow=1)
ggsave("Fig.1C.tiff", plot = my_plot_1C, width = 7, height = 2.5, dpi = 600, units = "in")


VAR_PPH_graph_function_2 <- function (datatable, viral_species, y_axis_var, x_lab, y_lab, ncol_var, DPI_var){
  ggplot(subset(datatable, taxon_species == viral_species), aes(x = Treatment, y = {{y_axis_var}}, fill = Treatment)) +
    geom_boxplot()+
    geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
    scale_fill_manual(values = c("#F8766D", "#ffcc00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
    # stat_compare_means(method = "wilcox.test", comparison = list (c("WTMeV", "1LAMV"), c("WTMeV", "2MMR"), c("WTMeV", "2rMeVH"), c("WTMeV", "1MMR/1rMeVH")), size = 3, label = "p.signif", hide.ns = TRUE) +
    stat_compare_means(method = "wilcox.test", comparison = list (c("WTMeV", "1LAMV"), c("1LAMV", "2MMR"), c("1LAMV", "2rMeVH"), c("1LAMV", "1MMR/1rMeVH"), c("2MMR", "2rMeVH"), c("2MMR", "1MMR/1rMeVH"), c("2rMeVH", "1MMR/1rMeVH")), size = 3, label = "p.signif", hide.ns = TRUE) +
    theme_classic(base_size = 10) +    
    facet_wrap(~DPI, ncol = ncol_var)+
    labs(x = x_lab,
         y = y_lab,
         fill = x_lab)+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))+
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),     # Remove x-axis text (tick labels)
      axis.ticks.x = element_blank()     # Remove x-axis ticks
    )+
    labs(title = paste0(y_lab, " of ", viral_species, " - ", DPI_var))  
}
plot1 <- VAR_PPH_graph_function_2 (VARscore_1_DPV, "Measles virus", VAR, "Treatment groups", "VarScores", 5, "DPV")

my_plot_1B <- plot1
ggsave("Fig.1B.tiff", plot = my_plot_1B, width = 7, height = 2.5, dpi = 600, units = "in")

# Fig2A and Fig4A ------------
Viral_peptide_hits_datatable_function <- function(datatable, virus) {
  peptide_hits <-datatable %>%
    filter(taxon_species == virus)%>%
    filter_at(vars(gene_symbol), all_vars(!is.na(.)))%>%
    mutate(hit_or_not = ifelse(hfc > 1, 1,0))%>%
    group_by(Macaques_ID, gene_symbol, DPI, Treatment, UniProt_acc)%>%
    mutate(total_peptide_hits= sum(hit_or_not)) %>%
    group_by(gene_symbol, DPI, Treatment)%>%
    mutate(average_peptide_hits = mean(total_peptide_hits))
  return (peptide_hits)
}
Average_peptide_hits_per_protein_function <- function (datatable, virus, color_pattern, DPI_var){
  print(ggplot(datatable, aes(x=DPI , y = gene_symbol))+
          geom_raster(aes(fill = average_peptide_hits))+
          scale_fill_gradient2(low ="white", high=color_pattern, 
                               midpoint=0) +
          labs(x = DPI_var,  y= "Protein",  title = paste0("Antibody binding epitopes per " , virus, " viral protein"), fill = "Average peptide hit per protein") +
          facet_wrap( ~Treatment, ncol = 6, scales = "fixed")+
          theme_classic(base_size=10)+
          theme(legend.position="bottom")+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
}
Measles_peptide_hits_DPV <- Viral_peptide_hits_datatable_function(Merged_hfc_all_DPV, "Measles virus") %>%
  mutate (gene_symbol = if_else(gene_symbol == "P/V", "P/V/C", gene_symbol))%>% #change problematic gene_symbol naming 
  mutate (gene_symbol = if_else(gene_symbol == "P", "P/V/C", gene_symbol))

Measles_peptide_hits_DPV_test <- Measles_peptide_hits_DPV %>%
  filter (average_peptide_hits >1 & gene_symbol == "H") 
plot1 <-Average_peptide_hits_per_protein_function(Measles_peptide_hits_DPV, "Measles", "springgreen4", "Days post vaccination (DPV)")

my_plot_2A <- plot1
ggsave("Fig.2A.tiff", plot = my_plot_2A, width = 7, height = 3, dpi = 600, units = "in")


Measles_peptide_hits_DPC <- Viral_peptide_hits_datatable_function(Merged_hfc_all_DPC, "Measles virus") %>%
  mutate (gene_symbol = if_else(gene_symbol == "P/V", "P/V/C", gene_symbol))%>% #change problematic gene_symbol naming 
  mutate (gene_symbol = if_else(gene_symbol == "P", "P/V/C", gene_symbol))
plot1 <- Average_peptide_hits_per_protein_function(Measles_peptide_hits_DPC, "Measles", "springgreen4", "Days post challenge (DPC)") + theme (legend.position = "right")

my_plot_4A <- plot1
ggsave("Fig.4A.tiff", plot = my_plot_4A, width = 7, height = 2.5, dpi = 600, units = "in")

# Fig2B and Fig4B ------------------------
Viral_peptide_hits_datatable_function_2 <- function(datatable, virus) {
  peptide_hits <-datatable %>%
    filter(taxon_species == virus)%>%
    filter_at(vars(gene_symbol), all_vars(!is.na(.)))%>%
    mutate(hit_or_not = ifelse(hfc > 1, 1,0))%>%
    group_by(Macaques_ID, gene_symbol,DPI,Treatment, UniProt_acc)%>%
    mutate(total_peptide_hits= sum(hit_or_not)) %>%
    group_by(gene_symbol,DPI,Treatment)%>%
    summarize(average_peptide_hits = mean(total_peptide_hits))
  return (peptide_hits)
}
Viral_peptide_hits_datatable_function_3 <- function(datatable, virus) {
  peptide_hits <-datatable %>%
    filter(taxon_species == virus)%>%
    filter_at(vars(gene_symbol), all_vars(!is.na(.)))%>%
    mutate(hit_or_not = ifelse(hfc > 1, 1,0))%>%
    group_by(Macaques_ID, gene_symbol,DPI,Treatment, UniProt_acc)%>%
    mutate(total_peptide_hits= sum(hit_or_not)) %>%
    mutate(gene_symbol_UniProt_acc = paste0(gene_symbol, "_", UniProt_acc)) %>%
    group_by(gene_symbol, gene_symbol_UniProt_acc,DPI,Treatment)%>%
    summarize(average_peptide_hits = mean(total_peptide_hits))
  return (peptide_hits)
}

Measles_peptide_hits_DPV <- Viral_peptide_hits_datatable_function_2(Merged_hfc_all_DPV, "Measles virus") %>%
  mutate (gene_symbol = if_else(gene_symbol == "P/V", "P/V/C", gene_symbol))%>% 
  mutate (gene_symbol = if_else(gene_symbol == "P", "P/V/C", gene_symbol))

Measles_peptide_hits_DPC <- Viral_peptide_hits_datatable_function_2(Merged_hfc_all_DPC, "Measles virus") %>%
  mutate (gene_symbol = if_else(gene_symbol == "P/V", "P/V/C", gene_symbol))%>%
  mutate (gene_symbol = if_else(gene_symbol == "P", "P/V/C", gene_symbol))

Measles_peptide_hits_DPV_per_peptide_species <- Viral_peptide_hits_datatable_function_3(Merged_hfc_all_DPV, "Measles virus")
Measles_peptide_hits_DPC_per_peptide_species <- Viral_peptide_hits_datatable_function_3(Merged_hfc_all_DPC, "Measles virus") 

Average_peptide_hits_per_protein_function_selected_lineplot <- function (datatable, virus, protein_var, DPI_var){
  print(ggplot(subset(datatable, gene_symbol == protein_var), aes(x=DPI , y = average_peptide_hits))+
          geom_line(aes(group = Treatment))+
          geom_point(aes(group = Treatment))+
          labs(x = DPI_var, y= "average peptide hits",  title = paste0(virus, " ", protein_var, " protein - antibody binding")) +
          facet_wrap( ~Treatment, ncol = 6, scales = "fixed")+
          theme_classic(base_size=10)+
          theme(legend.position="bottom")+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)))
  
}
plot1 <- Average_peptide_hits_per_protein_function_selected_lineplot (Measles_peptide_hits_DPV, "Measles", "H", "Days Post Vaccination (DPV)")
plot2 <- Average_peptide_hits_per_protein_function_selected_lineplot (Measles_peptide_hits_DPC, "Measles", "H", "Days Post Challenge (DPC)")
grid.arrange(grobs = list (plot1, plot2), ncol=1)

my_plot_2B <-grid.arrange(grobs = list (plot1, plot2), ncol=1)
ggsave("Fig.2B.tiff", plot = my_plot_2B, width = 7, height = 4.5, dpi = 600, units = "in")

my_plot_4B <-grid.arrange(grobs = list (plot1, plot2), ncol=1)
ggsave("Fig.4B.tiff", plot = my_plot_4B, width = 7, height = 4.5, dpi = 600, units = "in")

# Fig2C and Fig4C ------------------------
Epitopes_list_and_dataset_function <- function (datatable, virus, hfc_value){
  Epitopes_list <- datatable %>%
    filter(taxon_species == virus) %>%
    mutate(Viral_epitope_tile = paste0 (gene_symbol, pos_start, "_", pos_end )) %>%
    filter (hfc > hfc_value) 
  
  Epitopes_list <- unique(Epitopes_list$Viral_epitope_tile)
  
  Epitopes_dataset<- datatable %>%
    filter(taxon_species == virus) %>%
    mutate(Viral_epitope_tile = paste0 (gene_symbol, pos_start, "_", pos_end ))%>%
    filter(Viral_epitope_tile %in% Epitopes_list) %>%
    group_by(Viral_epitope_tile, DPI, Treatment, UniProt_acc)%>%
    mutate(averaged_hfc = mean (hfc))
  
  return (Epitopes_dataset)
}
Epitopes_per_group_graph_function <- function (datatable, virus, DPI_var){
  
  ggplot(datatable, aes(x = DPI , y = Viral_epitope_tile))+
    geom_raster(aes(fill = averaged_hfc)) +
    scale_fill_gradientn(colors= c("white", "darkblue", "black")) +
    labs(x = DPI_var, y = "Peptide tile", title = paste0("Antibody binding epitopes per " , virus, " viral peptide tile"),fill = "Average hit-fold change" ) +
    facet_wrap(~Treatment, ncol = 6, scales = "fixed") +
    theme_classic(base_size=10) +
    theme(axis.text = element_text(size=8))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(legend.position="bottom")
  
}

Measles_epitopes_dataset <- Epitopes_list_and_dataset_function (Merged_hfc_all_DPV, "Measles virus", 20)
plot1 <- Epitopes_per_group_graph_function(Measles_epitopes_dataset, "Measles", "Days post vaccination (DPV)")

my_plot_2C <- plot1
ggsave("Fig.2C.tiff", plot = my_plot_2C, width = 7, height = 6, dpi = 600, units = "in")

Measles_epitopes_dataset <- Epitopes_list_and_dataset_function (Merged_hfc_all_DPC, "Measles virus", 20)
plot1 <- Epitopes_per_group_graph_function(Measles_epitopes_dataset, "Measles", "Days post challenge (DPC)") + theme (legend.position = "right")

my_plot_4C <- plot1
ggsave("Fig.4C.tiff", plot = my_plot_4C, width = 7, height = 5.5, dpi = 600, units = "in")

# Fig3 ---------------
plot1 <-VAR_PPH_graph_function_1 (VARscore_1_DPC, "Measles virus", VAR, "Days Post Challenge (DPC)", "VarScores", 4, -0.03)
plot2 <-VAR_PPH_graph_function_1 (subset(VARscore_1_DPC, Treatment %in% c("WTMeV", "2MMR", "1MMR/1rMeVH")), "Mumps virus", VAR, "Days Post Challenge (DPC)", "VarScores", 4, 0.013)
plot3 <-VAR_PPH_graph_function_1 (subset(VARscore_1_DPC, Treatment %in% c("WTMeV", "2MMR", "1MMR/1rMeVH")), "Rubella virus", VAR, "Days Post Challenge (DPC)", "VarScores", 4, 0.288)
grid.arrange(grobs = list (plot1, plot2, plot3), ncol=1)

my_plot_3 <- grid.arrange(grobs = list (plot1, plot2, plot3), ncol=1)
ggsave("Fig.3.tiff", plot = my_plot_3, width = 5, height = 6, dpi = 600, units = "in")

# Fig5 ----------
# Latest dataset without D0
Varscore_ELISA_titers_WTMeV <- read_excel("Varscore_ELISA_titers_WTMeV_without_d0.xlsx")

plot_1 <- ggscatter(Varscore_ELISA_titers_WTMeV, 
          x = "Varscore", 
          y = "IgG_ELISA_titer", 
          add = "reg.line", 
          conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightgray")) +
  scale_y_continuous(trans = "log10", labels = scales::scientific) + 
  labs(x = "VarScore", y = "log10 MeV-specific IgG titers") +
  ggtitle("Spearman coefficient of \nMeV VarScore and \nMeV-specific Ab \nELISA titers") +
  stat_cor(method = "spearman", 
           label.y = 5.3)+
  theme_classic(base_size = 9)

plot_2 <- ggscatter(Varscore_ELISA_titers_WTMeV, 
          x = "Varscore", 
          y = "IgG_avidity_index", 
          add = "reg.line", 
          conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightgray")) +
  labs(x = "VarScore", y = "Avidity Index") +
  ggtitle("Spearman coefficient of \nMeV VarScore and \nMeV-specific Ab \navidity index") +
  stat_cor(method = "spearman",
           label.y = 3.5)+
  theme_classic(base_size = 9)

plot_3 <- ggscatter(Varscore_ELISA_titers_WTMeV, 
          x = "Varscore", 
          y = "PRNT", 
          add = "reg.line", 
          conf.int = TRUE, 
          add.params = list(color = "blue", fill = "lightgray")) +
  scale_y_continuous(trans = "log10", labels = scales::scientific) + 
  labs(x = "VarScore", y = "PRNT") +
  ggtitle("Spearman coefficient of \nMeV VarScore and \nMeV-specific PRNT") +
  stat_cor(method = "spearman")+
  theme_classic(base_size = 9)

my_plot_5 <- grid.arrange(grobs = list (plot_1, plot_2, plot_3), nrow=1)
ggsave("Fig.5.tiff", plot = my_plot_5, width = 8, height = 4, dpi = 600, units = "in")
