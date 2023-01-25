#combined analysis of ivip and MPXV infection in HeLa
#msglm run by Sabri
#Author: Yiqi Huang
###############################################################################

project_id <- 'mpxv'
message('Project ID=', project_id)
analysis_version = "20230105"

library(tidyverse)
library(ggnewscale)
library(eulerr)
library(RColorBrewer)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))

combined.df <- read_tsv(file.path(analysis_path, "reports", "downstream_analysis", "mpxv_HeLa_combined_contrasts_report_wide.txt"))

#summary plots for host and viral proteins----
combined_stats <- combined.df %>% 
  filter(ci_target == "average",
         !is_viral_ivip & !is_viral_mpxv & !is_contaminant) %>% 
  select(object_id, contains("hit"), contains("change")) %>% 
  mutate(is_hit_consistent.MVA_dE9L_6h = is_hit.MVA_dE9L_6h_VS_mock_6h & is_hit.MVA_dE9L_12h_VS_mock_12h & (change.MVA_dE9L_6h_VS_mock_6h == change.MVA_dE9L_12h_VS_mock_12h),
         is_hit_consistent.MVA_dE9L_24h = is_hit.MVA_dE9L_24h_VS_mock_24h & is_hit.MVA_dE9L_12h_VS_mock_12h & (change.MVA_dE9L_24h_VS_mock_24h == change.MVA_dE9L_12h_VS_mock_12h),
         is_hit_consistent.MVA_dE9L_12h = is_hit_consistent.MVA_dE9L_6h|is_hit_consistent.MVA_dE9L_24h,
         
         is_hit_consistent.MVA_F_6h = is_hit.MVA_F_6h_VS_mock_6h & is_hit.MVA_F_12h_VS_mock_12h & (change.MVA_F_6h_VS_mock_6h == change.MVA_F_12h_VS_mock_12h),
         is_hit_consistent.MVA_F_24h =  is_hit.MVA_F_24h_VS_mock_24h & is_hit.MVA_F_12h_VS_mock_12h & (change.MVA_F_24h_VS_mock_24h == change.MVA_F_12h_VS_mock_12h),
         is_hit_consistent.MVA_F_12h = is_hit_consistent.MVA_F_6h|is_hit_consistent.MVA_F_24h,
         
         is_hit_consistent.CVA_152_6h = is_hit.CVA_152_6h_VS_mock_6h & is_hit.CVA_152_12h_VS_mock_12h & (change.CVA_152_6h_VS_mock_6h == change.CVA_152_12h_VS_mock_12h),
         is_hit_consistent.CVA_152_24h = is_hit.CVA_152_24h_VS_mock_24h & is_hit.CVA_152_12h_VS_mock_12h & (change.CVA_152_24h_VS_mock_24h == change.CVA_152_12h_VS_mock_12h),
         is_hit_consistent.CVA_152_12h = is_hit_consistent.CVA_152_6h|is_hit_consistent.CVA_152_24h,
         
         is_hit_consistent.VACV_WR_6h = is_hit.VACV_WR_6h_VS_mock_6h & is_hit.VACV_WR_12h_VS_mock_12h & (change.VACV_WR_6h_VS_mock_6h == change.VACV_WR_12h_VS_mock_12h),
         is_hit_consistent.VACV_WR_24h = is_hit.VACV_WR_24h_VS_mock_24h & is_hit.VACV_WR_12h_VS_mock_12h & (change.VACV_WR_24h_VS_mock_24h == change.VACV_WR_12h_VS_mock_12h),
         is_hit_consistent.VACV_WR_12h = is_hit_consistent.VACV_WR_6h|is_hit_consistent.VACV_WR_24h,
         
         is_hit_consistent.MPXV_6h = is_hit.MPXV_6h_VS_mock_6h & is_hit.MPXV_12h_VS_mock_12h & (change.MPXV_6h_VS_mock_6h == change.MPXV_12h_VS_mock_12h),
         is_hit_consistent.MPXV_24h = is_hit.MPXV_24h_VS_mock_24h & is_hit.MPXV_12h_VS_mock_12h & (change.MPXV_24h_VS_mock_24h == change.MPXV_12h_VS_mock_12h),
         is_hit_consistent.MPXV_12h = is_hit_consistent.MPXV_6h|is_hit_consistent.MPXV_24h) %>% 
  rename_with(~ str_remove(., "_VS_.*"), everything()) %>% 
  pivot_longer(-object_id, 
               names_to = c( ".value","Var"), 
               names_sep="\\." ) %>% 
  filter(!str_detect(Var, "0h")) %>%
  mutate(Var = str_remove(Var, "h")) %>% 
  separate(Var, into =c("virus", "other", "time")) %>% 
  mutate(virus = ifelse(other == "dE9L", paste0(virus,"_dE3L"), virus),
         time = ifelse(is.na(time), other, time))

combined_stats_summary <- combined_stats %>% 
  filter(change == "+"|change == "-") %>% 
  group_by(virus, time, change) %>% 
  summarise(sum_hit = sum(is_hit),
            sum_consistent = sum(is_hit_consistent)) %>% 
  mutate(percentage = sum_consistent/sum_hit,
         plot_sum_hit = ifelse(change == "+", sum_hit, -sum_hit),
         plot_sum_consistent = ifelse(change == "+", sum_consistent, -sum_consistent),
         time = factor(time, levels = c("6", "12", "24")),
         virus = factor(virus, levels = c("MPXV", "VACV", "CVA", "MVA", "MVA_dE3L")))

combined_stats_summary2 <- combined_stats %>% 
  filter(change == "+"|change == "-") %>% 
  group_by(virus, time) %>% 
  summarise(sum_hit = sum(is_hit),
            sum_consistent = sum(is_hit_consistent)) %>% 
  mutate(percentage = sum_consistent/sum_hit,
         time = factor(time, levels = c("6", "12", "24")),
         virus = factor(virus, levels = c("MPXV", "VACV", "CVA", "MVA", "MVA_dE3L")))

combined_stats_virus <- combined.df %>% 
  filter(ci_target == "average",
         is_viral_ivip | is_viral_mpxv,
         !is_contaminant) %>% 
  select(object_id, is_viral_ivip, is_viral_mpxv, contains("hit")) %>% 
  rename_with(~ str_remove(., "_VS_.*"), everything()) %>% 
  pivot_longer(contains("hit"), 
               names_to = c( ".value","Var"), 
               names_sep="\\." ) %>% 
  filter(!str_detect(Var, "0h")) %>%
  mutate(Var = str_remove(Var, "h")) %>% 
  separate(Var, into =c("virus", "other", "time")) %>% 
  mutate(virus = ifelse(other == "dE9L", paste0(virus,"_dE3L"), virus),
         time = ifelse(is.na(time), other, time),
         is_hit_viral = ifelse(virus == "MPXV", is_hit&is_viral_mpxv, is_hit&is_viral_ivip))

combined_stats_virus_summary <- combined_stats_virus %>% 
  group_by(virus, time) %>% 
  summarise(sum_viral_hit = sum(is_hit_viral)) %>% 
  mutate(time = factor(time, levels = c("6", "12", "24")),
         virus = factor(virus, levels = c("MPXV", "VACV", "CVA", "MVA", "MVA_dE3L")))

 

virus_palette <- c("MPXV" ="#F9CB40", "VACV" = "#f74a0f", "CVA" = "#83e000", "MVA" ="#2DA8FF", "MVA_dE3L" = "#2D3FFF")
virus_c_palette <- c("MPXV" ="#e5ae07", "VACV" = "#b43206", "CVA" = "#569400", "MVA" ="#0083e0", "MVA_dE3L" = "#0013e0")
brewer.pal(5, "Set2")
virus_palette2 <- c("MPXV" ="#FC8D62", "VACV" = "#E78AC3", "CVA" = "#8DA0CB", "MVA" ="#66C2A5", "MVA_dE3L" = "#A6D854")

(p1 <- ggplot(combined_stats_summary)+
  geom_bar(aes(x = time, y = plot_sum_hit, fill = virus), stat = "identity")+
  scale_fill_manual(values = virus_palette)+
  geom_text(aes(x = time, y = plot_sum_hit, label = sum_hit, vjust = ifelse(plot_sum_hit >0,-0.25, 1.25)), size = 4)+
  new_scale_fill() +
  
  geom_bar(aes(x = time, y = plot_sum_consistent, fill = virus), stat = "identity")+
  scale_fill_manual(values = virus_c_palette)+
  
  geom_line(data= combined_stats_virus_summary, aes(x = time, y = sum_viral_hit, group = 1), colour = "dark grey", size = 1.5)+
    geom_text(data= combined_stats_virus_summary, aes(x = time, y = sum_viral_hit, group = 1, label = sum_viral_hit), vjust = -0.25,
              colour = "dark grey", size = 4)+
  
  facet_wrap(~ virus, nrow = 1)+
  geom_hline(yintercept = 0)+
  theme_void())

ggsave(p1, filename = file.path(analysis_path, "plots", paste0("ivip_summary_plot_", analysis_version, ".pdf")),
       width = 8, height = 4)  


(p2 <- ggplot(combined_stats_summary)+
    geom_bar(aes(x = time, y = plot_sum_hit, fill = virus), stat = "identity")+
    scale_fill_manual(values = virus_palette2)+
    #scale_fill_brewer(palette = "Set2")+
    geom_text(aes(x = time, y = plot_sum_hit, label = sum_hit, vjust = ifelse(plot_sum_hit >0,-0.25, 1.25)), size = 4)+

    facet_wrap(~ virus, nrow = 1)+
    geom_hline(yintercept = 0)+
    theme_classic())

ggsave(p2, filename = file.path(analysis_path, "plots", paste0("ivip_summary_plot_host_", analysis_version, ".pdf")),
       width = 8, height = 4) 

(p3 <- ggplot(combined_stats_virus_summary)+
    geom_bar(aes(x = time, y = sum_viral_hit, fill = virus), stat = "identity")+
    scale_fill_manual(values = virus_palette)+
    geom_text(aes(x = time, y = sum_viral_hit, label = sum_viral_hit), vjust = -0.25, size = 4)+
    
    facet_wrap(~ virus, nrow = 1)+
    #geom_hline(yintercept = 0)+
    scale_y_continuous(expand = c(0, 0))+
    theme_classic())

ggsave(p3, filename = file.path(analysis_path, "plots", paste0("ivip_summary_plot_viral_", analysis_version, ".pdf")),
       width = 8, height = 4) 

#euler plot for three poxviruses----
ivip_24h_euler_list <- list(
  VACV = combined.df %>% filter(ci_target == "average", !is_viral_ivip, ! is_viral_mpxv, is_hit.VACV_WR_24h_VS_mock_24h)#, change.VACV_WR_24h_VS_mock_24h == "-") 
  %>% select(object_id) %>% unlist,
  CVA = combined.df %>% filter(ci_target == "average", !is_viral_ivip, ! is_viral_mpxv, is_hit.CVA_152_24h_VS_mock_24h)#, change.CVA_152_24h_VS_mock_24h == "-") 
  %>% select(object_id) %>% unlist,
  MVA = combined.df %>% filter(ci_target == "average", !is_viral_ivip, ! is_viral_mpxv, is_hit.MVA_F_24h_VS_mock_24h)#, change.MVA_F_24h_VS_mock_24h == "-") 
  %>% select(object_id) %>% unlist
)
ivip_24h_euler <- euler(ivip_24h_euler_list)
(ivip_24h_euler_plot <- plot(ivip_24h_euler, fills = c("#E78AC3", "#8DA0CB", "#66C2A5"), quantities = T))

ggsave(filename = file.path(analysis_path, "plots", paste0("ivip_all_euler_host_", analysis_version, ".pdf")),
       plot = ivip_24h_euler_plot, width=6, height=8, device=cairo_pdf, family="Arial")


View(combined.df %>% filter(ci_target == "average", !is_viral_ivip, ! is_viral_mpxv, 
                            is_hit.VACV_WR_24h_VS_mock_24h, is_hit.CVA_152_24h_VS_mock_24h, is_hit.MVA_F_24h_VS_mock_24h,
                            change.MVA_F_24h_VS_mock_24h == change.VACV_WR_24h_VS_mock_24h,
                            change.VACV_WR_24h_VS_mock_24h == change.CVA_152_24h_VS_mock_24h))
