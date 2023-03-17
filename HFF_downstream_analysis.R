#downstream analysis of HFF multiomics datasets
#Author: Yiqi Huang
###############################################################################

project_id <- 'mpxv'
message('Project ID=', project_id)
analysis_version = "20230127"

library(tidyverse)
library(viridis)
library(RColorBrewer)
library(pheatmap)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))

tx.df <- read_tsv(file.path(analysis_path, "reports", "sup_tables", "Supplementary table X Transcriptome of HFF cells infected with MPXV.txt"))
fp.df <- read_tsv(file.path(analysis_path, "reports", "sup_tables","Supplementary table X Total proteome of HFF cells infected with MPXV.txt"))
pp.df <- read_tsv(file.path(analysis_path, "reports", "sup_tables", "Supplementary table X Phosphoproteome of HFF cells infected with MPXV.txt"))

#summary bar plots for each omics----
tx_summary.df <- tx.df %>% 
  select(gene_name, change.6h, change.12h, change.24h) %>% 
  pivot_longer(-gene_name, names_to = "time", values_to = "change") %>% 
  filter(change != ".", gene_name != "MPXV") %>% 
  mutate(time = as.numeric(str_remove_all(time, "change\\.|h"))) %>% 
  group_by(time, change) %>% 
  summarise(n = n()) %>% 
  mutate(plot_n = ifelse(change == "+", n, -n))

nrow(tx.df %>% filter(is_hit.6h|is_hit.12h|is_hit.24h))

fp_summary.df <- fp.df %>% 
  filter(!is_viral, !is_contaminant) %>% 
  select(object_id, change.6h, change.12h, change.24h) %>% 
  pivot_longer(-object_id, names_to = "time", values_to = "change") %>% 
  filter(change != ".") %>% 
  mutate(time = as.numeric(str_remove_all(time, "change\\.|h"))) %>% 
  group_by(time, change) %>% 
  summarise(n = n()) %>% 
  mutate(plot_n = ifelse(change == "+", n, -n))

nrow(fp.df %>% filter( !is_contaminant) %>% filter(is_hit.6h|is_hit.12h|is_hit.24h))

fp_summary_viral.df <- fp.df %>% 
  filter(is_viral, !is_contaminant) %>% 
  select(object_id, change.6h, change.12h, change.24h) %>% 
  pivot_longer(-object_id, names_to = "time", values_to = "change") %>% 
  filter(change != ".") %>% 
  mutate(time = as.numeric(str_remove_all(time, "change\\.|h"))) %>% 
  group_by(time, change) %>% 
  summarise(n = n()) %>% 
  mutate(plot_n = ifelse(change == "+", n, -n))

pp_summary.df <- pp.df %>% 
  filter(!is_viral, !is_contaminant) %>% 
  select(ptmngroup_id, change.6h, change.12h, change.24h) %>% 
  pivot_longer(-ptmngroup_id, names_to = "time", values_to = "change") %>% 
  filter(change != ".") %>% 
  mutate(time = as.numeric(str_remove_all(time, "change\\.|h"))) %>% 
  group_by(time, change) %>% 
  summarise(n = n()) %>% 
  mutate(plot_n = ifelse(change == "+", n, -n))

nrow(pp.df %>% filter( !is_contaminant) %>% filter(is_hit.6h|is_hit.12h|is_hit.24h))

pp_summary_viral.df <- pp.df %>% 
  filter(is_viral, !is_contaminant) %>% 
  select(ptmngroup_id, change.6h, change.12h, change.24h) %>% 
  pivot_longer(-ptmngroup_id, names_to = "time", values_to = "change") %>% 
  filter(change != ".") %>% 
  mutate(time = as.numeric(str_remove_all(time, "change\\.|h"))) %>% 
  group_by(time, change) %>% 
  summarise(n = n()) %>% 
  mutate(plot_n = ifelse(change == "+", n, -n))

(tx <- ggplot(tx_summary.df, aes(x = factor(time), y = plot_n))+
    geom_bar( fill = "#F9CB40", stat = "identity")+
    geom_text(aes(label = n, vjust = ifelse(plot_n >0,-0.25, 1.25)), size = 4)+
    geom_hline(yintercept = 0)+
    theme_classic())

ggsave(tx, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("hff_summary_plot_tx_", analysis_version, ".pdf")),
       width = 4, height = 4)

(fp_viral <- ggplot(fp_summary_viral.df, aes(x = factor(time), y = plot_n))+
    geom_bar( fill = "#F9CB40", stat = "identity")+
    geom_text(aes(label = n, vjust = ifelse(plot_n >0,-0.25, 1.25)), size = 4)+
    #geom_hline(yintercept = 0)+
    scale_y_continuous(expand = c(0, 0))+
    theme_classic())

ggsave(fp_viral, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("hff_summary_plot_fp_viral_", analysis_version, ".pdf")),
       width = 4, height = 4)

(pp <- ggplot(pp_summary.df, aes(x = factor(time), y = plot_n))+
    geom_bar( fill = "#F9CB40", stat = "identity")+
    geom_text(aes(label = n, vjust = ifelse(plot_n >0,-0.25, 1.25)), size = 4)+
    geom_hline(yintercept = 0)+
    #scale_y_continuous(expand = c(0, 0))+
    theme_classic())

ggsave(pp, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("hff_summary_plot_pp_", analysis_version, ".pdf")),
       width = 4, height = 4)

#supplementary table of the enrichment analysis----
enrichment.df <- read_tsv(file.path(analysis_path, "reports", "mpxv_united_hit_covers_20221117.txt"))

enrichment4paper.df <- enrichment.df %>% 
  filter(str_detect(coll_id, "GO|Reactome")) %>% 
  mutate(dataset = factor(dataset, levels = c("rnaseq", "fp", "phospho")),
         intersect_genes = str_remove_all(intersect_genes, fixed("..."))) %>% 
  select(dataset, time = timepoint_lhs, change, term_collection, term_id, term_name, term_description = term_descr,
         set_overlap_log10pvalue, set_overlap_pvalue, n_overlap = nmasked, n_dataset = nhit, n_non_overlapping_term = nunmasked, genes_overlap = intersect_genes) %>% 
  arrange(term_collection, term_id, dataset, change, time)

write_tsv(enrichment4paper.df,
          file.path(analysis_path, "reports", "sup_tables", paste0("Supplementary table X Biological functions and pathways enriched in HFF cells infected by MPXV.txt")))

phospho_enrichment.df <- read_tsv(file.path(analysis_path, "reports", "mpxv_phospho_20221105_hit_covers_20221117.txt"))

pp_enrichment4paper.df <- phospho_enrichment.df %>% 
  mutate(intersect_genes = str_remove_all(intersect_genes, fixed("..."))) %>% 
  select(time = timepoint_lhs, change, term_collection, term_name,
         set_overlap_log10pvalue, set_overlap_pvalue, n_overlap = nmasked, n_dataset = nhit, n_non_overlapping_term = nunmasked, sites_overlap = intersect_genes) %>% 
  arrange(term_name, change, time)

write_tsv(pp_enrichment4paper.df,
          file.path(analysis_path, "reports", "sup_tables", paste0("Supplementary table X_2 Biological functions and pathways enriched in HFF cells infected by MPXV.txt")))

#plot the kinase enrichment results----
uniprot_gene2id <- read_tsv(file = "/pool/analysis/yhuang/mpxv/data/rnaseq_20221116/uniprot_gene2id_20220531.tab") %>% 
  dplyr::rename(gene_name = `Gene names`) %>% 
  separate_rows(gene_name, sep = " ") %>% 
  mutate(gene_name = str_remove_all(gene_name, ";")) %>% 
  group_by(Entry) %>% 
  mutate(rank = 1:n())

kinase_viral.df <- read_tsv(file.path(analysis_path, "reports", "downstream_analysis", "all_viral_vs_all_psites_correct.txt"))

kinase_viral_hits.df <- kinase_viral.df %>% 
  filter(adjusted_p_value <= 0.01, enrichment_value_log2 >0) %>% 
  mutate(kinase = str_replace(kinase, "STLK3", "STK39")) %>% 
  left_join(select(uniprot_gene2id, kinase = gene_name, Entry)) %>% 
  left_join(select(uniprot_gene2id, gene_name, Entry, rank)) %>% 
  filter(rank == 1, 
         gene_name != "NR2C2") %>% #note that not all the kinase names are standard Uniprot gene names... NR2C2 is not a kinase
  mutate(kinase_group = factor(kinase_group, levels = c("STE", "TKL", "CK1", "Other")))

(kinase_viral <- ggplot(kinase_viral_hits.df, aes(x = enrichment_value_log2, y = reorder(gene_name, enrichment_value_log2)))+
  geom_bar(aes(fill = adjusted_p_value_log10_abs), stat = "identity")+
    scale_fill_distiller(palette = "YlOrRd", direction = 1)+
    #scale_fill_viridis(option = "inferno", limits = c(0,6), direction = -1)+
    #scale_size_continuous(range = c(1,8))+
  geom_text(aes(label = gene_name), hjust = -0.25, size = 3)+
  labs(
    x='log2(enrichment value)', y=NULL,
    #color='kinase group',
    fill='-log10(adj.p-value)'
  ) +
    expand_limits(x = (2.5))+
  theme_classic() +
    facet_grid(kinase_group ~ .,  scales = "free_y", space = "free")+
    theme(axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="bottom",
          #panel.background = element_rect(fill = "grey95", color = NA)
  )
  )


ggsave(kinase_viral, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("kinase_enrichment_plot_viral_bar_correct_", analysis_version, ".pdf")),
       width = 4, height = 6)

#as a comparison, plot the CoV2 viral psites----
kinase_viral_CoV2.df <- read_tsv(file.path(analysis_path, "reports", "downstream_analysis", "all_viral_vs_all_psites_CoV2.txt"))

kinase_viral_CoV2_hits.df <- kinase_viral_CoV2.df %>% 
  filter(p_value <= 0.05, enrichment_value_log2 >0) #%>% 
  mutate(kinase = str_replace(kinase, "STLK3", "STK39")) %>% 
  left_join(select(uniprot_gene2id, kinase = gene_name, Entry)) %>% 
  left_join(select(uniprot_gene2id, gene_name, Entry, rank)) %>% 
  filter(rank == 1, 
         gene_name != "NR2C2") %>% #note that not all the kinase names are standard Uniprot gene names... NR2C2 is not a kinase
  mutate(kinase_group = factor(kinase_group, levels = c("STE", "TKL", "CK1", "Other")))

(kinase_viral <- ggplot(kinase_viral_hits.df, aes(x = enrichment_value_log2, y = reorder(gene_name, enrichment_value_log2)))+
    geom_point(aes(size = adjusted_p_value_log10_abs))+
    scale_size_continuous(range = c(1,8))+
    geom_text(aes(label = gene_name), hjust = -0.5, size = 3)+
    labs(
      x='log2(enrichment value)', y=NULL,
      color='kinase group',size='-log10(adj.p-value)'
    ) +
    expand_limits(x = (2.5))+
    theme_classic() +
    facet_grid(kinase_group ~ .,  scales = "free_y")+
    theme(axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="bottom",
          panel.background = element_rect(fill = "grey95",
                                          color = NA))
)


ggsave(kinase_viral, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("kinase_enrichment_plot_viral_", analysis_version, ".pdf")),
       width = 4, height = 6)

#plot the kinase motif prediction enrichment for host----
kinase_host_6h.df <- read_tsv(file.path(analysis_path, "reports", "downstream_analysis", "host_6h_vs_all_psites.txt"))
kinase_host_12h.df <- read_tsv(file.path(analysis_path, "reports", "downstream_analysis", "host_12h_vs_all_psites.txt"))
kinase_host_24h.df <- read_tsv(file.path(analysis_path, "reports", "downstream_analysis", "host_24h_vs_all_psites.txt"))

kinase_host_6h_sig.df <- kinase_host_6h.df %>% 
  mutate(up_6h = ifelse(upregulated_adjusted_p_value <= 0.01, upregulated_enrichment_value_log2, 0),
         down_6h = ifelse(downregulated_adjusted_p_value <= 0.01, downregulated_enrichment_value_log2, 0)) %>% 
  select(kinase, kinase_group, up_6h, down_6h) #%>% 
  #mutate(down_6h = -down_6h)

kinase_host_12h_sig.df <- kinase_host_12h.df %>% 
  mutate(up_12h = ifelse(upregulated_adjusted_p_value <= 0.01, upregulated_enrichment_value_log2, 0),
         down_12h = ifelse(downregulated_adjusted_p_value <= 0.01, downregulated_enrichment_value_log2, 0)) %>% 
  select(kinase, kinase_group, up_12h, down_12h) #%>% 
  #mutate(down_12h = -down_12h)
  
kinase_host_24h_sig.df <- kinase_host_24h.df %>% 
  mutate(up_24h = ifelse(upregulated_adjusted_p_value <= 0.01, upregulated_enrichment_value_log2, 0),
         down_24h = ifelse(downregulated_adjusted_p_value <= 0.01, downregulated_enrichment_value_log2, 0)) %>% 
  select(kinase, kinase_group, up_24h, down_24h) #%>% 
  #mutate(down_24h = -down_24h)

kinase_host_wide_dominant <- kinase_host_6h.df %>% 
  select(kinase, kinase_group, direction_6h = dominant_direction, enrichment_value_log2_6h = dominant_enrichment_value_log2, adj_p_6h = dominant_adjusted_p_value) %>% 
  left_join(select(kinase_host_12h.df, kinase, kinase_group, direction_12h = dominant_direction, enrichment_value_log2_12h = dominant_enrichment_value_log2, adj_p_12h = dominant_adjusted_p_value)) %>% 
  left_join(select(kinase_host_24h.df, kinase, kinase_group, direction_24h = dominant_direction, enrichment_value_log2_24h = dominant_enrichment_value_log2, adj_p_24h = dominant_adjusted_p_value))

kinase_host_wide_dominant_sig <- kinase_host_wide_dominant %>% 
  filter(adj_p_6h <= 0.01|adj_p_12h <= 0.01|adj_p_24h <= 0.01) %>% 
  mutate(kinase = str_replace(kinase, "STLK3", "STK39"),
         kinase = str_replace(kinase, "SBK", "SBK1"),
         kinase = str_replace(kinase, "PKAC", "PRKAC"),
         kinase = str_replace(kinase, "PKG", "PRKG"),
         kinase = str_replace(kinase, "PKC", "PRKC"),
         kinase = str_replace(kinase, "P38A", "MAPK14"),
         kinase = str_replace(kinase, "DNAPK", "PRKDC"),
         kinase = str_replace(kinase, "P70S6KB", "RPS6KB2"),
         ) %>% 
  left_join(select(uniprot_gene2id, kinase = gene_name, Entry)) %>% 
  left_join(select(uniprot_gene2id, gene_name, Entry, rank)) %>% 
  filter(rank == 1, 
         #gene_name != "NR2C2",
         gene_name != "PRKCA",
         gene_name != "ANTXR1",
         gene_name != "SERPINA2",
         gene_name != "PAK6",
         gene_name != "ALPK3") %>% #note that not all the kinase names are standard Uniprot gene names... and there are overlapping names too
  mutate(kinase_group = forcats::fct_relevel(kinase_group, "Other", after = Inf),
         general_direction = rowSums(. == "upregulated set")) %>% 
  arrange(desc(general_direction), kinase_group, gene_name)

kinase_host_wide_only_sig <- kinase_host_6h_sig.df %>% 
  left_join(kinase_host_12h_sig.df) %>% 
  left_join(kinase_host_24h_sig.df) %>% 
  mutate(kinase_old = kinase,
    kinase = str_replace(kinase, "STLK3", "STK39"),
         kinase = str_replace(kinase, "SBK", "SBK1"),
         kinase = str_replace(kinase, "PKAC", "PRKAC"),
         kinase = str_replace(kinase, "PKG", "PRKG"),
         kinase = str_replace(kinase, "PKC", "PRKC"),
         kinase = str_replace(kinase, "P38A", "MAPK14"),
         kinase = str_replace(kinase, "DNAPK", "PRKDC"),
         kinase = str_replace(kinase, "P70S6KB", "RPS6KB2"),
         kinase_group = forcats::fct_relevel(kinase_group, "Other", after = Inf),
         general_direction = rowSums(.[3:8])
  ) %>% 
  filter(general_direction != 0) %>% 
  left_join(select(kinase_host_wide_dominant_sig, kinase, gene_name)) %>% 
  filter(!is.na(gene_name)) %>% 
  arrange(desc(general_direction), kinase_group, gene_name) %>% 
  select(gene_name, kinase, kinase_old,kinase_group, down_6h, down_12h, down_24h, up_6h, up_12h, up_24h)

kinase_host_wide_sig <- kinase_host_wide_only_sig %>% 
  select(gene_name, kinase_old, kinase_group) %>% 
  left_join(select(kinase_host_6h.df, kinase_old = kinase, up_6h = upregulated_enrichment_value_log2, down_6h = downregulated_enrichment_value_log2)) %>% 
  left_join(select(kinase_host_12h.df, kinase_old = kinase, up_12h = upregulated_enrichment_value_log2, down_12h = downregulated_enrichment_value_log2)) %>% 
  left_join(select(kinase_host_24h.df, kinase_old = kinase, up_24h = upregulated_enrichment_value_log2, down_24h = downregulated_enrichment_value_log2)) %>% 
  select(gene_name, kinase_old,kinase_group, down_6h, down_12h, down_24h, up_6h, up_12h, up_24h)


kinase_only_sig.mtx <- as.matrix(select(kinase_host_wide_only_sig, -c(gene_name:kinase_group)))
rownames(kinase_only_sig.mtx) <- kinase_host_wide_only_sig$gene_name
row_anno <- data.frame("kinase_group" = kinase_host_wide_only_sig$kinase_group)
rownames(row_anno) <-  rownames(kinase_only_sig.mtx)

paletteLength <- 256
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(-2.5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(2.5/paletteLength, 2.5, length.out=floor(paletteLength/2)))

pheatmap(kinase_only_sig.mtx,
         cluster_cols = F, cluster_rows = T,
         color = myColor,
         border_color = "white",
         breaks = myBreaks,
         cellwidth = 20, cellheight = 10,
         annotation_row = row_anno,
         legend_breaks = c(-2, -1, 0, 1, 2,  max(kinase_only_sig.mtx)), 
         main = "", legend_labels = c("-2", "-1", "0", "1", "2", "log2(enrichment)\n"),
         legend = TRUE,
         width = 6, height = 15,
         filename = file.path(analysis_path, "plots", "hff_downstream", paste0("kinase_enrichment_plot_host_only_sig_", analysis_version, ".pdf")))

kinase_sig.mtx <- as.matrix(select(kinase_host_wide_sig, -c(gene_name:kinase_group)))
rownames(kinase_sig.mtx) <- kinase_host_wide_sig$gene_name
pheatmap(kinase_sig.mtx,
         cluster_cols = F, cluster_rows = T,
         color = myColor,
         border_color = "white",
         breaks = myBreaks,
         cellwidth = 20, cellheight = 10,
         annotation_row = row_anno,
         legend_breaks = c(-2, -1, 0, 1, 2,  max(kinase_sig.mtx)), 
         main = "", legend_labels = c("-2", "-1", "0", "1", "2", "log2(enrichment)\n"),
         legend = TRUE,
         width = 6, height = 15,
         filename = file.path(analysis_path, "plots", "hff_downstream", paste0("kinase_enrichment_plot_host_sig_", analysis_version, ".pdf")))

#load all individual results for viral sites----
site_anno.df <- read_tsv("/pool/pub3rdparty/PhosphoSitePlus/20230203/Kinase_Substrate_Dataset.gz", col_names = T, skip = 2)
viral_sites.df <- read_tsv(file.path(analysis_path, "reports", "sup_tables", "Supplementary table X_2 Phosphoproteome of HFF cells infected with MPXV.txt")) %>% 
  mutate(short_seq = str_extract(flanking_15AAs, ".{6}\\*.{6}"), 
         short_seq = str_remove_all(short_seq, "\\_"))

site_anno_human.df <- site_anno.df %>% 
  filter(KIN_ORGANISM == "human") %>% 
  select(gene_name = GENE, kinase = KINASE, protein_ac = KIN_ACC_ID) %>% 
  unique() %>% 
  mutate(kinase = ifelse(kinase != "AlphaK3", toupper(kinase), kinase))

viral_site_kinase_all.df <- list.files(path= file.path(analysis_path, "reports", "downstream_analysis", "viral_site_scores"), full.names = TRUE) %>% 
  lapply(read_tsv) %>% 
  bind_rows() 
write_tsv(viral_site_kinase_all.df,
          file.path(analysis_path, "reports", "downstream_analysis", paste0("viral_kinase_sites_all.txt")))
viral_site_kinase_all.df <- read_tsv(file.path(analysis_path, "reports", "downstream_analysis", paste0("viral_kinase_sites_all.txt")))

pre_viral_site_kinase_sig.df <- viral_site_kinase_all.df %>% 
  filter(score_log2 > 0, site_percentile > 95) %>% 
  left_join(site_anno_human.df)

na_anno.df <- pre_viral_site_kinase_sig.df %>% 
  filter(is.na(protein_ac)) %>% 
  select(kinase, kinase_group) %>% 
  unique() %>% 
  mutate(old_kinase = kinase,
         kinase = str_replace_all(kinase, c(ALK2="ACVR1", smMLCK = "MYLK", caMLCK = "MYLK3", SBK = "SBK1", skMLCK = "MYLK2", PKACG = "PRKACG"))) %>%  
  left_join(select(uniprot_gene2id, kinase = gene_name, Entry)) %>% 
  left_join(select(uniprot_gene2id, gene_name, Entry, rank)) %>% 
  filter(rank == 1, gene_name != "NEK9")

viral_site_kinases_sig.df <- pre_viral_site_kinase_sig.df %>% 
  full_join(select(na_anno.df, kinase = old_kinase, protein_ac = Entry, gene_name), by = "kinase") %>% 
  mutate(gene_name = coalesce(gene_name.x, gene_name.y),
         protein_ac = coalesce(protein_ac.x, protein_ac.y)) 

viral_sites_kinases.df <- viral_sites.df %>% 
  left_join(select(viral_site_kinases_sig.df, kinase, kinase_group, kinase_gene_name = gene_name, kinase_ac = protein_ac, short_seq = sequence )) %>% 
  mutate(VACV_homo = str_extract(protein_description, "(?<=\\()Cop.*(?=\\))"),
         VACV_homo = str_remove(VACV_homo, " "))
write_tsv(viral_sites_kinases.df,
          file.path(analysis_path, "reports", "downstream_analysis", paste0("viral_kinase_sites_sig_annotated.txt")))

viral_sites_timeXkinases.df <- viral_sites_kinases.df %>% 
  filter(ptm_status != "potential") %>% 
  group_by(gene_name, protein_description, time, short_seq) %>% 
  summarise(count = sum(!is.na(kinase))) %>% 
  group_by(gene_name, protein_description, short_seq) %>% 
  slice_head() %>% 
  replace_na(list(time = 0))

(timeXsites <- ggplot(viral_sites_timeXkinases.df, aes(factor(time), count))+
  geom_boxplot(aes(group = time))+
  geom_jitter())

ggsave(timeXsites, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("timeXsigmotifs_plot_viral_", analysis_version, ".pdf")),
       width = 4, height = 4)

viral_proteins_kinase_stat.df <- viral_sites_kinases.df %>% 
  group_by(gene_name, VACV_homo, short_seq, kinase_gene_name) %>% 
  arrange(time) %>% 
  slice_head() %>% 
  group_by(gene_name, VACV_homo) %>% 
  #ungroup() %>% 
  summarise(count = sum(!is.na(kinase))) %>% 
  mutate(gene_name = str_remove_all(gene_name, "MPXV_|L$|R$"), 
         VACV_homo = str_remove(VACV_homo, "L$|R$"))

(stat <- ggplot(viral_proteins_kinase_stat.df, aes(x = count, y = reorder(paste0(gene_name, " (", VACV_homo, ")"), count))) +
  geom_bar(stat = "identity")+
  labs(y = NULL)+
  theme(axis.ticks.y = element_blank())+
    theme_classic()+
    scale_x_continuous(expand = c(0, 0)))
  
ggsave(stat, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("proteinXmotif_plot_viral_", analysis_version, ".pdf")),
       width = 4, height = 6)

viral_sites_kinase_stat.df <- viral_sites_kinases.df %>% 
  group_by(gene_name, VACV_homo, short_seq, kinase_gene_name) %>% 
  arrange(time) %>% 
  slice_head() %>% 
  group_by(gene_name, VACV_homo, ptm_pos, ptm_AA, short_seq) %>% 
  summarise(count = sum(!is.na(kinase))) %>% 
  mutate(gene_name = str_remove_all(gene_name, "MPXV_|L$|R$"), 
         VACV_homo = str_remove(VACV_homo, "L$|R$"))

#present the MAPKs and their sites as a network----
MAPK_sites.df <- viral_sites_kinases.df %>% 
  filter(str_detect(kinase_gene_name, "MAP|RAF")) %>% 
  mutate(MAPK_group = str_extract(kinase_gene_name, "(?<=MAP)."),
         MAPK_group = str_replace(MAPK_group, "K", "1")) %>% 
  replace_na(list(MAPK_group = "3")) %>% #BRAF and RAF are MAP3Ks
  mutate(MAPK_group = as.numeric(MAPK_group))

MAPK_sites.df %>% filter(ptm_status != "potential") %>% group_by(gene_name, ptm_pos, ptm_AA) %>% arrange(time) %>% slice_head() %>% View()

MAPK_sites.df %>% group_by(gene_name, ptm_pos, ptm_AA, kinase_gene_name) %>% arrange(time) %>% slice_head() %>% group_by(MAPK_group) %>% summarise(n())

MAPK_sites_by_viral_protein.df <- MAPK_sites.df %>% 
  filter(ptm_status != "potential") %>%
  group_by(gene_name, ptm_pos, ptm_AA, kinase_gene_name) %>% 
  arrange(time) %>% 
  slice_head() %>% 
  ungroup() %>% 
  left_join(select(viral_site_kinases_sig.df, kinase_gene_name = gene_name, short_seq = sequence, score_log2, site_percentile)) %>% 
  mutate(gene_name = str_remove_all(gene_name, "MPXV_|L$|R$"), 
         VACV_homo = str_remove(VACV_homo, "L$|R$"))
write_tsv(MAPK_sites_by_viral_protein.df,
          file.path(analysis_path, "reports", "downstream_analysis", paste0("viral_MAPK_sites_sig_annotated.txt")))

sources <- MAPK_sites_by_viral_protein.df %>% 
  distinct(kinase_gene_name) %>% 
  rename(label = kinase_gene_name)

destinations <- MAPK_sites_by_viral_protein.df %>%
  mutate(gene_name = paste0(gene_name, " (", VACV_homo,")")) %>% 
  distinct(gene_name) %>% 
  rename(label = gene_name)

nodes <- full_join(sources, destinations, by = "label") %>% 
  mutate(id = 1:nrow(.)) %>%
  select(id, everything())

edges <- MAPK_sites_by_viral_protein.df %>% 
  mutate(gene_name = paste0(gene_name, " (", VACV_homo,")")) %>% 
  rename(weight = score_log2) %>% 
  left_join(nodes, by = c("kinase_gene_name" = "label")) %>% 
  rename(from = id) %>% 
  left_join(nodes, by = c("gene_name" = "label")) %>% 
  rename(to = id) %>% 
  select(from, to, weight)

library(tidygraph)
library(ggraph)

net.tidy <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
ggraph(net.tidy, layout = "grid") + 
  geom_node_point() +
  geom_edge_link(aes(width = weight), alpha = 0.8) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_text(aes(label = label), repel = TRUE, max.overlaps = Inf) +
  labs(edge_width = "log2(enrichment)") +
  theme_graph(base_family="sans")

#plot a dotplot for the MAPKs----
MAPK_dotplot.df <- MAPK_sites_by_viral_protein.df %>% 
  select(gene_name, kinase_gene_name, MAPK_group, score_log2, site_percentile) 

MAPK_wide.df <- MAPK_dotplot.df %>% 
  select(-c(site_percentile, MAPK_group)) %>% 
  group_by(gene_name, kinase_gene_name) %>% 
  slice_max(score_log2) %>% 
  ungroup() %>% 
  pivot_wider(names_from = gene_name, values_from = score_log2) %>% 
  column_to_rownames(var="kinase_gene_name")

MAPK_wide.df[is.na(MAPK_wide.df)] <- 0
d <- dist(scale(MAPK_wide.df), method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1, cex = 0.6, hang = -1)
hc1$labels[hc1$order]

d2 <- dist(scale(t(MAPK_wide.df)), method = "euclidean")
hc2 <- hclust(d2, method = "complete" )
plot(hc2, cex = 0.6, hang = -1)

MAPK_dotpot_ordered.df <- MAPK_dotplot.df %>% 
  mutate(kinase_gene_name = factor(kinase_gene_name, levels = hc1$labels[hc1$order]),
         gene_name = factor(gene_name, levels = hc2$labels[hc2$order]))
         
(dot_plot <- ggplot(MAPK_dotpot_ordered.df, aes(x = gene_name, y = kinase_gene_name, size = score_log2, colour = site_percentile))+
  geom_point()+
  scale_colour_distiller(palette = "YlOrRd", direction = 1)+
  theme_classic()+
  theme(axis.line  = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x= "", y = ""))

(MAPK_label <- ggplot(MAPK_dotpot_ordered.df, 
       aes(x = 1, y = kinase_gene_name, fill = factor(MAPK_group, levels = c("4", "3", "2", "1")))) + 
  geom_tile() + 
  scale_fill_brewer(palette = 'Pastel2')+
  theme_minimal())
(total_MAPK <- ggplot(MAPK_dotpot_ordered.df, aes(x = gene_name))+
    geom_bar()+
    theme_classic()
    )

ggsave(dot_plot, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("MAPK_viral_sites_", analysis_version, ".pdf")),
       width = 8, height = 8)

ggsave(MAPK_label, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("MAPK_viral_sites_label_", analysis_version, ".pdf")),
       width = 5, height = 8)

ggsave(total_MAPK, filename = file.path(analysis_path, "plots", "hff_downstream", paste0("MAPK_viral_sites_count_", analysis_version, ".pdf")),
       width = 8, height = 5)

#tidy up the drug enrichment hits from Valter----
drug_targets_wide.df <- read_tsv(file.path(analysis_path, "reports", "downstream_analysis", paste0("drug_targets_from_Valter_20230306b.txt")))

drug_targets_long.df <- drug_targets_wide.df %>% 
  rename(gene_name = node_1_name) %>% 
  mutate(protein_id = 1:nrow(.)) %>% 
  left_join(uniprot_gene2id %>% select(gene_name, protein_ac = Entry, rank)) %>% 
  group_by(gene_name) %>% 
  arrange(rank) %>% 
  slice_head() %>% 
  ungroup() %>% 
  rename(original_gene_name = gene_name) %>% 
  filter(!is.na(protein_ac)) %>% 
  left_join(uniprot_gene2id %>% filter(rank == 1) %>% select(protein_ac = Entry, gene_name)) %>%
  select(-original_gene_name, -rank) %>% 
  pivot_longer(-c(protein_id, gene_name, protein_ac), names_to = "comparison", values_to = "is_hit")

write_tsv(drug_targets_long.df, file.path(analysis_path, "reports", "downstream_analysis", paste0("drug_targets_from_Valter_long_20230306b.txt")))
  
