#Differential Expression with DESeq2 for mxpv infected HFFs
#Yiqi Huang
#2022.08.13

project_id <- "mpxv"
datatype <- "rnaseq"
data_version <- 20220812
analysis_version <- 20220813

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, "ggplot_ext.R"))

library(DESeq2)
library(Glimma)
library(tidyverse)
library(pheatmap)
library(msglm)
library(BiocParallel)
library(factoextra)
library(ggrepel)
library(ashr)
library(parallel)
library(ggrastr)
library(rlang)
library(writexl)

#Setup parallel environment
register(SnowParam(8, exportglobals = FALSE))#For parallel computing with functions from DESeq2

mop.max_nprocesses <- 8
mop.nprocesses <- 8

if ( !exists("mop.nprocesses") ) {
  # identify the number of CPU cores available for the job
  mop.nprocesses <- ( if( Sys.getenv('SLURM_CPUS_PER_TASK') != '' ) as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) # SLURM way
                      else if ( Sys.getenv('NSLOTS') != '' ) as.integer(Sys.getenv('NSLOTS')) else 1 ) * # SGE way
    ( if( Sys.getenv('SLURM_NTASKS') != '' ) as.integer(Sys.getenv('SLURM_NTASKS')) else 1 )
  message( mop.nprocesses, ' CPU(s) available' )
  if ( exists("mop.max_nprocesses") ) mop.nprocesses <- min( mop.nprocesses, mop.max_nprocesses )
}

if ( mop.nprocesses > 1 ) {
  message( 'Creating ', mop.nprocesses, '-SOCKET cluster' )#Windows doesn't support forking, use PSOCK instead
  require( parallel )
  mop.cluster <- makeCluster( mop.nprocesses, outfile = "")#The default is "PSOCK" if we remove type = "FORK"
  clusterSetRNGStream( mop.cluster, 1293143 )
} else {
  message( 'Not using parallel computation' )
  mop.cluster <- NULL
}

clusterEvalQ(mop.cluster, library(tidyverse))
clusterEvalQ(mop.cluster, library(DESeq2))
clusterEvalQ(mop.cluster, library(ashr))

rna_data_path <- file.path(data_path, str_c(datatype, "_", data_version))
data <- read_tsv(file.path(rna_data_path, "487VG_Human_DGE_Matrix_noheader.txt"), guess_max = 1000) 
annotation <- read_tsv(file.path(rna_data_path, "487VG_Human_SampleAnnotation_notail.txt"))
plot_path <- file.path(analysis_path, "plots", str_c(datatype, "_", data_version, "_", analysis_version))
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

#Data pre-processing ----
counts <- data %>% 
  pivot_longer(-GENE, names_to = "sample", values_to = "count") %>% 
  left_join(select(annotation, Barcode, UniqueSampleID), by = c("sample" = "Barcode")) %>% 
  select(-sample) %>% 
  pivot_wider(names_from = UniqueSampleID, values_from = count) %>% 
  remove_rownames %>% 
  column_to_rownames(var="GENE") 

coldata <- filter(annotation, Loaded == 1) %>% 
  column_to_rownames(var = "UniqueSampleID") %>% 
  mutate(treatment = str_extract(rownames(.), "(?<=h\\-)\\w+(?=\\-TR)"),
         treatment = str_replace(treatment, "MPX", "MPXV"),
         timepoint = str_extract(rownames(.), "(?<=\\-)\\d+(?=h)"),
         replicate = str_extract(rownames(.), "(?<=\\-TR)\\d+$"),
         group = interaction(treatment, timepoint)) %>%
  droplevels()

#Creating a design matrix and contrasts----
#condition = treatment X timepoint (adapted from Alexey's CoV2 timecourse infection analyis)
conditionsXsamples.df <- dplyr::select(coldata, condition = group, treatment, timepoint) %>% 
  mutate(timepoint_num = as.numeric(timepoint),
         timepoint = factor(timepoint),
         treatment = factor(treatment, levels = c("mock", "MPXV")),
         sample = rownames(.)) %>%  #note that for DESeq2 analysis, each sample needs a row, not only each condition
  remove_rownames() %>% 
  dplyr::distinct() %>% 
  left_join(dplyr::select(., treatment, timepoint_after=timepoint_num), by = "treatment") %>% 
  dplyr::mutate(is_after = timepoint_num >= timepoint_after) %>%
  distinct() %>% 
  tidyr::pivot_wider(all_of(c("condition", "treatment", "timepoint", "timepoint_num", "sample")),
                     names_prefix = "after", names_from = timepoint_after, values_from = is_after) %>%
  rename_at(vars(starts_with("after")), ~str_c(., "h")) %>% 
  mutate(after0h = NULL)

#The following df is just conditions like in the original script for msglm
conditions.df <- dplyr::select(coldata, condition = group, treatment, timepoint) %>% 
  mutate(timepoint_num = as.numeric(timepoint),
         timepoint = factor(timepoint),
         treatment = factor(treatment, levels = c("mock", "MPXV"))) %>%  
  dplyr::distinct() %>% 
  left_join(dplyr::select(., treatment, timepoint_after=timepoint_num), by = "treatment") %>% 
  dplyr::mutate(is_after = timepoint_num >= timepoint_after) %>%
  distinct() %>% 
  tidyr::pivot_wider(all_of(c("condition", "treatment", "timepoint", "timepoint_num")),
                     names_prefix = "after", names_from = timepoint_after, values_from = is_after) %>%
  rename_at(vars(starts_with("after")), ~str_c(., "h")) %>% 
  mutate(after0h = NULL)

all_conditions <- as.character(conditions.df$condition)

sampleXeffect_orig.mtx <- model.matrix(
  formula(str_c("~ (" , str_c(str_subset(colnames(conditionsXsamples.df), "^after\\d+h"), collapse =" + "), " ) * treatment")),
  conditionsXsamples.df
)

conditionXeffect_orig.mtx <- model.matrix(
  formula(str_c("~ (" , str_c(str_subset(colnames(conditions.df), "^after\\d+h"), collapse =" + "), " ) * treatment")),
  conditions.df
)

sampleXeffect.mtx <- sampleXeffect_orig.mtx[, str_detect(colnames(sampleXeffect_orig.mtx),
                                                         "^treatment[^:]+$", negate=TRUE)] #&
                                              #(apply(conditionXeffect_orig.mtx, 2, function(x) min(abs(x))) == 0)] 
#Remove the treatment (we don't need the general effect of the treatment)
conditionXeffect.mtx <- conditionXeffect_orig.mtx[, str_detect(colnames(conditionXeffect_orig.mtx),
                                                               "^treatment[^:]+$", negate=TRUE)] #&
                                                    #(apply(conditionXeffect_orig.mtx, 2, function(x) min(abs(x))) == 0)] 

dimnames(sampleXeffect.mtx) <- list(condition = as.character(conditionsXsamples.df$sample),
                                    colnames(sampleXeffect.mtx))
dimnames(conditionXeffect.mtx) <- list(condition = as.character(conditions.df$condition),
                                       colnames(conditionXeffect.mtx))

pheatmap(sampleXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(plot_path,
                              paste0(project_id,"_",datatype, "_exp_design_DESeq2_", analysis_version, ".pdf")),
         width = 8, height = 12)
pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(plot_path,
                              paste0(project_id,"_",datatype, "_exp_design_condition_", analysis_version, ".pdf")),
         width = 8, height = 8)

compound_metaconditions <- c()
all_metaconditions <- c(levels(conditions.df$condition), compound_metaconditions)
conditionXmetacondition.mtx <- msglm::constant_matrix(FALSE, list(condition = all_conditions,
                                                                  metacondition = all_metaconditions))

for (cname in levels(conditions.df$condition)) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
} #In this setting, each metacondition corresponds to one condition 

contrasts.df <- bind_rows(
  # all treatment pairs at each timepoint
  left_join(dplyr::mutate(dplyr::select(conditions.df, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint_num)),
            dplyr::mutate(dplyr::select(conditions.df, metacondition_rhs = condition, treatment_rhs = treatment, timepoint_lhs = timepoint_num))) %>%
    filter(as.integer(treatment_lhs) > as.integer(treatment_rhs)) %>%
    mutate(timepoint_rhs = timepoint_lhs,
           contrast_kind = "treatment_vs_treatment"),
  # all time points of the same treatment 
  left_join(dplyr::mutate(dplyr::select(conditions.df, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint_num)),
            dplyr::mutate(dplyr::select(conditions.df, metacondition_rhs = condition, treatment_lhs = treatment, timepoint_rhs = timepoint_num))) %>%
    filter(as.integer(timepoint_lhs) > as.integer(timepoint_rhs),
           treatment_lhs != "mock") %>%
    mutate(treatment_rhs = treatment_lhs,
           contrast_kind = "timepoint_vs_timepoint")
) %>%
  mutate(contrast = case_when(contrast_kind == "treatment_vs_treatment" ~ str_c(treatment_lhs, "_vs_", treatment_rhs, "@", timepoint_rhs, "h"),
                              contrast_kind == "timepoint_vs_timepoint" ~ str_c(timepoint_lhs, "h_vs_", timepoint_rhs, "h@", treatment_rhs)),
         contrast_type = "comparison") %>% 
  distinct()
all_contrasts = contrasts.df$contrast

metaconditionXcontrast.mtx <- msglm::constant_matrix(0, list(metacondition = all_metaconditions,
                                                             contrast = all_contrasts))

for (i in 1:nrow(contrasts.df)) {
  contr <- contrasts.df$contrast[[i]]
  metaconditionXcontrast.mtx[contrasts.df$metacondition_lhs[[i]], contr] <- 1.0
  metaconditionXcontrast.mtx[contrasts.df$metacondition_rhs[[i]], contr] <- -1.0
}

contrastXmetacondition.mtx <- t(metaconditionXcontrast.mtx)

pheatmap(contrastXmetacondition.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plot_path,
                              paste0(project_id,"_", datatype, "_exp_design_contrasts_", analysis_version, ".pdf")),
         width = 8, height = 16)

conditionXeffect_reorder.mtx <- conditionXeffect.mtx[colnames(contrastXmetacondition.mtx), ,drop=FALSE]

pheatmap(conditionXeffect_reorder.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         # filename = file.path(plot_path, paste0(project_id,"_",datatype, "_exp_design_reorder_", analysis_version, ".pdf")),
         width = 8, height = 12)

contrastXeffect.mtx <- contrastXmetacondition.mtx %*%  conditionXeffect_reorder.mtx
contr.matrix <- t(contrastXeffect.mtx)
pheatmap(contrastXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plot_path,
                              paste0(project_id,"_", datatype, "_exp_design_contrastsXeffect_", analysis_version, ".pdf")),
         width = 8, height = 12, fontsize_row = 5)

#Assembly the matrices for DESeq2 analysis----
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = sampleXeffect.mtx)

keep <- rowSums(counts(dds) >= 10) >= 2 #Very very mild filtering
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)#This information of normalised count is useful for visualisation, but DESeq2 itself doesn't need it.
normalized_counts.df <- as.data.frame(normalized_counts) %>%  
  mutate(GeneID = rownames(.)) %>% 
  pivot_longer(-GeneID, names_to = "sample", values_to = "normalizedcount") %>% 
  mutate(log2normalizedcount = log2(normalizedcount+0.5)) %>% #plus 0.5 to avoid -Inf values for the 0s
  left_join(conditionsXsamples.df, by = "sample") 

counts_by_condition.df <- normalized_counts.df %>%
  group_by(GeneID, condition) %>%
  summarise(count_mean = mean(normalizedcount), .groups = "drop")

#Run DESeq2 analysis!----
dds <- DESeq(dds, parallel = TRUE)
plotDispEsts(dds)

clusterExport(mop.cluster, varlist=c("dds", "conditionXeffect.mtx", "contrastXeffect.mtx"))
DESeq_conditions <- clusterApplyLB(mop.cluster, 1:nrow(conditionXeffect.mtx), 
                                   function(i) results(dds, contrast = conditionXeffect.mtx[i,], tidy = TRUE) %>% mutate(condition = rownames(conditionXeffect.mtx)[i]))
DESeq_conditions.df <- bind_rows(DESeq_conditions) %>% 
  left_join(conditions.df) %>% 
  dplyr::rename("GeneID" = "row") %>% 
  mutate("q2.5" = log2FoldChange - qt(0.975,df=2)*lfcSE,
         "q97.5" = log2FoldChange + qt(0.975,df=2)*lfcSE,
         "q25" = log2FoldChange - qt(0.75,df=2)*lfcSE,
         "q75" = log2FoldChange + qt(0.75,df=2)*lfcSE) %>% 
  select(GeneID, condition, "q2.5", "q25", "q50" = log2FoldChange, "q75", "q97.5",lfcSE, treatment, timepoint, timepoint_num )

DESeq_contrasts <- clusterApplyLB(mop.cluster, 1:nrow(contrastXeffect.mtx), 
                                  function(i) results(dds, contrast = contrastXeffect.mtx[i,], tidy = TRUE) %>% mutate(contrast = rownames(contrastXeffect.mtx)[i]))
DESeq_contrasts.df <- bind_rows(DESeq_contrasts) %>% 
  left_join(contrasts.df) %>% 
  dplyr::rename("GeneID" = "row") 

#Shrink the contrasts, no need to shrink the conditions (strange line plots)
DESeq_ashr_contrasts<- clusterApplyLB(mop.cluster, 1:nrow(contrastXeffect.mtx), 
                                      function(i) lfcShrink(dds, contrast = contrastXeffect.mtx[i,], type = "ashr") %>%
                                        as.data.frame() %>% 
                                        mutate(contrast = rownames(contrastXeffect.mtx)[i],
                                               GeneID = rownames(.)))

DESeq_ashr_contrasts.df <- bind_rows(DESeq_ashr_contrasts) %>% 
  mutate(log2MeanNormCount = log2(baseMean)) %>% 
  #replace_na(list(status = FALSE)) %>% 
  left_join(contrasts.df, by = "contrast") %>%
  left_join(counts_by_condition.df %>% dplyr::rename(count_mean_lhs = count_mean), by = c("GeneID", "metacondition_lhs" = "condition" )) %>%
  left_join(counts_by_condition.df %>% dplyr::rename(count_mean_rhs = count_mean), by = c("GeneID", "metacondition_rhs" = "condition"))

save(dds, normalized_counts, normalized_counts.df,# rld, rld.df, rld_design, rld_design.df, 
     DESeq_contrasts.df, DESeq_ashr_contrasts.df, DESeq_conditions.df, #DESeq_ashr_conditions.df, 
     conditions.df, conditionsXsamples.df, contrasts.df,
     conditionXeffect.mtx,sampleXeffect.mtx, contrastXeffect.mtx, contrastXmetacondition.mtx,
     file = file.path(results_path, str_c(project_id, '_DESeq2_data_', analysis_version, '_full.RData')))

#Saving just the raw model information for the combi analysis later
save(dds, 
     file = file.path(results_path, str_c(project_id, '_DESeq2_data_', analysis_version, '.RData')))

#Make plots ----
plot_path <- file.path(analysis_path, "plots", str_c(datatype,"_", data_version, "_", analysis_version))
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
treatment_palette <- c("mock" = "gray", "MPXV" = "#F9CB40")

#Plot modelled counts for all genes 
sel_objects.df <- DESeq_conditions.df %>% 
  filter(!str_detect(GeneID,"-AS" ),
         !str_detect(GeneID, "-\\d{4,}"),
         !str_detect(GeneID, "AC\\d{3,}"),
         !str_detect(GeneID, "LINC")) %>% 
  inner_join(filter(DESeq_contrasts.df, baseMean > 1) %>% select(GeneID) %>% unique())


group_by(sel_objects.df, GeneID) %>% 
  do({
    sel_obj.df <- .
    #sel_obj_id <- sel_obj.df$GeneID[1]
    timecourse_plots_path <- file.path(plot_path,"rnaseq_timecourse" )
    if (!dir.exists(timecourse_plots_path)) { dir.create(timecourse_plots_path, recursive=TRUE) }
    plot_filename_prefix <- sel_obj.df$GeneID[1]
    plot_title_prefix <- sel_obj.df$GeneID[1]
    message("Plotting ", plot_title_prefix)
    
    if (nrow(sel_obj.df) > 0) {
      p <-
        ggplot(data=sel_obj.df, aes(x = timepoint_num, color=treatment, fill = treatment)) +
        geom_ribbon(aes(x = timepoint_num, ymin = `q2.5`, ymax=`q97.5`),
                    alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
        geom_ribbon(aes(x = timepoint_num, ymin = `q25`, ymax=`q75`),
                    alpha=0.5, stat = "identity", size=0.5) +
        geom_line(aes(x = timepoint_num, y = `q50`), alpha = 0.5, size=1, linetype = "solid", stat="identity") +
        geom_line(data = sel_obj.df, aes(x = timepoint_num, y = `50%`), 
                  alpha=0.5,linetype = "longdash", size=1, stat="identity") +
        geom_point(data=filter(normalized_counts.df, GeneID %in% sel_obj.df$GeneID),
                   aes(x = timepoint_num, y = log2normalizedcount), position = position_jitter(width = 0.1, height = 0), size=1) +
        #geom_point(data=filter(rld.df, GeneID == sel_obj.df$object_label[1]),
        #           aes(x = timepoint_num, y = rld), position = position_jitter(width = 0.75, height = 0), size=1) +
        theme_bw_ast(base_family = "", base_size = 8) +
        scale_x_continuous(breaks=unique(conditions.df$timepoint_num)) +
        scale_color_manual(values=treatment_palette2) +
        scale_fill_manual(values=treatment_palette2) +
        ylab("log2 normalized counts")+
        xlab("timepoint/h")+
        ggtitle(str_c(plot_title_prefix, " rna count timecourse"))+
        #subtitle=sel_obj.df$protein_descriptions) +
        facet_wrap( ~ GeneID, scales = "free")
      NULL
      ggsave(plot=p, filename=file.path(timecourse_plots_path, str_c(plot_filename_prefix, "_rnaseq_timecourse", ".pdf")),
             device=cairo_pdf, width=12, height=8)
      
    }
    tibble()
  })

# plot treatment contrasts per timepoint volcano
object_rna_contrast_thresholds.df <- tibble(
  #contrast_kind = "treatment_vs_treatment",
  contrast_type = "comparison",
  p_value_threshold = 0.00001,
  #log2FC_threshold = log2(1.5),
  log2FC_threshold = 1,
  log2FC_max = 3,
  contrast_offset_log2 = 0
)

object_contrasts_4show.df <- DESeq_ashr_contrasts.df %>%
  dplyr::inner_join(object_rna_contrast_thresholds.df) %>%
  drop_na() %>% 
  dplyr::mutate(is_signif = padj <= p_value_threshold & abs(log2FoldChange) >= log2FC_threshold,
                is_hit = is_signif & (count_mean_lhs >= 50 | count_mean_rhs >= 50),
                p_value_compressed = 10^(-sapply(-log10(padj), mlog10_pvalue_compress)),
                show_label = coalesce(is_hit, FALSE),
                log2FC_trunc = pmax(-log2FC_max, pmin(log2FC_max, log2FoldChange)),
                truncation = scatter_truncation(log2FoldChange, log2FC_trunc, padj, padj, is_hit | !is_signif),
                truncation_type = point_truncation_type(truncation, is_signif),
                status = case_when(is_hit & log2FoldChange > 0 ~ 1,
                is_hit & log2FoldChange < 0 ~ -1,
                TRUE ~ 0),
                change = case_when(status == 1 ~ "+", 
                                   status == -1 ~ "-",
                                   TRUE ~ ".")) %>%
  dplyr::group_by(contrast) %>%
  dplyr::mutate(show_label = if_else(rep.int(sum(show_label) >= 400L, n()), is_hit, show_label)) %>%
  dplyr::ungroup()

contrasts_stats.df <- object_contrasts_4show.df %>%
  filter(!str_detect(GeneID,"-AS" ),
         !str_detect(GeneID, "-\\d{4,}"),
         !str_detect(GeneID, "AC\\d{3,}"),
         !str_detect(GeneID, "LINC")) %>% 
  group_by(contrast) %>% 
  summarise(sum_sig = sum(is_signif), sum_hit = sum(is_hit), 
  sum_hit_plus = sum(status == 1), sum_hit_minus = sum(status == -1))

object_contrasts_4show_test.df <- object_contrasts_4show.df %>% 
  filter(contrast == "MPXV_vs_mock@6h")

#
#object_thomas.df <- read.table(file.path(results_path, "rnaseq_Thomas_20220818", "DESeq2_Descr_Dummy_mpx24hrs_vs_mock24hrs_N5.txt"), 
#dec = ",", sep = "\t", header = TRUE) %>%
#replace_na(list(log2FoldChange = 0, padj = 1)) %>%
#mutate(contrast_type = "comparison",
#contrast = "MPXV_vs_mock@24h",
#log2FoldChange = as.numeric(log2FoldChange),
#padj = as.numeric(padj)) %>%
#dplyr::inner_join(object_rna_contrast_thresholds.df) %>%
#  dplyr::mutate(is_signif = padj <= p_value_threshold,
#                is_hit = is_signif & abs(log2FoldChange) >= log2FC_threshold,
#                p_value_compressed = 10^(-sapply(-log10(padj), mlog10_pvalue_compress)),
#                show_label = coalesce(is_hit, FALSE),
#                log2FC_trunc = pmax(-log2FC_max, pmin(log2FC_max, log2FoldChange)),
#                truncation = scatter_truncation(log2FoldChange, log2FC_trunc, padj, padj, is_hit | !is_signif),
#                truncation_type = point_truncation_type(truncation, is_signif)
#                ) %>%
#  dplyr::mutate(show_label = if_else(rep.int(sum(show_label) >= 400L, n()), is_hit, show_label)) 

#object_thomas.df %>%
object_contrasts_4show.df %>%
#object_contrasts_4show_test.df %>% 
  group_by(contrast) %>% do({
    sel_object_contrast.df <- .
    sel_object_contrast_thresholds.df <- object_rna_contrast_thresholds.df
    message("Plotting ", sel_object_contrast.df$contrast[[1]])
    nlabels <- nrow(dplyr::filter(sel_object_contrast.df, is_signif & show_label))
    
    p <- ggplot(sel_object_contrast.df,
                aes(x=log2FC_trunc, y=p_value_compressed, shape=truncation, size=truncation_type)) +
      geom_hline(data=sel_object_contrast_thresholds.df,
                 aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
      #geom_hline(data=sel_object_contrast_thresholds.df,
      #           aes(yintercept = p_value_max), linetype=1, color="darkgray") +
      geom_vline(data=sel_object_contrast_thresholds.df,
                 aes(xintercept = contrast_offset_log2), linetype=1, color="darkgray") +
      geom_vline(data=sel_object_contrast_thresholds.df,
                 aes(xintercept = contrast_offset_log2 + log2FC_threshold), linetype=2, color="darkgray") +
      geom_vline(data=sel_object_contrast_thresholds.df,
                 aes(xintercept = contrast_offset_log2 - log2FC_threshold), linetype=2, color="darkgray") +
      geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_signif),
                      alpha=0.1, size=0.5, color="darkgray") +
      geom_point_rast(data=dplyr::filter(sel_object_contrast.df, is_signif & !is_hit),
                      size=0.5, color="darkgray") +
      geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & is_hit)) +
      geom_text_repel(data=dplyr::filter(sel_object_contrast.df, is_hit & show_label),
                       aes(label = GeneID),
                       #vjust=-0.01,
                       size=if_else(nlabels > 20, 2.5, 3.5),
                       force=if_else(nlabels > 20, 0.25, 1),
                       label.padding=if_else(nlabels > 20, 0.1, 0.25),
                       max.overlaps = Inf,
                       show.legend = FALSE, segment.color = "gray") +
      scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
      #scale_fill_gradient(low="gray75", high="black") +
      #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
      scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
      scale_size_manual(values=point_truncation_size_palette, guide="none") +
      #facet_grid(p_value_range ~ contrast, scales = "free_y") +
      ggtitle(sel_object_contrast.df$contrast[[1]]) +
      theme_bw_ast()
    volcano_plot_path <- file.path(plot_path, "volcanos_rna_contrasts")
    if (!dir.exists(volcano_plot_path)) { dir.create(volcano_plot_path) }
    ggsave(filename = file.path(volcano_plot_path,
                                str_c(project_id, '_', analysis_version, '_volcano_',
                                      str_replace_all(sel_object_contrast.df$contrast[[1]],":|@", "_"), '.pdf')),
           plot = p, width=15, height=18, device=cairo_pdf, family="Segoe UI Symbol")
    tibble()
  })

#Example to visualise certain contrasts via glimma
object_contrasts_4show.df %>% 
  filter(contrast_kind == "treatment_vs_treatment", treatment_rhs == "mock") %>% 
  group_by(contrast) %>% do({
    sel_object_contrast.df <- .
    message("Plotting ", sel_object_contrast.df$contrast[[1]])
    
    glMDPlot(sel_object_contrast.df,
         main= paste0("MD plot:", sel_object_contrast.df$contrast[[1]]),
         xval="log2MeanNormCount",
         yval="log2FoldChange",
         counts=counts(dds),
         anno=select(sel_object_contrast.df, GeneID),
         groups=dds$timepoint,
         samples=colnames(dds),
         status= sel_object_contrast.df$status,
         sample.cols = case_when(dds$treatment == "mock" ~ "gray", 
                                 dds$treatment == "MPXV" ~"#F9CB40"),
         launch = FALSE,
         display.columns =c("GeneID", "padj", "is_signif", "is_hit"), 
         folder=paste0("MD_plot_", str_replace_all(sel_object_contrast.df$contrast[[1]],":|@", "_")), path = file.path(plot_path, "MD_plots"))
    tibble()
         })

#Export a table of hits (the full table is too large)
treatment_vs_treatment_long <- object_contrasts_4show.df %>% 
  filter(contrast_kind == "treatment_vs_treatment")
time_vs_time_long <- object_contrasts_4show.df %>% 
  filter(contrast_kind == "timepoint_vs_timepoint") 

hits2export_vsmock_summary <- treatment_vs_treatment_long %>% 
  filter(treatment_rhs == "mock", is_hit) %>% 
  select(GeneID, treatment = treatment_lhs, timepoint = timepoint_lhs, change) %>% 
  group_by(GeneID, treatment) %>% 
  summarise(first_sig_timepoint = min(timepoint), all_timepoints = paste0(timepoint, collapse = ","), change = paste0(change, collapse = ",")) %>% 
  ungroup() %>% 
  mutate(change_summary = case_when(str_detect(change, "\\+") & str_detect(change, "\\-") ~ "both",
                                    str_detect(change, "\\+") ~ "+",
                                    str_detect(change, "\\-") ~ "-",
                                    TRUE ~ NA_character_))

write_xlsx(treatment_vs_treatment_long, file.path(analysis_path, "reports", paste0(project_id, '_DESeq2_treatment_vs_treatment_long_', analysis_version, '.xlsx')))
write_xlsx(time_vs_time_long, file.path(analysis_path, "reports", paste0(project_id, '_DESeq2_timepoint_vs_timepoint_long_', analysis_version, '.xlsx')))
write_xlsx(hits2export_vsmock_summary, file.path(analysis_path, "reports", paste0(project_id, '_DESeq2_hits_vsmock_summary_', analysis_version, '.xlsx')))

#Prepare the data for enrichment analysis in Julia ----
#Match the gene_id to protein_acs
ensembl_xrefs.df <- read_tsv(file.path(data_path, str_c(datatype, "_", data_version),"ensembl_xrefs.txt")) #we need this file because the "Gene names" from Rad group are not always the conventional gene names...
uniprot_gene2id <- read_tsv(file.path(data_path,str_c(datatype, "_", data_version), "uniprot_gene2id_20220531.tab")) %>% 
  rename("xref" = "Cross-reference (OpenTargets)") %>% 
  mutate(xref = str_remove(xref, ";$")) %>% 
  separate_rows(xref, sep = ";")

rnaseq_genes_xrefs.df <- DESeq_conditions.df %>%
  select(GeneID) %>% 
  unique() %>% 
  filter(!str_detect(GeneID,"-AS" ),
         !str_detect(GeneID, "-\\d{4,}"),
         !str_detect(GeneID, "AC\\d{3,}"),
         !str_detect(GeneID, "LINC")) %>%
  left_join(select(ensembl_xrefs.df, GeneID = xref_id, gene_ensembl_id)) %>% 
  left_join(select(ensembl_xrefs.df, gene_ensembl_id, xref_id, xref_db)) %>% 
  distinct()

rnaseq_genes2protacs_full.df <- rnaseq_genes_xrefs.df %>% 
  filter(xref_db == "Uniprot") %>% 
  inner_join(uniprot_gene2id, by = c("xref_id" = "Entry")) %>% 
  mutate(lead_gene_name = word(str_remove_all(`Gene names`, ";"), 1),
         ensembl_match = gene_ensembl_id == xref,
         gene_name_match = lead_gene_name == GeneID) %>%
  group_by(GeneID) %>% 
  slice_max(ensembl_match+gene_name_match) %>%
  group_by(xref_id) %>% 
  slice_max(ensembl_match+gene_name_match) %>% 
  ungroup()

rnaseq_genes2protacs.df <- rnaseq_genes2protacs_full.df %>% 
  select(GeneID, "protein_acs" = xref_id) %>% 
  unique() %>% 
  group_by(GeneID) %>% 
  summarise_at(vars(-group_cols()),  str_c, collapse="; ") %>% 
  mutate(lead_protein_ac = word(str_remove_all(protein_acs, ";"), 1))

DESeq_ashr_contrastsXprotacs_trVtr.df <- DESeq_ashr_contrasts.df %>% 
  inner_join(rnaseq_genes2protacs.df) %>% 
  filter(contrast_kind == "treatment_vs_treatment", treatment_rhs == "mock")

#Saving the RData for Julia
save(DESeq_contrasts.df, DESeq_ashr_contrasts.df, DESeq_conditions.df,
    object_contrasts_4show.df,
     conditions.df, conditionsXsamples.df, contrasts.df,
     conditionXeffect.mtx,sampleXeffect.mtx, contrastXeffect.mtx, contrastXmetacondition.mtx,
     rnaseq_genes2protacs.df,
     file = file.path(results_path, str_c(project_id, '_DESeq2_data_', analysis_version, '_julia.RData')))

