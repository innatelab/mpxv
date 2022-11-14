# plotting the volcanos for contrasts
# peptide heatmaps
# how model fits the data, phospho
# Yiqi Huang 2022.08.12

#The following part is only needed if starting from a fresh environment----
project_id <- "mpxv"
data_version <- "20221105"
fit_version <- "20221111"
mstype <- "phospho"
message("Project ID=", project_id, " data version=", data_version)

require(tidyverse)
require(rlang)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

message('Loading data...') 
load(file.path(scratch_path, str_c(project_id, '_msdata_full_', mstype, "_", data_version, '.RData')))
load(file.path(scratch_path, str_c(project_id, '_msglm_data_', mstype, "_", fit_version, '.RData')))
#load(file.path(scratch_path, str_c(project_id, '_msglm_fit_meanfield_', mstype, "_", fit_version, '.RData')))
load(file.path(results_path, str_c(project_id, '_msglm_fit_', mstype, "_", fit_version, '.RData')))

modelobj <- msdata$msentities['object']
quantobj <- msdata$msentities['quantobject']
modelobj_idcol <- paste0(modelobj, "_id")
quantobj_idcol <- paste0(quantobj, "_id")

contrasts.df <- dplyr::ungroup(msglm_def$contrasts) %>%
  dplyr::inner_join(tidyr::pivot_wider(dplyr::mutate(as.data.frame.table(msglm_def$metaconditionXcontrast, responseName="w"),
                                                     side = if_else(w > 0, "lhs", "rhs")) %>% dplyr::filter(w != 0),
                                       c(contrast), names_from = "side", values_from = "metacondition",
                                       names_glue = "{.value}_{side}")) %>%
  dplyr::mutate(offset = 0, offset_prior = 0) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_lhs = condition, treatment_lhs = treatment)) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_rhs = condition, treatment_rhs = treatment))

object_contrasts_thresholds.df <- dplyr::select(contrasts.df, offset, offset_prior, contrast, contrast_type) %>%
  dplyr::mutate(
    p_value_threshold = case_when(contrast_type == "filtering" ~ 1E-2,
                                  contrast_type == "comparison" ~ 1E-2,
                                  TRUE ~ NA_real_),
    median_threshold = case_when(contrast_type == "filtering" ~ pmax(2.0, 2.0 + abs(offset - offset_prior)),
                                 contrast_type == "comparison" ~ pmax(1.0, 0.25 + abs(offset - offset_prior)),
                                 TRUE ~ NA_real_),
    median_max = case_when(contrast_type == "filtering" ~ 12,
                           contrast_type == "comparison" ~ 6,
                           TRUE ~ NA_real_)
  )

#plotting starts from here----
require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)
require(ggpubr)

source(file.path(misc_scripts_path, 'ggplot_ext.R'))
source(file.path(misc_scripts_path, 'furrr_utils.R'))

modelobj_suffix <- "_ptmn"
modelobjs_df <- msdata$objects %>% 
  left_join(select(msdata_full$proteins, protein_ac, protein_description = protein_name)) %>% 
  mutate(protac_label = protein_ac)
obj_labu_shift <- msdata[[str_c(quantobj, "_mscalib")]]$zShift

treatment_palette <- c("mock" = "gray", "MPXV" = "#F9CB40")
hit_palette <- c("non-hit" = "grey", "hit" = "black", "oxidation hit" = "dark grey", "viral hit" = "#F9CB40", "only sig" = "light blue", "is hit nofp" = "dark blue")
fp_treatment_palette <- c(mock="lightgray", MPXV = "gold")
base_font_family <- "Segoe UI Symbol"
base_plot_path <- file.path(analysis_path, 'plots', str_c(data_info$msfolder, "_", fit_version))
sel_ci_target <- "average"


#volcano for all contrasts (check assemble_fits_phospho for the missing steps if starting from the meanfield results)----
object_contrasts_4show.df <- object_contrasts_nofp.df %>%
  select(-contains("threshold")) %>%
  mutate(fp_is_hit = FALSE, fp_median = 0) %>% 
  dplyr::inner_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - offset) >= median_threshold),
                is_hit_nomschecks = is_signif & !is_contaminant,
                is_hit_nofp = is_hit_nomschecks & ((nmsruns_quanted_lhs_max >= 3) | (nmsruns_quanted_rhs_max >= 3)),
                is_hit = is_hit_nofp & (is_viral | !coalesce(fp_is_hit, FALSE) | sign(median) != sign(coalesce(fp_median, 0))), # |
                #(abs(median - fp_median) >= 2)),
                hit_type = case_when(is_signif & !is_hit_nofp ~ "only sig",
                                     is_hit_nofp & !is_hit ~ "is hit nofp",
                                     is_hit & str_detect(object_label, "Oxidation") ~ "oxidation hit",
                                     is_hit & is_viral ~ "viral hit",
                                     is_hit ~ "hit", 
                                     TRUE ~ "non-hit"), 
                is_sel = is_hit_nomschecks,
                mean_trunc = pmax(-median_max, pmin(median_max, mean - offset)) + offset,
                median_trunc = pmax(-median_max, pmin(median_max, median - offset)) + offset,
                truncation = scatter_truncation(median, median_trunc, -log10(p_value), -log10(p_value),
                                                is_hit_nomschecks | !is_signif),
                truncation_type = point_truncation_type(truncation, is_signif),
                show_label = coalesce(is_hit_nomschecks, FALSE),
                short_label = str_remove_all(object_label, "Oxidation_|Phospho_|_M\\d($|\\.\\.\\.)")) %>%
  tidyr::separate(contrast, into = c("group1", "group2"), sep = "_vs_", remove = FALSE) %>%
  dplyr::group_by(contrast, ci_target) %>%
  dplyr::mutate(show_label = if_else(rep.int(sum(show_label & median > 0) >= 400L, n()), is_hit, show_label)) %>% #what does this do?
  dplyr::ungroup()

require(furrr)
plan(multicore, workers = 32)
plot_furrr_opts <- furrr_options(globals = c("base_plot_path", "base_font_family",
                                             "project_id", "data_version", "fit_version", "obj_labu_shift",
                                             "mlog10_trans", "mlog_pow_trans", "mlog_breaks",
                                             "theme_bw_ast",
                                             "point_truncation_shape_palette", "point_truncation_size_palette",
                                             "treatment_palette",
                                             "sel_ci_target", "modelobj", "quantobj", "modelobj_idcol", "quantobj_idcol", "modelobjs_df",
                                             "msglm_def", "fit_stats", "msdata", "msdata_full", "object_contrasts_4show.df"),
                                 packages = c("dplyr", "ggplot2", "Cairo", "ggrepel", "ggrastr", "stringr"),
                                 stdout = TRUE)

group_by(object_contrasts_4show.df, ci_target, contrast,
         offset, median_threshold, p_value_threshold) %>%
  group_walk(.keep = TRUE,
             #future_group_walk(.progress=TRUE, .keep=TRUE, .options=plot_furrr_opts,
             function(sel_object_contrast.df, contrast_info) {
               message("Plotting ", contrast_info$contrast, " ci_target=", contrast_info$ci_target)
               labels_lhs.df <- dplyr::filter(sel_object_contrast.df, (median > offset) & (is_sel | is_signif & show_label))
               if (nrow(labels_lhs.df) > 300) {
                 labels_lhs.df <- dplyr::filter(labels_lhs.df, is_sel)
               }
               labels_rhs.df <- dplyr::filter(sel_object_contrast.df, (median < offset) & (is_sel | is_signif & show_label))
               if (nrow(labels_rhs.df) > 300) {
                 labels_rhs.df <- dplyr::filter(labels_rhs.df, is_sel)
               }
               labels.df <- bind_rows(labels_lhs.df, labels_rhs.df)
               nlabels <- nrow(labels.df)
               
               p <- ggplot(sel_object_contrast.df,
                           aes(x=median_trunc, y=p_value, shape=truncation, size=truncation_type, color=hit_type)) +
                 geom_hline(data=contrast_info, aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
                 #geom_hline(data=contrast_info, aes(yintercept = p_value_max), linetype=1, color="darkgray") +
                 geom_vline(data=contrast_info, aes(xintercept = offset), linetype=1, color="darkgray") +
                 geom_vline(data=contrast_info, aes(xintercept = offset + median_threshold), linetype=2, color="darkgray") +
                 geom_vline(data=contrast_info, aes(xintercept = offset - median_threshold), linetype=2, color="darkgray") +
                 geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_signif),
                                 alpha = 0.1, size = 0.5, color = "darkgray") +
                 geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & !is_hit), shape=1) +
                 geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & is_hit)) +
                 geom_text_repel(data = labels.df,
                                 aes(label = short_label),
                                 size = if_else(nlabels > 20, 2.5, 3.5),
                                 force = if_else(nlabels > 20, 0.25, 1.0),
                                 nudge_y = -0.12,
                                 point.padding = 0.02,
                                 box.padding = if_else(nlabels > 20, 0.1, 0.25),
                                 show.legend = FALSE, segment.color = "gray",
                                 max.overlaps = Inf) +
                 scale_y_continuous(trans = mlog10_trans(), limits = c(1.0, NA)) +
                 #scale_fill_gradient(low="gray75", high="black") +
                 #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
                 scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
                 scale_size_manual(values=point_truncation_size_palette, guide="none") +
                 scale_color_manual(values = hit_palette, na.value = "magenta") +
                 #facet_grid(p_value_range ~ contrast, scales = "free_y") +
                 ggtitle(contrast_info$contrast, subtitle=str_c("ci_target=", contrast_info$ci_target)) +
                 theme_bw_ast(base_family = base_font_family)
               plot_path <- file.path(base_plot_path, str_c("volcanos_contrasts_", contrast_info$ci_target))
               if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
               
               ggsave(filename = file.path(plot_path,
                                           str_c(project_id, '_', fit_version, '_volcano_',
                                                 str_replace_all(contrast_info$contrast, ":|@", "_"), '.pdf')),
                      plot=p, width=15, height=18, device=cairo_pdf, family=base_font_family) #For windows users only: beware of what font you're using! If it's the first time you use this code, you need to import the fonts with "extrafont" package
             })

#Making timecourse for all proteins. The dots are from the original peptide intensity values of MQ.----
#sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(gene_name, "ZC3H13") & ptm_type == "Phospho")  #This is used for debugging
#sel_objects.df <- dplyr::semi_join(modelobjs_df, dplyr::select(fit_diff, object_id), by = "object_id")
sel_objects.df <- dplyr::semi_join(modelobjs_df, dplyr::select(fit_stats$objects, object_id), by="object_id")#This is all!

msdata_full$pepmodstates <- dplyr::mutate(msdata_full$pepmodstates,
                                          pepmodstate_seq = str_c(pepmod_seq, ".", charge))

dplyr::left_join(sel_objects.df, dplyr::select(msdata_full$ptmngroup_idents, ptmngroup_id, n_pepmodstates)) %>% unique() %>% 
  dplyr::group_by(object_id) %>% slice_max(order_by = n_pepmodstates) %>% do({
    sel_obj.df <- .
    sel_ptm_type <- sel_obj.df$ptm_type
    message("Plotting ", sel_obj.df$object_label, " time course")
    
    sel_var <- if_else(sel_ci_target == "average", "obj_cond_labu", "obj_cond_labu_replCI")
    sel_obj_conds.df <- dplyr::inner_join(dplyr::select(sel_obj.df, object_id),
                                             dplyr::filter(fit_stats$object_conditions, var == sel_var)) %>%
      dplyr::inner_join(dplyr::select(msglm_def$conditions, condition, treatment, timepoint, timepoint_num), by="condition") %>%
      dplyr::arrange(timepoint_num) %>% 
      dplyr::mutate(dplyr::across(c(mean, median, starts_with("q")),
                                  ~2^(.x + obj_labu_shift))) %>%
      dplyr::mutate(q97.5_limit = max(median + (q75 - q25)*5))
    
    sel_obj_contrasts.df <- dplyr::semi_join(dplyr::filter(object_contrasts_4show.df, ci_target == sel_ci_target & str_starts(var, "obj_cond_labu"),
                                                           group1 == "MPXV"),
                                             sel_obj.df, by="object_id") %>% 
      mutate(y.position = max(sel_obj_conds.df$`q97.5`)*1.1,
             treatment = group1,
             timepoint = as.numeric(str_extract(group2, "(?<=@)\\d+")),
             #p_value = formatC(p_value, format = "e", digits = 2),
             p_label = case_when(p_value < 0.001 ~ "***",
                                 p_value < 0.01 ~ "**",
                                 p_value < 0.05 ~ "*",
                                 TRUE ~ ""))# for the pvalue labels in the graph
    
    sel_intensity_range.df <- dplyr::group_by(sel_obj_conds.df, object_id) %>% dplyr::summarise(
      intensity_min = 0.5*quantile(q25, 0.05, na.rm=TRUE),
      intensity_max = 1.5*quantile(q75, 0.95, na.rm=TRUE),
      .groups="drop")
    sel_obj_msdata.df <- sel_obj.df %>%
      dplyr::inner_join(sel_intensity_range.df) %>%
      dplyr::inner_join(msdata_full$ptmngroup2pepmodstate) %>%
      dplyr::inner_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id, pepmodstate_seq)) %>% 
      dplyr::inner_join(msdata_full$pepmodstate_intensities) %>%
      dplyr::inner_join(msdata_full$ptmn_locprobs) %>%
      dplyr::inner_join(dplyr::select(msdata$msrun_shifts, msrun, total_msrun_shift)) %>%
      dplyr::mutate(intensity_norm = intensity*2^(-total_msrun_shift),
                    intensity_norm_trunc = pmin(pmax(intensity_norm,intensity_min), intensity_max),
                    intensity_used = intensity_norm_trunc
      ) %>%
      dplyr::inner_join(msdata_full$msruns) %>%
      dplyr::mutate(locprob_valid = coalesce(ptm_locprob, 0) >= data_info$locprob_min,
                    pvalue_valid = coalesce(psm_pvalue, 1) <= data_info$pvalue_ident_max,
                    ms_status = case_when(locprob_valid & pvalue_valid ~ "valid",
                                          !locprob_valid & pvalue_valid ~ "bad PTM loc",
                                          locprob_valid & !pvalue_valid ~ "bad ident",
                                          TRUE ~ 'bad ident&loc')) %>% 
      group_by(object_id) %>% 
      mutate(total_pepmodstates = n_distinct(pepmodstate_id)) %>% 
      ungroup()
    #print(sel_obj_msdata.df)
    if (nrow(sel_obj_conds.df) > 0) {
      p <-
        ggplot(data=sel_obj_conds.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
        geom_ribbon(aes(x = timepoint_num, ymin = `q2.5`, ymax=`q97.5`),
                    alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
        geom_ribbon(aes(x = timepoint_num, ymin = q25, ymax= q75),
                    alpha=0.5, stat = "identity", size=0.5) +
        geom_path(aes(x = timepoint_num, y = median), alpha=0.5, size=1, stat="identity") +
        geom_point(data=sel_obj_msdata.df,
                   aes(y = intensity_used, shape=ms_status),
                   position = position_jitter(width = 0.75, height = 0, seed=12323), size=1.5) +
        geom_text(data=sel_obj_msdata.df,
                  aes(y = intensity_used, label=replicate),
                  position = position_jitter(width = 0.75, height = 0, seed=12323), size=1, color="lightgray", show.legend=FALSE) +
        geom_text(data = sel_obj_contrasts.df, aes(x = timepoint, y = y.position, label = p_label), colour = "black")+
        theme_bw_ast(base_family = "", base_size = 8) +
        scale_x_continuous(breaks=unique(msdata$msruns$timepoint_num)) +
        scale_color_manual(values=treatment_palette) +
        scale_shape_manual(values=c("valid"=19, "bad PTM loc"=8, "bad ident"=1, "bad ident&loc"=4)) +
        scale_fill_manual(values=treatment_palette) +
        scale_y_log10("PTM Intensity") +
        ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
                subtitle=str_c(sel_obj_msdata.df$ptmns, "\n",
                sel_obj.df$protein_description, " (total npepmodstates=", unique(sel_obj_msdata.df$total_pepmodstates), ")")) +
        facet_wrap( ~ object_label+pepmodstate_seq, ncol =1, scales = "free")
      plot_path <- file.path(base_plot_path,
                             str_c("timecourse_", sel_ptm_type, "_", sel_ci_target,
                                   modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")))
      if (!dir.exists(plot_path)) {dir.create(plot_path, recursive = TRUE)}
      ggsave(p, file = file.path(plot_path, str_c(project_id, "_", data_info$msfolder, '_', fit_version, "_",
                                                  str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], ".pdf")),
             width=8, height=6, device = cairo_pdf)
    }
    tibble()
  })

#include fp data for the timecourse----
require(multidplyr)
plot_cluster <- new_cluster(16)
cluster_library(plot_cluster, c("dplyr", "ggplot2", "stringr", "Cairo", "readr", "rlang"))
cluster_copy(plot_cluster, c("data_info", "analysis_path", "base_plot_path", "base_font_family",
                             "project_id", "data_version", "fit_version", "obj_labu_shift",
                             "mlog10_trans", "mlog_pow_trans", "mlog_breaks",
                             "theme_bw_ast",
                             "point_truncation_shape_palette", "point_truncation_size_palette",
                             "treatment_palette",
                             "sel_ci_target", "modelobj", "quantobj", "modelobj_idcol", "quantobj_idcol", "modelobjs_df", "modelobj_suffix",
                             #"fp.env", "ptmngroup2protregroup.df", "fp_treatment_palette",
                             "msglm_def", "fit_stats", "msdata", "msdata_full", "object_contrasts_4show.df"))

dplyr::left_join(sel_objects.df, dplyr::select(msdata_full$ptmngroup_idents, ptmngroup_id, n_pepmodstates)) %>% unique() %>% 
  dplyr::group_by(object_id) %>% slice_max(order_by = n_pepmodstates) %>% partition(plot_cluster) %>%
  do({
    sel_obj.df <- .
    sel_ptm_type <- sel_obj.df$ptm_type
    message("Plotting ", sel_obj.df$object_label, " time course")
    
    sel_var <- if_else(sel_ci_target == "average", "obj_cond_labu", "obj_cond_labu_replCI")
    sel_obj_conds.df <- dplyr::inner_join(dplyr::select(sel_obj.df, object_id),
                                          dplyr::filter(fit_stats$object_conditions, var == sel_var)) %>%
      dplyr::inner_join(dplyr::select(msglm_def$conditions, condition, treatment, timepoint, timepoint_num), by="condition") %>%
      dplyr::arrange(timepoint_num) %>% 
      dplyr::mutate(dplyr::across(c(mean, median, starts_with("q")),
                                  ~2^(.x + obj_labu_shift))) %>%
      dplyr::mutate(q97.5_limit = max(median + (q75 - q25)*5))
    
    sel_obj_contrasts.df <- dplyr::semi_join(dplyr::filter(object_contrasts_4show.df, ci_target == sel_ci_target & str_starts(var, "obj_cond_labu"),
                                                           group1 == "MPXV"),
                                             sel_obj.df, by="object_id") %>% 
      mutate(y.position = log10(max(sel_obj_conds.df$`q97.5`)*1.1),
             treatment = group1,
             time = str_extract(group2, "(?<=@)\\d+"),
             p_value = formatC(p_value, format = "e", digits = 2))# for the pvalue labels in the graph
    
    sel_intensity_range.df <- dplyr::group_by(sel_obj_conds.df, object_id) %>% dplyr::summarise(
      intensity_min = 0.5*quantile(q25, 0.05, na.rm=TRUE),
      intensity_max = 1.5*quantile(q75, 0.95, na.rm=TRUE),
      .groups="drop")
    sel_obj_msdata.df <- sel_obj.df %>%
      dplyr::inner_join(sel_intensity_range.df) %>%
      dplyr::inner_join(msdata_full$ptmngroup2pepmodstate) %>%
      dplyr::inner_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id, pepmodstate_seq)) %>% 
      dplyr::inner_join(msdata_full$pepmodstate_intensities) %>%
      dplyr::inner_join(msdata_full$ptmn_locprobs) %>%
      dplyr::inner_join(dplyr::select(msdata$msrun_shifts, msrun, total_msrun_shift)) %>%
      dplyr::mutate(#intensity_norm_orig = intensity,
        #intensity_norm_orig_scaled = intensity_norm_orig*exp(-qobj_shift),
        intensity_norm = intensity*exp(-total_msrun_shift),
        #intensity_norm_scaled = intensity_norm*exp(-qobj_shift),
        #intensity_norm_orig_scaled_trunc = pmin(pmax(intensity_norm_orig_scaled,
        #                                             intensity_min), intensity_max),
        intensity_norm_trunc = pmin(pmax(intensity_norm,intensity_min), intensity_max),
        intensity_used = intensity_norm_trunc
      ) %>%
      dplyr::inner_join(msdata_full$msruns) %>%
      dplyr::mutate(locprob_valid = coalesce(ptm_locprob, 0) >= data_info$locprob_min,
                    pvalue_valid = coalesce(psm_pvalue, 1) <= data_info$pvalue_ident_max,
                    ms_status = case_when(locprob_valid & pvalue_valid ~ "valid",
                                          !locprob_valid & pvalue_valid ~ "bad PTM loc",
                                          locprob_valid & !pvalue_valid ~ "bad ident",
                                          TRUE ~ 'bad ident&loc')) %>% 
      group_by(object_id) %>% 
      mutate(total_pepmodstates = n_distinct(pepmodstate_id)) %>% 
      ungroup()
    #print(sel_obj_msdata.df)
    if (nrow(sel_obj_conds.df) > 0) {
      p <-
        ggplot(data=sel_obj_conds.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
        geom_ribbon(aes(x = timepoint_num, ymin = `q2.5`, ymax=`q97.5`),
                    alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
        geom_ribbon(aes(x = timepoint_num, ymin = q25, ymax= q75),
                    alpha=0.5, stat = "identity", size=0.5) +
        geom_path(aes(x = timepoint_num, y = median), alpha=0.5, size=1, stat="identity") +
        geom_point(data=sel_obj_msdata.df,
                   aes(y = intensity_used, shape=ms_status),
                   position = position_jitter(width = 0.75, height = 0, seed=12323), size=1.5) +
        geom_text(data=sel_obj_msdata.df,
                  aes(y = intensity_used, label=replicate),
                  position = position_jitter(width = 0.75, height = 0, seed=12323), size=1, color="lightgray", show.legend=FALSE) +
        theme_bw_ast(base_family = "", base_size = 8) +
        scale_x_continuous(breaks=unique(msdata$msruns$timepoint_num)) +
        scale_color_manual(values=treatment_palette) +
        scale_shape_manual(values=c("valid"=19, "bad PTM loc"=8, "bad ident"=1, "bad ident&loc"=4)) +
        scale_fill_manual(values=treatment_palette) +
        scale_y_log10("PTM Intensity") +
        facet_wrap( ~ object_label+pepmodstate_seq, ncol =1, scales = "free")
      
      h <- as.integer(2+3*n_distinct(sel_obj_msdata.df$pepmodstate_id))
      
      if (exists("fp.env")) {
        sel_fp_conds.df <- semi_join(fp.env$fit_stats$object_conditions, semi_join(ptmngroup2protregroup.df, sel_obj.df), by = c("protregroup_id" = "fp_protregroup_id")) %>%
          filter(var == if_else(sel_ci_target == "average", "obj_cond_labu", "obj_cond_labu_replCI")) %>%
          dplyr::inner_join(dplyr::select(msglm_def$conditions, condition, treatment, timepoint, timepoint_num), by="condition")
        if (nrow(sel_fp_conds.df) > 0L) {
          fp_plot <- ggplot(data=sel_fp_conds.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
            geom_ribbon(aes(x = timepoint_num, ymin = `q2.5`, ymax=`q97.5`),
                        alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
            geom_ribbon(aes(x = timepoint_num, ymin = q25, ymax=q75),
                        alpha=0.5, stat = "identity", size=0.5) +
            geom_path(aes(x = timepoint_num, y = median), alpha=0.5, size=1, stat="identity") +
            theme_bw_ast(base_family = "", base_size = 8) +
            scale_x_continuous(breaks=unique(msdata_full$msruns$timepoint_num)) +
            scale_color_manual(values=fp_treatment_palette) +
            scale_fill_manual(values=fp_treatment_palette) +
            ylab("Protein Intensity")
          p <- ggpubr::ggarrange(p, fp_plot, ncol=1, heights=c(0.6, 0.3))
          h <- h + 2L
        }
      }
      p <- p + ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
                       subtitle=str_c(sel_obj_msdata.df$ptmns, "\n",
                       sel_obj.df$protein_description, " (total npepmodstates=", unique(sel_obj_msdata.df$total_pepmodstates), ")"))
      plot_path <- file.path(base_plot_path,
                             str_c("timecourse_", sel_ptm_type, "_", sel_ci_target,
                                   modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")))
      if (!dir.exists(plot_path)) {dir.create(plot_path, recursive = TRUE)}
      ggsave(p, file = file.path(plot_path, str_c(project_id, "_", data_info$msfolder, '_', fit_version, "_",
                                                  str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], "_long.pdf")),
             width=8, height=h, limitsize=FALSE, device = cairo_pdf)
    }
    tibble()
  })

#plot phosphosites on viral proteins----
plot_version <- 20221113
ptm_pvalue_ident_max <- 1E-3
ptm_pvalue_quant_max <- 1E-2
ptm_locprob_ident_min <- 0.75
ptm_locprob_quant_min <- 0.5

ptm_status_levels <- c("N/A", "potential", "low conf.", "observed")

ptmngroupXcondition_stats.df <- dplyr::select(msdata_full$ptmngroups, is_viral, ptm_type, ptmngroup_id, ptmn_id) %>% dplyr::filter(is_viral) %>%
  dplyr::left_join(msdata_full$ptmn2pepmodstate) %>%
  dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
  dplyr::left_join(dplyr::select(msdata_full$ptmn_locprobs, evidence_id, ptmn_id, ptm_locprob, msrun, ptm_AA_seq)) %>%
  dplyr::left_join(msdata_full$msruns) %>%
  dplyr::group_by(ptm_type, ptmn_id, ptmngroup_id, condition) %>%
  dplyr::summarise(ptm_pvalue_min = min(psm_pvalue, na.rm=TRUE),
                   ptm_locprob_max = max(ptm_locprob, na.rm=TRUE),
                   n_idented_and_localized  = n_distinct(str_c(msrun, " ", pepmodstate_id)[(coalesce(psm_pvalue, 1) <= ptm_pvalue_ident_max) &
                                                                                             (coalesce(ptm_locprob, 0) >= ptm_locprob_ident_min)]),
                   n_quanted  = n_distinct(str_c(msrun, " ", pepmodstate_id)[(coalesce(psm_pvalue, 1) <= ptm_pvalue_quant_max) &
                                                                               (coalesce(ptm_locprob, 0) >= ptm_locprob_quant_min)]),
                   .groups="drop") %>% 
  left_join(msdata_full$ptm2protein) %>% 
  left_join(select(msdata_full$proteins, protein_ac, protein_name, seq)) %>% 
  mutate(seq_length = str_length(seq)) %>% 
  filter(!is.na(ptm_AA_seq))


viral_phospho_intensities.df <- dplyr::filter(fit_stats$object_conditions, var == "obj_cond_labu", is_viral, ptm_type == "Phospho", var == "obj_cond_labu") %>% 
  dplyr::inner_join(dplyr::select(msglm_def$conditions, condition, treatment, timepoint, timepoint_num), by="condition") %>%
  filter(treatment == "MPXV") %>% 
  dplyr::left_join(ptmngroupXcondition_stats.df) %>% 
  dplyr::mutate(ptm_correct = (ptm_type == "Phospho" & ptm_AA_seq %in% c("S","T","Y")),
                is_observed = coalesce(n_idented_and_localized, 0) > 0,
                is_observed_lowconf = coalesce(n_quanted) > 0,
                ptm_status = case_when(is_observed ~ "observed",
                                       is_observed_lowconf ~ "low conf.",
                                       ptm_correct ~ "potential",
                                       TRUE ~ "N/A") %>% factor(ordered=TRUE, levels=ptm_status_levels),
                identifier = protein_ac) %>%
    dplyr::mutate(timepoint_label = factor(str_c(timepoint, "h p.i."), levels=str_c(sort(unique(msdata_full$msruns$timepoint_num)), "h p.i.")),
                shown_ptm_pos =  ptm_pos,
                shown_median_log2 = pmax(if_else(!is.na(median) & ptm_status %in% c("observed", "low conf."),
                                                 obj_labu_shift + median, 3.5), 3.5),
                ptm_label = str_c(ptm_AA_seq, ptm_pos),
                ptm_status = as.character(ptm_status),
                ptm_typeXstatus = case_when(ptm_status == "observed" ~ ptm_type,
                                            ptm_status == "N/A" ~ ptm_status,
                                            TRUE ~ str_c(ptm_type, " (", ptm_status, ")"))) %>% 
  filter(!is.na(ptm_label)) %>% 
  dplyr::arrange(timepoint_num)

ptm_palette <- c("Phospho" = "#5e268f", "N/A" = "gray") 
ptm_typeXstatus_shape_palette <- c("Phospho"=22L, "Phospho (low conf.)"=22L, "Phospho (potential)" = 22L)#, "N/A" = 4L)
ptm_typeXstatus_fill_palette <- c("Phospho"="#5e268f", "Phospho (low conf.)"="#c29ae5", "Phospho (potential)" = "white")#, "N/A" = NA)
ptm_status_width_palette <- c("observed"=0.5, "low conf."=0.5, "potential" = 0.5)#, "N/A" = 0.5)
domain_type_palette <- c("localization" = "lemonchiffon", "function" = "palegreen3")
domain_type_color_palette <- c("localisation" = "lemonchiffon4", "function" = "palegreen4")

dplyr::group_by(viral_phospho_intensities.df, protein_ac) %>%
  group_walk(function(sel_obj_conds.df, ...){
    viral_gene <- sel_obj_conds.df$genename[[1]]
    identifier <- sel_obj_conds.df$identifier[[1]]
    message("Plotting PTMs of ", viral_gene)
    p <- ggplot(sel_obj_conds.df, aes(x = ptm_pos, y = shown_median_log2)) +
      geom_segment(aes(xend = ptm_pos, y = 1, yend = shown_median_log2,
                       color= ptm_typeXstatus))+
      geom_point(aes(shape=ptm_typeXstatus, fill=ptm_typeXstatus), color="gray", size=2.5) +
      geom_hline(yintercept = 1, colour = "black", size = 3)+
      #scale_color_manual("Virus", values = ptm_typeXstatus_color_palette) + new_scale_color() +
      geom_text_repel(data=sel_obj_conds.df,
                      aes(label = ptm_label, color=ptm_type),
                      min.segment.length=0.2, segment.alpha=0.5, segment.size=0.25, nudge_y = 1, size=4, max.overlaps = Inf) +
      scale_color_manual("PTM", values = ptm_palette) +
      scale_fill_manual("PTM", values = ptm_typeXstatus_fill_palette) +
      scale_shape_manual("PTM", values = ptm_typeXstatus_shape_palette) +
      scale_y_continuous(expression(log[2](Intensity), parse=TRUE)) +
      scale_x_continuous("Aligned Position", limits=c(0L, sel_obj_conds.df$seq_length[[1]])) +
      facet_grid(timepoint_label ~ .) +
      theme_bw_ast(base_family = "") +
      theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line = element_blank()) +
      guides(colour = "none")+
      ggtitle(str_c("PTMs on ", viral_gene, " protein, ", identifier),
              subtitle=sel_obj_conds.df$protein_name)
    plot_path <- file.path(analysis_path, "plots", str_c(mstype, '_', data_version, "_", fit_version), str_c("viral_ptm_", plot_version))
    if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
    ggsave(p, file = file.path(plot_path, str_c(project_id, "_", mstype, '_', fit_version, "_viral_",
                                                viral_gene, "_", plot_version,
                                                ".pdf")),
           width=12+nrow(sel_obj_conds.df)/15, height=7, device = cairo_pdf, family="Arial")
  })


