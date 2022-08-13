# plotting the volcanos for contrasts
# peptide heatmaps
# how model fits the data, phospho
# Yiqi Huang 2022.08.12

#The following part is only needed if starting from a fresh environment----
project_id <- "mpxv"
data_version <- "20220812"
fit_version <- "20220813"
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
modelobj_suffix <- "_ptmn"
modelobjs_df <- msdata$objects %>% 
  left_join(select(msdata_full$proteins, protein_ac, protein_description = protein_name)) %>% 
  mutate(protac_label = protein_ac)
modelobj_idcol <- paste0(modelobj, "_id")
quantobj_idcol <- paste0(quantobj, "_id")
obj_labu_shift <- msdata[[str_c(quantobj, "_mscalib")]]$zShift

#plotting starts from here----
require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)
require(ggpubr)

source(file.path(misc_scripts_path, 'ggplot_ext.R'))
source(file.path(misc_scripts_path, 'furrr_utils.R'))

treatment_palette <- c("mock" = "gray", "MPXV" = "#F9CB40")
hit_palette <- c("non-hit" = "grey", "hit" = "black", "oxidation hit" = "dark grey", "viral hit" = "#F9CB40" )
base_font_family <- "Segoe UI Symbol"
base_plot_path <- file.path(analysis_path, 'plots', str_c(data_info$msfolder, "_", fit_version))
sel_ci_target <- "average"

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
    p_value_threshold = case_when(contrast_type == "filtering" ~ 1E-3,
                                  contrast_type == "comparison" ~ 1E-3,
                                  TRUE ~ NA_real_),
    median_threshold = case_when(contrast_type == "filtering" ~ pmax(2.0, 2.0 + abs(offset - offset_prior)),
                                 contrast_type == "comparison" ~ pmax(0.5, 0.25 + abs(offset - offset_prior)),
                                 TRUE ~ NA_real_),
    median_max = case_when(contrast_type == "filtering" ~ 12,
                           contrast_type == "comparison" ~ 6,
                           TRUE ~ NA_real_)
  )

#volcano for all contrasts----
object_contrasts_4show.df <- fit_contrasts$object_conditions %>%
  dplyr::filter(var %in% c("obj_cond_labu", "obj_cond_labu_replCI")) %>%
  dplyr::ungroup() %>%
  select(-contains("threshold")) %>%
  dplyr::inner_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median - offset) >= median_threshold),
                is_sel = FALSE,
                is_hit_nomschecks = is_signif & !is_contaminant,
                is_hit = is_hit_nomschecks,
                hit_type = case_when(is_hit & str_detect(object_label, "Oxidation") ~ "oxidation hit",
                                     is_hit & is_viral ~ "viral hit",
                                     is_hit ~ "hit", TRUE ~ "non-hit"),
                mean_trunc = pmax(-median_max, pmin(median_max, mean - offset)) + offset,
                median_trunc = pmax(-median_max, pmin(median_max, median - offset)) + offset,
                truncation = scatter_truncation(median, median_trunc, -log10(p_value), -log10(p_value),
                                                is_hit_nomschecks | !is_signif),
                truncation_type = point_truncation_type(truncation, is_signif),
                show_label = coalesce(is_hit_nomschecks, FALSE),
                short_label = str_remove_all(object_label, "Oxidation_|Phospho_|_M\\d$")) %>%
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

#Making timecourse for all proteins. The dots are from the original LFQ values of MQ.----
sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(gene_name, "RBMX") & ptm_type == "Phospho")  #This is used for debugging
sel_objects.df <- dplyr::semi_join(modelobjs_df, dplyr::select(fit_stats$objects, object_id), by="object_id")#This is all!

msdata_full$pepmodstates <- dplyr::mutate(msdata_full$pepmodstates,
                                          pepmodstate_seq = str_c(pepmod_seq, ".", charge))
fp_treatment_palette <- c(mock="lightgray", MPXV = "gold")

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
                                          TRUE ~ 'bad ident&loc'))
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
        scale_y_log10("log10(PTM Intensity)") +
        ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
                subtitle=str_c(sel_obj.df$protein_description, " (npepmodstates=", sel_obj.df$n_pepmodstates, ")")) +
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






require(multidplyr)
plot_cluster <- new_cluster(16)
cluster_library(plot_cluster, c("dplyr", "ggplot2", "stringr", "Cairo", "readr", "rlang"))
cluster_copy(plot_cluster, c("conditions.df", "fit_stats", "msdata_full", "theme_bw_ast",
                             "sel_std_type", "treatment_palette", "fp_treatment_palette",
                             "global_labu_shift", "project_id", "data_info",
                             "analysis_path", "msfolder", "fit_version", "modelobj_suffix", "fp.env", "ptm2protregroup.df"
                             
                             "base_plot_path", "base_font_family",
                             "project_id", "data_version", "fit_version", "obj_labu_shift",
                             "mlog10_trans", "mlog_pow_trans", "mlog_breaks",
                             "theme_bw_ast",
                             "point_truncation_shape_palette", "point_truncation_size_palette",
                             "treatment_palette",
                             "sel_ci_target", "modelobj", "quantobj", "modelobj_idcol", "quantobj_idcol", "modelobjs_df",
                             "msglm_def", "fit_stats", "msdata", "msdata_full", "object_contrasts_4show.df"))

dplyr::left_join(sel_objects.df, dplyr::select(msdata_full$ptmn_stats, ptmn_id, n_pepmodstates)) %>%
  dplyr::left_join(dplyr::select(msdata_full$proteins, protein_ac, protein_description=protein_name)) %>%
  group_by(object_id) %>% partition(plot_cluster) %>%
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
                                          TRUE ~ 'bad ident&loc'))
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
        scale_y_log10("log10(PTM Intensity)") +
        ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
                subtitle=str_c(sel_obj.df$protein_description, " (npepmodstates=", sel_obj.df$n_pepmodstates, ")")) +
        facet_wrap( ~ object_label+pepmodstate_seq, ncol =1, scales = "free")
      
      h <- 2+3*n_distinct(sel_obj_iactions.df$pepmodstate_id)
      
      if (exists("fp.env")) {
        sel_fp_conds.df <- semi_join(fp.env$fit_stats$object_conditions, semi_join(ptm2protregroup.df, sel_obj.df)) %>%
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
            scale_y_log10("log10(Protein Intensity)")
          p <- ggpubr::ggarrange(p, fp_plot, ncol=1, heights=c(0.6, 0.3))
          h <- h + 2L
        }
      }
      p <- p + ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
                       subtitle=str_c(sel_obj.df$protein_description, " (npepmodstates=", sel_obj.df$n_pepmodstates, ")"))
      if (!dir.exists(plot_path)) {dir.create(plot_path, recursive = TRUE)}
      ggsave(p, file = file.path(plot_path, str_c(project_id, "_", data_info$msfolder, '_', fit_version, "_",
                                                  str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], "_long.pdf")),
             width=8, height=h, limitsize=FALSE, device = cairo_pdf)
    }
    tibble()
  })
