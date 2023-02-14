# assembly and analysis of mpxc phospho data, phospho
# Experiments done in Aug 2022
# Author: Yiqi Huang
###############################################################################

project_id <- 'mpxv'
message('Project ID=', project_id)
data_version <- "20221105"
fit_version <- "20221111"
mstype <- "phospho"
message("Assembling fit results for project ", project_id,
        " (dataset v", data_version, ", fit v", fit_version, ")")

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))

require(msglm)
require(tidyverse)
require(furrr)

message('Loading data...')
load(file.path(scratch_path, paste0(project_id, '_msglm_data_', mstype, "_", fit_version, '.RData')))
load(file.path(scratch_path, paste0(project_id, '_msdata_full_', mstype, "_", data_version, '.RData')))

message('Loading MSGLM model fit results...')

modelobj <- msdata$msentities[["object"]]
quantobj <- msdata$msentities[["quantobject"]]
modelobjs_df <- msdata$objects
modelobj_idcol <- paste0(modelobj, "_id")
quantobj_idcol <- paste0(quantobj, "_id")

fit_path <- file.path(scratch_path, paste0("msglm_", project_id, "_", mstype))#_', fit_version))
fit_files <- list.files(fit_path, paste0("msglm_", project_id, "_", mstype, "_", fit_version, '_\\d+\\.RData'))
message('Found ', length(fit_files), ' model file(s)')
fit_files.df <- tibble(filename = as.character(fit_files)) %>%
  tidyr::extract(filename, "chunk", ".+_(\\d+).RData$", convert=TRUE, remove=FALSE) %>%
  dplyr::left_join(dplyr::select(modelobjs_df, chunk, object_id), by="chunk") %>%
  dplyr::arrange(chunk)
require(RMySQL)
chunk_dispatcher_conn <- dbConnect(RMySQL::MySQL(),
                                   dbname="inlab_computing", user="inlab_dispatcher",
                                   host="tumevi4-websrv1.srv.mwn.de",
                                   password=Sys.getenv("INLAB_DISPATCHER_PASSWD"),
                                   timeout=300)
chunk_statuses.df <- dbGetQuery(chunk_dispatcher_conn,
                                str_c("SELECT * FROM jobchunks WHERE user='ge54heq2' AND job_id='",
                                      project_id, "_", mstype, "_", fit_version, "'")) %>%
  dplyr::mutate(start_time = as.POSIXlt(start_time),
                end_time = as.POSIXlt(end_time),
                fit_time = end_time - start_time)
dbDisconnect(chunk_dispatcher_conn)
table(chunk_statuses.df$status)
fit_files.df <- dplyr::inner_join(fit_files.df,
                                  dplyr::filter(chunk_statuses.df, status == "complete") %>%  #& fit_time > 0) %>%
                                    dplyr::select(chunk, status, fit_time), by="chunk")
#fit_diff <- anti_join(fit_files.df, select(fit_stats$global, object_id))

# load fit results in parallel ----
plan(multisession, workers = 32)
fit_chunks <- seq_len(nrow(fit_files.df)) %>%
  furrr::future_map(.progress = TRUE, .options = furrr_options(stdout=FALSE, packages=c("msglm"), globals=c("fit_path", "fit_files.df")),
                    ~load_fit_chunk(.x))
names(fit_chunks) <- purrr::map_chr(fit_chunks, ~paste0(.$results_info$fit_version, '_',
                                                        .$msglm_results$objects$stats$object_id[1]))

fit_stats <- combine_fit_chunks(fit_chunks, 'stats')
fit_contrasts <- combine_fit_chunks(fit_chunks, 'contrast_stats')

rm(fit_chunks)

message('Done.')

# additional filter for hits ----
obj_conditions.df <- tidyr::expand(fit_stats$object_conditions, ptmngroup_id, condition) %>%
  dplyr::inner_join(msdata$msruns) %>%
  dplyr::left_join(msdata_full$ptmngroup2pepmodstate) %>%
  dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
  dplyr::mutate(is_quanted = !is.na(intensity),
                is_idented = (coalesce(psm_pvalue, 1) <= data_info$pvalue_ident_max)) %>%
  dplyr::group_by(condition, ptmngroup_id) %>%
  dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_quanted]), # count msruns
                   nmsruns_idented = n_distinct(msrun[is_idented])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(object_id = ptmngroup_id)

contrastXcondition.df <- as_tibble(as.table(msglm_def$conditionXmetacondition)) %>% dplyr::filter(n != 0) %>%
  dplyr::select(-n) %>%
  dplyr::inner_join(as_tibble(as.table(msglm_def$metaconditionXcontrast))) %>% dplyr::filter(n != 0) %>% 
  dplyr::arrange(contrast, metacondition, condition)

pre_object_contrasts.df <- dplyr::inner_join(obj_conditions.df, contrastXcondition.df) %>%
  dplyr::mutate(is_lhs = n > 0) %>%
  dplyr::group_by(object_id, contrast, is_lhs) %>%
  dplyr::summarise(has_quanted = any(nmsruns_quanted != 0),
                   nmsruns_quanted_min = min(nmsruns_quanted, na.rm=TRUE),
                   nmsruns_quanted_max = max(nmsruns_quanted, na.rm=TRUE),
                   has_idented = any(nmsruns_idented != 0),
                   nmsruns_idented_min = min(nmsruns_idented, na.rm=TRUE),
                   nmsruns_idented_max = max(nmsruns_idented, na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(nmsruns_quanted_min = if_else(has_quanted, nmsruns_quanted_min, 0L),
                nmsruns_quanted_max = if_else(has_quanted, nmsruns_quanted_max, 0L),
                nmsruns_idented_min = if_else(has_idented, nmsruns_idented_min, 0L),
                nmsruns_idented_max = if_else(has_idented, nmsruns_idented_max, 0L)) %>%
  dplyr::group_by(object_id, contrast) %>%
  dplyr::summarise(nmsruns_quanted_lhs_min = nmsruns_quanted_min[is_lhs],
                   nmsruns_quanted_lhs_max = nmsruns_quanted_max[is_lhs],
                   nmsruns_idented_lhs_min = nmsruns_idented_min[is_lhs],
                   nmsruns_idented_lhs_max = nmsruns_idented_max[is_lhs],
                   nmsruns_quanted_rhs_min = nmsruns_quanted_min[!is_lhs],
                   nmsruns_quanted_rhs_max = nmsruns_quanted_max[!is_lhs],
                   nmsruns_idented_rhs_min = nmsruns_idented_min[!is_lhs],
                   nmsruns_idented_rhs_max = nmsruns_idented_max[!is_lhs]) %>%
  dplyr::ungroup() # In our case, the condition and metacondition is the same, so the nmsruns_quanted/idented will always have the same min and max, but it's not always the case if there're multiple conditions in one metacondition (e.g. in AP-MS analysis)

contrasts.df <- dplyr::ungroup(msglm_def$contrasts) %>%
  dplyr::inner_join(tidyr::pivot_wider(dplyr::mutate(as.data.frame.table(msglm_def$metaconditionXcontrast, responseName="w"),
                                                     side = if_else(w > 0, "lhs", "rhs")) %>% dplyr::filter(w != 0),
                                       c(contrast), names_from = "side", values_from = "metacondition",
                                       names_glue = "{.value}_{side}")) %>%
  dplyr::mutate(offset = 0, offset_prior = 0) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint)) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_rhs = condition, treatment_rhs = treatment, timepoint_rhs = timepoint)) %>% 
  dplyr::mutate(contrast_kind = ifelse(treatment_lhs == treatment_rhs, "timepoint_vs_timepoint", "treatment_vs_treatment"))

object_contrasts_thresholds.df <- contrasts.df %>%
  dplyr::mutate(
    p_value_threshold = case_when(contrast_type == "filtering" ~ 1E-2,
                                  contrast_type == "comparison" ~ 1E-2,
                                  TRUE ~ NA_real_),
    median_threshold = case_when(contrast_type == "filtering" ~ pmax(2.0, 2.0 + abs(offset - offset_prior)),
                                 contrast_type == "comparison" ~ pmax(0.5, 0.25 + abs(offset - offset_prior)),
                                 TRUE ~ NA_real_),
    median_max = case_when(contrast_type == "filtering" ~ 12,
                           contrast_type == "comparison" ~ 6,
                           TRUE ~ NA_real_)
  )

object_contrasts_nofp.df <- fit_contrasts$object_conditions %>%
  left_join(pre_object_contrasts.df) %>% 
  dplyr::filter(var %in% c("obj_cond_labu", "obj_cond_labu_replCI")) %>%
  dplyr::inner_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_valid_comparison = (nmsruns_quanted_lhs_max >= 3) | (nmsruns_quanted_rhs_max >= 3),
                is_signif = (p_value <= p_value_threshold) & (abs(median - offset) >= median_threshold),
                is_hit_nomschecks = is_signif & !is_contaminant,
                is_hit = is_hit_nomschecks & is_valid_comparison, 
                mean_trunc = pmax(-median_max, pmin(median_max, mean - offset)) + offset,
                median_trunc = pmax(-median_max, pmin(median_max, median - offset)) + offset,
                short_label = str_remove_all(object_label, "Oxidation_|Phospho_|_M\\d$|_M\\d\\.\\.\\.$"),
                change = if_else(is_signif, if_else(median < 0, "-", "+"), "."))

#include fp information
require(rlang)
fp.env <- new_environment()
load(file.path(results_path, str_c(project_id, "_msglm_fit_", "fp_20221104", ".RData")), envir = fp.env)
load(file.path(scratch_path, str_c(project_id, "_msdata_full_", "fp_20221104", ".RData")), envir = fp.env)
load(file.path(scratch_path, str_c(project_id, "_msglm_data_", "fp_20221104", ".RData")), envir = fp.env)
fp.env$obj_labu_shift <- fp.env$msdata$pepmodstate_mscalib$zShift

fp_object_contrasts.df <- dplyr::select(fp.env$object_contrasts.df,
                                        fp_protregroup_id=protregroup_id, ci_target,
                                        contrast, #treatment_lhs, treatment_rhs, timepoint_lhs, timepoint_rhs,
                                        fp_median=median, fp_p_value=p_value,
                                        fp_is_hit=is_hit, 
                                        fp_change=change) 

ptmngroup2protregroup.df <- dplyr::select(msdata_full$ptm2gene, ptm_id, protein_ac, ptm_pos) %>%
  dplyr::left_join(fp.env$msdata_full$protein2protregroup) %>% #note: here all the detected protein(re)groups including the ones from fractionation will be counted, but not all of them had enough peptides to be modelled in the real samples!
  dplyr::left_join(dplyr::select(fp.env$msdata_full$protregroups, protregroup_id, npepmods_unique)) %>%
  dplyr::arrange(desc(npepmods_unique)) %>%
  dplyr::distinct() %>%
  dplyr::arrange(ptm_id, desc(is_majority), desc(npepmods_unique), protein_ac_rank) %>%
  dplyr::group_by(ptm_id, ptm_pos) %>%
  dplyr::filter(row_number() == 1L) %>%
  dplyr::ungroup() %>%
  dplyr::select(ptm_id, ptm_pos, protein_ac, fp_protregroup_id=protregroup_id) %>% 
  dplyr::left_join(dplyr::select(msdata_full$ptmngroups, ptmngroup_id, protein_ac, ptm_pos))

object_contrasts.df <- dplyr::left_join(object_contrasts_nofp.df, dplyr::select(msdata_full$ptmngroups, ptmngroup_id, ptmn_id, ptm_id, ptm_pos, protein_ac)) %>% 
  left_join(ptmngroup2protregroup.df) %>%
  dplyr::left_join(fp_object_contrasts.df) %>%
  dplyr::mutate(is_hit_nofp = is_hit,
                is_hit = is_hit & (is_viral | !coalesce(fp_is_hit, FALSE) | sign(median) != sign(coalesce(fp_median, 0))), # |
                                     #(abs(median - fp_median) >= 2)),
                hit_type = case_when(is_signif & !is_hit_nofp ~ "only sig",
                                     is_hit_nofp & !is_hit ~ "is hit nofp",
                                     is_hit & str_detect(object_label, "Oxidation") ~ "oxidation hit",
                                     is_hit & is_viral ~ "viral hit",
                                     is_hit ~ "hit", 
                                     TRUE ~ "non-hit")) 

object_contrast_stats.df <- dplyr::group_by(object_contrasts.df, contrast, contrast_type, ci_target) %>%
  dplyr::summarise(p_value_001 = quantile(p_value, 0.001),
                   p_value_01 = quantile(p_value, 0.01),
                   p_value_05 = quantile(p_value, 0.05),
                   median_abs_50 = quantile(abs(median[p_value <= 0.1]), 0.5),
                   median_abs_95 = quantile(abs(median[p_value <= 0.1]), 0.95),
                   median_abs_99 = quantile(abs(median[p_value <= 0.1]), 0.99),
                   n_hits = sum(is_hit_nomschecks, na.rm = TRUE),
                   n_plus = sum(change == "+"),
                   n_minus = sum(change == "-")) %>%
  dplyr::ungroup() %>% filter(ci_target == "average")

object_contrasts_4enrichment.df <- object_contrasts.df %>% 
  group_by(ptmngroup_id) %>% 
  filter(any(is_valid_comparison)) %>% 
  ungroup()

object_contrasts_wide.df <- tidyr::pivot_wider(object_contrasts.df,
                                        id_cols = c("ci_target", "object_id", "object_label", "is_viral"),
                                        names_from = "contrast", names_sep = ".",
                                        values_from = c("median", "mean", "sd", "p_value", "is_hit", "is_valid_comparison", "change"))
ggplot(dplyr::filter(object_contrasts_wide.df, ci_target=="average",
                       `is_valid_comparison.MPXV_vs_mock@6h` & `is_valid_comparison.MPXV_vs_mock@12h` #&
                                                    #between(`median.MPXV_vs_mock@12h`, -3, 3) &
                                                    #between(`median.MPXV_vs_mock@6h`, -3, 3))
                                     ),
       aes(x = `median.MPXV_vs_mock@6h`, y = `median.MPXV_vs_mock@12h`)) +
  geom_point(alpha = 0.1, size=0.3) +
  #geom_smooth(method = "lm")+
  #geom_abline(slope = 1, intercept = 0, color = "firebrick") + 
  theme_bw_ast()

rfit_filepath <- file.path(results_path, paste0(project_id, '_msglm_fit_', mstype, "_", fit_version, '.RData'))
results_info <- list(project_id = project_id, data_version = data_version,
                     fit_version = fit_version)
message('Saving full analysis results to ', rfit_filepath, '...')
save(results_info, fit_stats, fit_contrasts, ptmngroup2protregroup.df,
     object_contrasts.df, object_contrasts_4enrichment.df, object_contrasts_wide.df, object_contrasts_thresholds.df,
     file = rfit_filepath)

#generate reports----
require(purrr)

ptm_annots.df <- read_tsv(file.path(data_path, str_c(mstype, "_", data_version), str_c("ptm_extractor_", data_version), str_c("ptm_annots_",data_version, ".txt")))

object_contrasts_report.df <- object_contrasts.df %>%
  filter(str_detect(contrast, "MPXV_vs_mock")) %>% 
  left_join(select(msdata$ptmngroups, ptmngroup_id, ptmn_id, ptmn_ids, ptmns, protein_ac, ptm_pos)) %>%
  left_join(msdata_full$ptm2gene) %>% 
  left_join(ptm_annots.df) %>%
  left_join(select(msdata_full$proteins, protein_ac, protein_name)) %>% 
  select(ptmngroup_id, ptmn_ids, gene_name = genename, protein_ac, protein_name, ptm_type, ptm_AA_seq, ptm_pos, 
         ptmngroup_label, ptmns, short_label, protein_ac,
         flanking_15AAs, ptm_is_known, domain, 
         kinase_gene_names, reg_function, reg_prot_iactions, reg_other_iactions, reg_pubmed_ids, diseases,
         ci_target, contrast,
         is_viral, is_contaminant, median, mean, sd, p_value,
         is_signif, is_valid_comparison, is_hit_nomschecks, is_hit_nofp, 
         is_hit, change) %>%
  tidyr::extract(contrast, c("timepoint"), ".*@(\\d+h)", remove = FALSE) %>%
  mutate(timepoint = factor(timepoint, levels = c("6h", "12h", "24h"))) %>% 
  arrange(timepoint) %>% 
  pivot_wider(c(ci_target, ptmngroup_id, gene_name, protein_name, protein_ac, is_viral, is_contaminant,
                ptm_AA_seq, ptm_pos, 
                ptmngroup_label, ptmns, short_label, flanking_15AAs, ptm_is_known, 
                kinase_gene_names, reg_function, reg_prot_iactions, reg_other_iactions, reg_pubmed_ids, diseases),
              names_from = "timepoint", values_from = c("is_valid_comparison", "change", "is_hit", "is_hit_nofp", 
                                                         "median", "p_value", "sd")) %>%
  replace_na(list(ptm_is_known = FALSE)) %>% 
  dplyr::arrange(gene_name, ptm_pos)

write_tsv(object_contrasts_report.df, file.path(analysis_path, "reports", paste0(project_id, '_phospho_contrasts_report_', fit_version, '_wide.txt')))

object_contrasts_report_long.df <- object_contrasts.df %>%
  left_join(select(msdata$ptmngroups, ptmngroup_id, ptmn_id, ptmn_ids, ptmns, protein_ac, ptm_pos)) %>%
  left_join(msdata_full$ptm2gene) %>% 
  left_join(ptm_annots.df) %>%
  left_join(select(msdata_full$proteins, protein_ac, protein_name)) %>% 
  select(ptmngroup_id, ptmn_ids, gene_name = genename, protein_ac, protein_name, ptm_type, ptm_AA_seq, ptm_pos, 
         ptmngroup_label, ptmns, short_label, protein_ac,
         flanking_15AAs, ptm_is_known, domain, 
         kinase_gene_names, reg_function, reg_prot_iactions, reg_other_iactions, reg_pubmed_ids, diseases,
         ci_target, contrast,
         is_viral, is_contaminant, median, mean, sd, p_value,
         is_valid_comparison, is_signif, is_hit_nomschecks, is_hit_nofp, 
         is_hit, change)

write_tsv(object_contrasts_report_long.df, file.path(analysis_path, "reports", paste0(project_id, '_phospho_contrasts_report_', fit_version, '_long.txt')))

#make supplementary table for publication----
object_contrast_average.df <- object_contrasts.df %>% 
  left_join(select(msdata$ptmngroups, ptmngroup_id, ptmn_id, ptmn_ids, ptmns, protein_ac, ptm_pos)) %>%
  left_join(msdata_full$ptm2gene) %>% 
  left_join(ptm_annots.df) %>%
  left_join(select(msdata_full$proteins, protein_ac, protein_name)) %>% 
  filter(ci_target == "average",
         treatment_rhs == "mock",
         ptm_code != "gl"|is.na(ptm_code)) %>% 
  mutate(ptm_AA = ifelse(is.na(ptm_AA), str_extract(object_label, "(?<=\\_)[STY](?=\\d)"),ptm_AA), #ptm_AA is NA if the ptm is unreported before 
         time = paste0(timepoint_lhs, "h"),
         time = factor(time, levels = c("6h", "12h", "24h")),
         change = ifelse(is_hit_nofp, change, "."),
         ptm_site = paste0(ptm_AA, ptm_pos),
         fp_protregroup_id = ifelse(is.na(fp_median), NA, fp_protregroup_id)) %>% 
  select(time, ptmngroup_id, ptmngroup_label, ptm_site, genename, protein_ac, protein_description = protein_name, is_contaminant, is_viral, ptm_pos, ptm_AA, flanking_15AAs,
         other_ptmn_labels = ptmns, protein_id = fp_protregroup_id, short_label,
         is_signif, is_hit, change, protein_change = fp_change, fold_change_log2 = median, protein_fold_change_log2 = fp_median,
         p_value, protein_p_value = fp_p_value, sd_log2 = sd) %>% 
  separate(other_ptmn_labels, sep = ";", into = c("first", "other_ptmn_labels"), extra = "merge") %>% 
  select(-first) %>% 
  arrange(time)

#for the final suppl. table, only include the multiplicity that contains the highest no. of hit (nofp), and lowest multiplicity
ptmn_fit_stats.df <- object_contrast_average.df %>% 
    dplyr::group_by(ptmngroup_id, short_label) %>%
    dplyr::summarise(n_hits = sum(is_signif), .groups="drop") %>%
    dplyr::group_by(short_label) %>%
  dplyr::arrange(desc(n_hits), ptmngroup_id) %>%
    dplyr::mutate(ptmnid_is_reference = row_number() == 1L) %>%
    dplyr::ungroup()

pivoted <- pivot_wider(object_contrast_average.df, ptmngroup_id:protein_id,
                       names_from = "time", values_from = c("is_hit", "change", "protein_change", "fold_change_log2", "protein_fold_change_log2",
                                                            "p_value", "protein_p_value", "sd_log2"),
                       names_sep=".")

names_to_order <- map(unique(object_contrast_average.df$time), ~ names(pivoted)[grep(paste0(".", .x), names(pivoted))]) %>% unlist
names_id <- setdiff(names(pivoted), names_to_order)

object_contrast_4paper.df <- ptmn_fit_stats.df %>% 
  filter(ptmnid_is_reference) %>% 
  select(ptmngroup_id) %>% 
  left_join(pivoted %>% select(names_id, names_to_order))%>% 
  mutate(flanking_15AAs = str_replace(flanking_15AAs, pattern = "(.{16})(.{15})", replacement = "\\1\\*\\2")) %>% 
  arrange(ptmngroup_id)
  
write_tsv(object_contrast_4paper.df,
          file.path(analysis_path, "reports", "sup_tables", paste0("Supplementary table X Phosphoproteome of HFF cells infected with MPXV.txt")))

#prepare a separate table for viral protein only----
ptm_pvalue_ident_max <- 1E-3
ptm_pvalue_quant_max <- 1E-2
ptm_locprob_ident_min <- 0.75
ptm_locprob_quant_min <- 0.5

ptm_status_levels <- c("N/A", "potential", "low conf.", "observed")

ptmngroupXcondition_stats.df <- dplyr::select(msdata_full$ptmngroups, is_viral, ptm_type, ptmngroup_id, ptmn_id) %>% 
  dplyr::filter(is_viral, ptm_type == "Phospho") %>%
  dplyr::left_join(msdata_full$ptmn2pepmodstate) %>%
  dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
  dplyr::left_join(dplyr::select(msdata_full$ptmn_locprobs, evidence_id, ptmn_id, ptm_id, ptm_locprob, msrun, ptm_AA_seq)) %>%
  dplyr::left_join(msdata_full$msruns) %>%
  replace_na(list(condition = "fractionation")) %>% 
  dplyr::group_by(ptm_type,ptm_id, ptmn_id, ptmngroup_id, condition) %>%
  dplyr::summarise(ptm_pvalue_min = min(psm_pvalue, 1, na.rm=TRUE),
                   ptm_locprob_max = max(ptm_locprob, 0, na.rm=TRUE),
                   n_idented_and_localized  = n_distinct(str_c(msrun, " ", pepmodstate_id)[(coalesce(psm_pvalue, 1) <= ptm_pvalue_ident_max) &
                                                                                             (coalesce(ptm_locprob, 0) >= ptm_locprob_ident_min)]),
                   n_quanted  = n_distinct(str_c(msrun, " ", pepmodstate_id)[(coalesce(psm_pvalue, 1) <= ptm_pvalue_quant_max) &
                                                                               (coalesce(ptm_locprob, 0) >= ptm_locprob_quant_min)]),
                   .groups="drop") %>% 
  left_join(msdata_full$ptm2protein) %>% 
  left_join(select(msdata_full$proteins, protein_ac, protein_name, seq)) %>% 
  mutate(seq_length = str_length(seq)) %>% 
  #filter(!is.na(ptm_AA_seq)) %>% 
  dplyr::mutate(ptm_correct = (ptm_type == "Phospho" & ptm_AA_seq %in% c("S","T","Y")),
                is_observed = coalesce(n_idented_and_localized, 0) > 0,
                is_observed_lowconf = (coalesce(ptm_pvalue_min, 1) <= ptm_pvalue_ident_max) |
                                         (coalesce(ptm_locprob_max, 0) >= ptm_locprob_ident_min),
                ptm_status = case_when(is_observed ~ "observed",
                                       is_observed_lowconf ~ "low conf.",
                                       ptm_correct ~ "potential",
                                       TRUE ~ "N/A") %>% factor(ordered=TRUE, levels=ptm_status_levels)) %>% 
  left_join(msdata_full$ptm2gene)

pre_viral_phospho4paper.df <- ptmngroupXcondition_stats.df %>% 
  select(gene_name = genename, protein_ac, protein_description = protein_name, ptm_pos, ptm_AA = ptm_AA_seq, condition, 
         pms_p_value_min = ptm_pvalue_min, ptm_locprob_max,
         n_idented_and_localized, n_quanted,
         ptm_status, flanking_15AAs) %>% 
  mutate(flanking_15AAs = str_replace(flanking_15AAs, pattern = "(.{16})(.{15})", replacement = "\\1\\*\\2")) %>% 
  separate(condition, into = c("source", "time")) %>% 
  mutate(time = as.numeric(time)) %>% 
  arrange(gene_name, ptm_pos, time)

viral_phospho4paper.df <- pre_viral_phospho4paper.df %>% 
  group_by(gene_name, ptm_pos, ptm_AA) %>% 
  mutate(count = n_distinct(source)) %>%
  ungroup() %>% 
  filter(count == 1|(count > 1 & source == "MPXV")) %>% 
  group_by(gene_name, ptm_pos, ptm_AA, source, time) %>% 
  arrange(desc(ptm_locprob_max), pms_p_value_min) %>% #when there are several multiplicities, select the best identified stats
  slice_head() %>% 
  select(-count)

write_tsv(viral_phospho4paper.df,
          file.path(analysis_path, "reports", "sup_tables", paste0("Supplementary table X_2 Phosphoproteome of HFF cells infected with MPXV.txt")))

  
