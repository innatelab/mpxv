#assembly and analysis of iVIP+MPXV infection in HeLa FP dataset
#msglm run by Sabri
#Author: Yiqi Huang
###############################################################################

project_id <- 'mpxv'
message('Project ID=', project_id)
data_version_ivip <- "20221129"
fit_version_ivip <- "20221130"

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))
require(msglm)
library(tidyverse)
library(rlang)

message('Loading data...')
load(file.path(scratch_path, paste0(project_id, '_ivip_msglm_data_', data_version_ivip, '.RData')))
load(file.path(scratch_path, paste0(project_id, '_ivip_msdata_full_',  data_version_ivip, '.RData')))
load(file.path(results_path, paste0(project_id, '_ivip_msglm_fit_', fit_version_ivip, '.RData')))

modelobj <- msdata$msentities[["object"]]
quantobj <- msdata$msentities[["quantobject"]]
modelobjs_df <- msdata$objects
modelobj_idcol <- paste0(modelobj, "_id")
quantobj_idcol <- paste0(quantobj, "_id")

#correct the msrun information (first two lines only used for fit ver 20220912)----
msruns.df <- msdata$msruns %>% 
  mutate(biological_rep = str_extract(msexperiment, "(?<=\\dh_)\\d(?=_)"))

msruns_summary <- msdata$msruns %>% group_by(condition, treatment, timepoint, replicate) %>% 
  summarise(count = n()) %>% 
  group_by(condition) %>% 
  mutate(total_msruns = sum(count))

ggplot(msruns_summary %>% select(condition, treatment, timepoint, total_msruns) %>% unique(), 
       aes(x = factor(timepoint), y = total_msruns)) + 
  geom_bar(stat = "identity")+
  theme_classic()+
  facet_wrap(~ treatment)

#generate reports---- 
if (modelobj == "protgroup") {
  # FIXME stats should be quantobj-dependent
  obj_conditions.df <- tidyr::expand(msdata_full$protgroup_intensities, protgroup_id, msrun) %>%
    dplyr::left_join(msdata_full$protgroup_intensities) %>%
    dplyr::inner_join(msdata$msruns) %>%
    dplyr::mutate(is_quanted = !is.na(intensity),
                  is_idented = replace_na(ident_type == "By MS/MS", FALSE)) %>%
    dplyr::group_by(condition, protgroup_id) %>%
    dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_quanted]),
                     nmsruns_idented = n_distinct(msrun[is_idented])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(object_id = protgroup_id)
  
  modelobjs_df <- dplyr::mutate(modelobjs_df,
                                is_msvalid_object = npepmods_unique >= 2L)#(nprotgroups_sharing_proteins == 1 || nproteins_have_razor > 0))
} else if (modelobj == "protregroup") {
  obj_conditions.df <- expand(msdata_full$pepmodstate_intensities, msrun, pepmod_id) %>%
    dplyr::inner_join(msdata$msruns) %>%
    dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
    dplyr::mutate(is_quanted = !is.na(intensity),
                  is_idented = is_quanted) %>%
    dplyr::left_join(filter(msdata$protregroup2pepmod, is_specific)) %>%
    dplyr::group_by(condition, protregroup_id) %>%
    dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_idented]), # count msruns
                     nmsruns_idented = n_distinct(msrun[is_quanted])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(object_id = protregroup_id)
  
  modelobjs_df <- dplyr::mutate(modelobjs_df,
                                is_msvalid_object = npepmods_unique >= 2L)
}

contrastXcondition.df <- as_tibble(as.table(msglm_def$conditionXmetacondition)) %>% dplyr::filter(n != 0) %>%
  dplyr::select(-n) %>%
  dplyr::inner_join(as_tibble(as.table(msglm_def$metaconditionXcontrast))) %>% dplyr::filter(n != 0) %>% 
  dplyr::arrange(contrast, metacondition, condition)

pre_object_contrasts.df <- obj_conditions.df %>% 
  group_by(object_id, condition) %>%
  mutate(obs = seq(n())) %>%
  ungroup %>%
  # complete the data
  tidyr::complete(object_id, condition, obs) %>%
  select(-obs) %>%
  mutate(protregroup_id = coalesce(protregroup_id, object_id)) %>% 
  replace_na(list( nmsruns_quanted = 0, nmsruns_idented = 0)) %>% 
  dplyr::inner_join(contrastXcondition.df) %>%
  left_join(select(msruns_summary, condition, total_msruns) %>% unique()) %>% 
  dplyr::mutate(is_lhs = n > 0) %>%
  dplyr::group_by(object_id, contrast, is_lhs, total_msruns
                  ) %>%
  dplyr::summarise(has_quanted = any(!is.na(nmsruns_quanted)),
                   nmsruns_quanted_min = min(nmsruns_quanted, na.rm=TRUE),
                   nmsruns_quanted_max = max(nmsruns_quanted, na.rm=TRUE),
                   has_idented = any(!is.na(nmsruns_idented)),
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
                   total_msruns_lhs = total_msruns[is_lhs],
                   nmsruns_quanted_rhs_min = nmsruns_quanted_min[!is_lhs],
                   nmsruns_quanted_rhs_max = nmsruns_quanted_max[!is_lhs],
                   nmsruns_idented_rhs_min = nmsruns_idented_min[!is_lhs],
                   nmsruns_idented_rhs_max = nmsruns_idented_max[!is_lhs],
                   total_msruns_rhs = total_msruns[!is_lhs]
                   ) %>%
  dplyr::ungroup() # In our case, the condition and metacondition is the same, so the nmsruns_quanted/idented will always have the same min and max, but it's not always the case if there're multiple conditions in one metacondition (e.g. in AP-MS analysis)

contrasts.df <- dplyr::ungroup(msglm_def$contrasts) %>%
  dplyr::inner_join(tidyr::pivot_wider(dplyr::mutate(as.data.frame.table(msglm_def$metaconditionXcontrast, responseName="w"),
                                                     side = if_else(w > 0, "lhs", "rhs")) %>% dplyr::filter(w != 0),
                                       c(contrast), names_from = "side", values_from = "metacondition",
                                       names_glue = "{.value}_{side}")) %>%
  dplyr::mutate(offset = 0, offset_prior = 0) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint)) %>%
  dplyr::left_join(dplyr::select(msglm_def$conditions, metacondition_rhs = condition, treatment_rhs = treatment, timepoint_rhs = timepoint))

object_contrasts_thresholds.df <- contrasts.df %>%
  mutate(p_value_threshold = case_when(TRUE ~ 0.001),
         #p_value_threshold_lesser = case_when(TRUE ~ 0.01),
         median_threshold = case_when(TRUE ~ 0.5)#,
         #median_threshold_lesser = case_when(TRUE ~ 0.125)
  )

object_contrasts.df <- fit_contrasts$object_conditions %>%
  inner_join(pre_object_contrasts.df) %>% 
  dplyr::filter(var %in% c("obj_cond_labu", "obj_cond_labu_replCI")) %>%
  dplyr::inner_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(#is_valid_comparison = (pmax(nmsruns_quanted_lhs_max, nmsruns_quanted_rhs_max)>=3), #quantified in 3/5 on either side
                is_valid_comparison = (nmsruns_quanted_lhs_max >= 0.5*total_msruns_lhs)|(nmsruns_quanted_rhs_max >= 0.5*total_msruns_rhs),
                is_signif = p_value <= p_value_threshold & abs(median) >= median_threshold,
                #is_signif_lesser = p_value <= p_value_threshold_lesser & abs(median) >= median_threshold_lesser,
                is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse,
                is_hit = is_hit_nomschecks & is_valid_comparison, 
                change = if_else(is_signif, if_else(median < 0, "-", "+"), ".")) 

object_contrast_stats.df <- dplyr::group_by(object_contrasts.df, contrast, treatment_lhs, timepoint_lhs, contrast_type, ci_target) %>%
  dplyr::summarise(p_value_001 = quantile(p_value, 0.001),
                   p_value_01 = quantile(p_value, 0.01),
                   p_value_05 = quantile(p_value, 0.05),
                   median_abs_50 = quantile(abs(median[p_value <= 0.1]), 0.5),
                   median_abs_95 = quantile(abs(median[p_value <= 0.1]), 0.95),
                   median_abs_99 = quantile(abs(median[p_value <= 0.1]), 0.99),
                   n_hits = sum(is_hit, na.rm = TRUE),
                   n_hits_nomschecks = sum(is_hit_nomschecks, na.rm = TRUE),
                   n_plus = sum(change == "+"),
                   n_minus = sum(change == "-")) %>%
  dplyr::ungroup() %>% 
  filter(ci_target == "average") %>% 
  arrange(treatment_lhs, timepoint_lhs)#%>% select(contrast, n_hits, n_plus, n_minus)

View(filter(object_contrast_stats.df, ci_target == "average") %>% dplyr::arrange(desc(n_hits)))

object_contrast_stats.df %>% select(contrast, treatment_lhs, timepoint_lhs, n_plus, n_minus) %>% 
  pivot_longer(cols = c(n_plus, n_minus), names_to = "change", values_to = "n") %>% 
  mutate(plot_n = ifelse(change == "n_plus", n, -n),
         change = relevel(factor(change), ref = "n_plus")) %>% 
  #arrange(timepoint_lhs) %>% View()
  ggplot(aes(timepoint_lhs, plot_n, fill = change))+ 
  geom_bar(stat = "identity")+
  facet_wrap(~treatment_lhs)

object_contrasts.df %>% filter(treatment_lhs == "MVA_dE9L", timepoint_lhs == 0, is_hit, ci_target == "average") %>% 
  ggplot(aes(nmsruns_quanted_lhs_min))+
  geom_bar()

report_cols <- c("object_label", "object_id", "gene_names",
                 "majority_protein_acs", "protein_descriptions",
                 "is_contaminant", "is_viral")

objects4report.df <- dplyr::select(msdata$objects, any_of(report_cols)) %>%
  dplyr::semi_join(dplyr::select(object_contrasts.df, object_id))

object_contrasts_report.df <- objects4report.df %>%
  dplyr::left_join(pivot_wider(object_contrasts.df, c(ci_target, object_id),
                               names_from = "contrast", values_from = c("change", "is_valid_comparison", "is_hit",  "median", "p_value", "sd"),
                               names_sep=".")) %>%
  dplyr::arrange(gene_names, majority_protein_acs, ci_target)

write_tsv(object_contrasts_report.df,
          file.path(analysis_path, "reports", paste0(project_id, '_ivip_contrasts_report_', fit_version_ivip, '_wide.txt')))

object_contrasts_wide.df <- pivot_wider(object_contrasts.df,
                                        id_cols = c("ci_target", "object_id", "object_label", "is_viral"),
                                        names_from = "contrast", names_sep = ".",
                                        values_from = c("median", "mean", "sd", "p_value", "is_hit", "change"))

rfit_filepath <- file.path(results_path, paste0(project_id, '_ivip_msglm_fit_', fit_version_ivip, '.RData'))
results_info <- list(project_id = project_id, mstype = "fp",
                     data_version = data_version_ivip, fit_version = fit_version_ivip,
                     modelobj = modelobj, quantobj = quantobj)
message('Saving full analysis results to ', rfit_filepath, '...')
save(results_info, fit_stats, fit_contrasts,
     object_contrasts.df, object_contrasts_wide.df,
     object_contrasts_thresholds.df,
     file = rfit_filepath)
message('Done.')

#make the tables for publication----
object_contrast_average.df <- object_contrasts.df %>% 
  filter(ci_target == "average", timepoint_lhs != 0) %>% 
  separate(treatment_lhs, into = c("virus", "other")) %>% 
  mutate(virus = ifelse(other == "dE9L", paste0(virus, "dE3L"), virus),
         virus = factor(virus, levels = c("VACV", "CVA", "MVA", "MVAdE3L")),
         condition = paste0(virus, "@", timepoint_lhs, "h"),
         change = ifelse(is_hit, change, ".")) %>% 
  arrange(timepoint_lhs, virus)

pivoted <- pivot_wider(object_contrast_average.df, object_id,
                        names_from = "condition", values_from = c("is_hit", "change",  "median", "p_value", "sd"),
                        names_sep=".")

names_to_order <- map(unique(object_contrast_average.df$condition), ~ names(pivoted)[grep(paste0(".", .x), names(pivoted))]) %>% unlist
names_id <- setdiff(names(pivoted), names_to_order)

object_contrast_4paper.df <- objects4report.df %>% 
  left_join(pivoted %>% select(names_id, names_to_order)) %>% 
  arrange(object_id)

write_tsv(object_contrast_4paper.df,
          file.path(analysis_path, "reports", "sup_tables", paste0("Supplementary table X Total proteome of HeLa cells infected with different poxviruses.txt")))
