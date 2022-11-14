# Prepare data for the phosphoproteome analysis for HFF cells infected by monkeypox
# Experiments done in Aug 2022 (new data generated in Nov 2022)
# Author: Yiqi Huang
###############################################################################

project_id <- 'mpxv'
message('Project ID=', project_id)
data_version <- "20221105"
fit_version <- "20221111"
mstype <- "phospho"
message('Dataset version is ', data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, "ggplot_ext.R"))

#library(maxquantUtils)
library(msimportr)
library(msglm)
library(tidyverse)
library(jsonlite)
library(pheatmap)

msfolder <- str_c(mstype, "_", data_version)
msdata_path <- file.path(data_path, msfolder, str_c("ptm_extractor_", data_version))

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version, mstype = mstype, msfolder = msfolder,
                  pepmodstate_mscalib_filename = "mscalib_eclipse_intensity_pepmodstate_eclipse_calib_20220812.json", 
                  quant_type = "intensity", quant_col_prefix = "intensity",
                  pvalue_max=1, pvalue_ident_max=1E-3,
                  locprob_min=0.5, locprob_ident_min=0.75,
                  empty_observation_sigmoid_scale=1/3,
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$pepmodstate_mscalib_filename, '...')
pepmodstate_mscalib <- read_mscalib_json(file.path(data_path, data_info$pepmodstate_mscalib_filename)) 
#the newer calib files are already log transformed. if you see a warning later about log transformation it means you're using an old file, which is not ideal!

fasta.dfs <- list(
  human = read_innate_uniprot_fasta(file.path(data_path, msfolder, "fasta/uniprot-reviewed_yes+taxonomy_9606.fasta")),
  mpxv = read_innate_uniprot_fasta(file.path(data_path, msfolder, "fasta/MPXV_reformatted_20220811.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, msfolder, "fasta/contaminants_20200405.fasta"))
)

msdata_full <- lapply(list(
  msruns = "rawfiles_info.txt",
  proteins = "proteins.txt.gz",
  ptmn_locprobs = "ptmn_locprobs.txt.gz",
  pepmodstate_intensities = "pms_intensities.txt.gz",
  #protgroups = "protgroups.txt.gz",
  #peptide2protgroup = "peptide_to_protgroup.txt.gz",
  peptides = "peptides.txt.gz",
  pepmodstates = "pepmodstates.txt.gz",
  #protein2protgroup = "protein_to_protgroup.txt.gz",
  peptide2protein = "peptide_to_protein.txt.gz",
  ptm2protein = "ptm_to_protein.txt.gz",
  ptm2gene = "ptm_to_gene.txt.gz",
  ptmns = "ptmns.txt.gz",
  ptmn2pepmodstate = "ptmn_to_pepmodstate.txt.gz",
  ptmn2ptmngroup = "ptmn_to_ptmngroup.txt.gz",
  #ptmngroups = "ptmngroups.txt.gz",
  ptmngroup2pepmodstate = "ptmngroup2pms.txt.gz"),
  function(fname) {
    message("Reading ", fname, "...")
    read_tsv(file.path(msdata_path, fname))
  }
)

msdata_full$msruns <- msdata_full$msruns %>%
  dplyr::mutate(is_ptm = replace_na(is_ptm, TRUE),
                treatment = factor(treatment, c("mock", "MPXV")),
                timepoint_num = as.integer(as.character(timepoint)),
                timepoint = factor(timepoint),
                replicate = as.integer(replicate)
                ) %>% 
  arrange(msrun)

msdata_full$pepmodstate_intensities <- dplyr::left_join(msdata_full$pepmodstate_intensities, dplyr::select(msdata_full$msruns, msrun, msexperiment)) %>% 
  dplyr::mutate(ident_type = factor(ident_type, c("ISO-MSMS", "MULTI-MSMS", "MSMS", "MULTI-SECPEP", "MULTI-MATCH", "MULTI-MATCH-MSMS"), ordered = TRUE),
                is_idented = coalesce(psm_pvalue, 1) <= data_info$pvalue_ident_max)

msdata_full$ptmn_locprobs <- dplyr::inner_join(msdata_full$ptmn_locprobs,
                                               dplyr::select(msdata_full$msruns, msrun, msexperiment))

msdata_full$ptmngroup2pepmodstate <- msdata_full$ptmngroup2pepmodstate %>% 
  mutate(is_specific = (ptm_type != "Oxidation")) %>% 
  group_by(ptmngroup_id, ptm_type) %>% 
  mutate(ptmngroup_label = dplyr::first(ptmn_label, order_by = ptmn_id),
         ptmngroup_label = ifelse(n_distinct(ptmn_label)>1, str_c(ptmngroup_label, "..."), ptmngroup_label)) %>% 
  dplyr::ungroup() #exclude all oxidation sites here
  
# add more properties to ptmngroups as it's the object being modeled, note that here the n_ident refers to the peaks not pepmodstates. There could be multiple peaks to one pepmodstate per msrun.
msdata_full$ptmngroups <- msdata_full$ptmn2ptmngroup %>%
  dplyr::left_join(dplyr::select(msdata_full$ptmns, ptmn_id, ptm_id, ptmn_label)) %>% 
  dplyr::inner_join(dplyr::select(dplyr::filter(msdata_full$ptm2gene, ptm_is_reference), ptm_id, ptm_pos, protein_ac, is_viral, is_contaminant)) %>%
  dplyr::left_join(dplyr::select(msdata_full$proteins, gene_name=genename, protein_ac)) %>%
  mutate(ptmn_label_short = str_remove_all(ptmn_label, "(^[^_]+_)|(_M\\d+$)")) %>% 
  group_by(ptmngroup_id) %>% 
  mutate(ptmn_ids = paste(ptmn_id, collapse = ";"),
         ptm_ids = paste(ptm_id, collapse = ";"),
         ptm_pos_all = paste(ptm_pos, collapse = ";"),
         protein_acs = paste(protein_ac, collapse = ";"),
         gene_names = paste(gene_name, collapse = ";"),
         ptmns = paste(ptmn_label_short, collapse = ";")
         ) %>% 
  left_join(msdata_full$ptmngroup2pepmodstate %>% select(ptmngroup_id, ptmngroup_label)) %>% 
  select(ptm_type, nselptms,ptmngroup_id, ptmngroup_label, ptmn_id, ptm_pos, protein_ac, gene_name, is_viral, is_contaminant, ptmn_ids, ptm_ids, ptm_pos_all, protein_acs, ptmns ) %>% 
  slice_head() %>% ungroup() %>% distinct() %>% 
  dplyr::mutate(ptmngroup_label_no_ptm_type = str_remove(ptmngroup_label, "^[^_]+_"))

msdata_full$ptmngroup_idents <- inner_join(msdata_full$ptmngroups, dplyr::select(msdata_full$ptmngroup2pepmodstate, ptmngroup_id, pepmodstate_id) %>% distinct()) %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmodstate_intensities, evidence_id, pepmodstate_id, msrun, msexperiment, intensity, ident_type, psm_pvalue)) %>%
  dplyr::inner_join(dplyr::select(msdata_full$ptmn_locprobs, evidence_id, ptm_locprob, ptmn_id)) %>%
  dplyr::mutate(is_quanted = !is.na(intensity),
                is_idented =  coalesce(ident_type, "") %in% c("ISO-MSMS", "MULTI-MSMS", "MSMS", "MULTI-SECPEP") &
                  coalesce(psm_pvalue, 1) <= data_info$pvalue_ident_max,
                is_localized = coalesce(ptm_locprob, 0) >= data_info$locprob_ident_min,
                is_valid_quant = is_quanted & coalesce(ptm_locprob, 0) >= data_info$locprob_min &
                  coalesce(psm_pvalue, 1) <= data_info$pvalue_max) %>%
  dplyr::group_by(msexperiment, ptmngroup_id, ptmngroup_label, ptm_type) %>%
  dplyr::summarise(n_quanted = sum(is_quanted),
                   n_idented = sum(is_idented),
                   #n_msruns = n_distinct(msrun),
                    n_pepmodstates = n_distinct(pepmodstate_id),
                   n_localized = sum(is_localized),
                   n_idented_and_localized = sum(is_localized & is_idented),
                   n_valid_quants = sum(is_valid_quant),
                   ptm_pvalue_min = min(psm_pvalue, na.rm=TRUE),
                   ptm_locprob_max = max(ptm_locprob, na.rm=TRUE),
                   ident_type = case_when(any(str_detect(ident_type, "MATCH")) ~ "By matching",
                                          any(str_detect(ident_type, "MS|SEC")) ~ "By MS/MS",
                                          TRUE ~ NA_character_),
                   .groups="drop")

all_proteins.df <- dplyr::bind_rows(
  dplyr::mutate(fasta.dfs$mpxv, is_viral = TRUE, is_contaminant = FALSE),
  dplyr::mutate(fasta.dfs$human, is_viral = FALSE, is_contaminant=FALSE),
  dplyr::mutate(fasta.dfs$contaminants, is_viral = FALSE, is_contaminant=TRUE))

msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

# condition = treatment X timepoint
conditions.df <- dplyr::select(msdata_full$msruns, condition, treatment, timepoint, timepoint_num) %>%
  dplyr::distinct()
conditions.df <- dplyr::left_join(conditions.df,
                                  dplyr::select(conditions.df, treatment, timepoint_after=timepoint_num)) %>%
  dplyr::mutate(is_after = timepoint_num >= timepoint_after) %>%
  tidyr::pivot_wider(all_of(c("condition", "treatment", "timepoint", "timepoint_num")),
                     names_prefix = "after", names_from = timepoint_after, values_from = is_after) %>%
  rename_at(vars(starts_with("after")), ~str_c(., "h")) %>%
  mutate(after0h = NULL, # not needed for exp_design
         infected = treatment != "mock")

# setup experimental design matrices
conditionXeffect_orig.mtx <- model.matrix(
  formula(str_c("~ (" , str_c(str_subset(colnames(conditions.df), "^after\\d+h"), collapse =" + "), " ) * treatment")),
  conditions.df
)
conditionXeffect.mtx <- conditionXeffect_orig.mtx[, str_detect(colnames(conditionXeffect_orig.mtx),
                                                               "^treatment[^:]+$", negate=TRUE) &
                                                    (apply(conditionXeffect_orig.mtx, 2, function(x) min(abs(x))) == 0)] 

dimnames(conditionXeffect.mtx) <- list(condition = as.character(conditions.df$condition),
                                       effect = colnames(conditionXeffect.mtx))
all_conditions <- as.character(conditions.df$condition)

plots_path <- file.path(analysis_path, "plots", str_c(msfolder, "_", fit_version))
if (!dir.exists(plots_path)) dir.create(plots_path)

(pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(plots_path,
                              paste0(project_id, "_exp_design_",fit_version, ".pdf")),
         width = 6, height = 6))

effects.df <- tibble(effect = colnames(conditionXeffect.mtx)) %>% 
  mutate(treatment = effect_factor(effect, "treatment", levels(conditions.df$treatment), NA),
         timepoint = effect_factor(effect, 'after', str_c(unique(conditions.df$timepoint_num), "hTRUE"), NA) %>%
           as.character() %>% str_remove_all("^after|hTRUE$") %>% factor(levels = levels(conditions.df$timepoint)),
         
         effect_type = case_when(!is.na(treatment) & !is.na(timepoint) ~ "treatmentXtimepoint",
                                 !is.na(timepoint) ~ "timepoint",
                                 TRUE ~ NA_character_),
         effect_label =  case_when(effect_type == "treatmentXtimepoint" ~ str_c(treatment, "+", timepoint, "h"),
                                   effect_type == "timepoint" ~ str_c(timepoint, "h"),
                                   TRUE ~ NA_character_),

         prior_mean = 0, #We don't expect any general trends due to any treatment, so the mean is 0 here.
         prior_tau = case_when(effect_type == "treatmentXtimepoint" ~ 0.25,
                               effect_type == "timepoint" ~ 0.5,
                               TRUE ~ NA_real_), #This part belongs to the horseshoe prior. Tau represents expected none-zero values. Controls how much the model would shrink the "low values". The higher it is, the more none-zero effect we expect and the less it will shrink!!! 
         prior_df1 = case_when(TRUE ~2.0), 
         prior_df2 = case_when(TRUE ~ 2.0), #This part controls the horseshoe+ prior. Controls how conservative the model is. We can leave everything at 2.0 (default value) or just delete it, and adjust it later if we see problems in the outcome.
         is_positive = FALSE)

#The metaconditions in this case correspond to the conditions
all_metaconditions <- c(all_conditions)
conditionXmetacondition.mtx <- msglm::constant_matrix(FALSE, list(condition = all_conditions,
                                                                  metacondition = all_metaconditions))

for (cname in unique(conditions.df$condition)) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
}

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>% 
  filter(n != 0) %>% 
  select(-n)

(pheatmap(ifelse(conditionXmetacondition.mtx, 1L, 0L), cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plots_path,
                              paste0(project_id, "_metaconditions_",fit_version, ".pdf")),
         width = 8, height = 6))

contrasts.df <- bind_rows(
  # all treatment pairs at each timepoint
  left_join(dplyr::mutate(dplyr::select(conditions.df, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint)),
            dplyr::mutate(dplyr::select(conditions.df, metacondition_rhs = condition, treatment_rhs = treatment, timepoint_lhs = timepoint))) %>%
    filter((as.integer(treatment_lhs) > as.integer(treatment_rhs)) &
             !((treatment_lhs == "infected") & (treatment_rhs != "mock"))) %>%
    mutate(timepoint_rhs = timepoint_lhs,
           contrast_kind = "treatment_vs_treatment",
           contrast = str_c(treatment_lhs, "_vs_", treatment_rhs, "@", timepoint_rhs, "h")),
  # all timepoints of the same treatment
  left_join(dplyr::mutate(dplyr::select(conditions.df, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint)),
            dplyr::mutate(dplyr::select(conditions.df, metacondition_rhs = condition, treatment_lhs = treatment, timepoint_rhs = timepoint))) %>%
    filter(as.integer(timepoint_lhs) > as.integer(timepoint_rhs),
           treatment_lhs != "mock") %>%
    mutate(treatment_rhs = treatment_lhs,
           contrast_kind = "timepoint_vs_timepoint",
           contrast = str_c(timepoint_lhs, "h_vs_", timepoint_rhs, "h@", treatment_lhs))
) %>%
  mutate(contrast_type = "comparison",
         offset = 0.0)
all_contrasts = contrasts.df$contrast

metaconditionXcontrast.mtx <- msglm::constant_matrix(0, list(metacondition = all_metaconditions,
                                                             contrast = all_contrasts))
for (i in 1:nrow(contrasts.df)) {
  contr <- contrasts.df$contrast[[i]]
  metaconditionXcontrast.mtx[contrasts.df$metacondition_lhs[[i]], contr] <- 1.0
  metaconditionXcontrast.mtx[contrasts.df$metacondition_rhs[[i]], contr] <- -1.0
} #do not rearrange the order in effects, otherwise it will create mismatch here

(pheatmap(metaconditionXcontrast.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(plots_path,
                              paste0(project_id, "_exp_design_contrasts_",fit_version, ".pdf")),
         width = 8, height = 12))

metaconditionXcontrast.df <- as_tibble(as.table(metaconditionXcontrast.mtx)) %>% 
  filter(n != 0) %>% 
  dplyr::rename(weight = n) %>% 
  mutate(contrast_type = "comparison",
         condition_role = if_else(contrast_type == 'filtering',
                                  if_else(weight > 0, 'signal', 'background'),
                                  'signal'))

contrasts.df <- dplyr::select(metaconditionXcontrast.df, contrast, contrast_type) %>%
  dplyr::distinct()

#no batch effects, everything done at once

msglm_def <- msglm_model(conditionXeffect.mtx, conditions.df, effects.df,
                         verbose=TRUE) %>%
  msglm::set_contrasts(metaconditionXcontrast.mtx, conditionXmetacondition.mtx,
                       contrasts.df)

msdata <- import_msglm_data(msdata_full, msglm_def, object="ptmngroup", quantobject="pepmodstate", verbose = TRUE) 

##########
## normalization
msdata4norm.df <- msdata_full$pepmodstate_intensities %>% ungroup() %>%
 #dplyr::filter(!is.na(intensity) & coalesce(psm_pvalue, 1) <= 1E-4) %>% # this filter is less relevant because we have effect regularization now 
  dplyr::inner_join(dplyr::filter(msdata_full$pepmodstates, str_detect(pepmod_seq, "Phospho|GlyGly")) %>% # select phospho or ubi pepmods, 
                      dplyr::select(pepmodstate_id, peptide_id)) %>%
  dplyr::semi_join(dplyr::filter(msdata_full$peptides, !is_contaminant & !is_reverse)) %>%
  dplyr::semi_join(dplyr::filter(msdata_full$proteins, !is_contaminant & !is_viral) %>%
                     dplyr::inner_join(msdata_full$peptide2protein)) %>%
  dplyr::group_by(pepmodstate_id, msrun) %>%
  dplyr::summarise(intensity = sum(intensity, na.rm = TRUE), .groups="drop") %>%
  dplyr::mutate(intensity = if_else(intensity != 0, intensity, NA_real_))

require(cmdstanr)
options(mc.cores=8L)
Sys.setenv(MKL_NUM_THREADS=1L)

# normalise msruns within the group of msruns with the same condition and timepoint
msruns_hnorm <- multilevel_normalize_experiments(pepmodstate_mscalib,
                                                 mutate(msdata$msruns, condition = condition, 
                                                        msdata$msruns, timepoint = timepoint),
                                                 msdata4norm.df,
                                                 quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "msrun",
                                                 #mcmc.iter = 2000L,
                                                 verbose=TRUE,
                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=500L),
                                                                    condition = list(cond_col="condition", max_objs=500L, missing_exp.ratio=0.2),
                                                                    timepoint = list(cond_col="timepoint", max_objs=500L, missing_exp.ratio=0.3)
                                                 ))

msdata$msrun_shifts <- msruns_hnorm$mschannel_shifts

rmsglmdata_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_data_', data_info$mstype, "_", data_info$fit_ver, '.RData'))
message('Saving data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata, msglm_def, msruns_hnorm,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_',  data_info$mstype, "_", data_info$data_ver,  '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full, fasta.dfs,
     file = rfulldata_filepath)

message('Done.')
