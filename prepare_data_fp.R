# Prepare data for the full proteome analysis for HFF cells infected by monkeypox
# Experiments done in Aug 2022 (new exp done in Nov 2022)
# Author: Yiqi Huang
###############################################################################

project_id <- 'mpxv'
message('Project ID=', project_id)
data_version <- "20221104"
fit_version <- "20221104"
mstype <- "fp"
message('Dataset version is ', data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, "ggplot_ext.R"))

library(msimportr)
library(msglm)
library(tidyverse)
library(jsonlite)
library(pheatmap)

msfolder <- str_c(mstype, "_", data_version)
msdata_path <- file.path(data_path, msfolder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version, mstype = mstype, msfolder = msfolder,
                  pepmodstate_mscalib_filename = "mscalib_eclipse_intensity_pepmodstate_eclipse_calib_20220812.json", 
                  quant_type = "LFQ", quant_col_prefix = "LFQ_Intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$pepmodstate_mscalib_filename, '...')
pepmodstate_mscalib <- read_mscalib_json(file.path(data_path, data_info$pepmodstate_mscalib_filename)) 
#the newer calib files are already log transformed. if you see a warning later about log transformation it means you're using an old file, which is not ideal!

msruns.df <- read_tsv(file.path(msdata_path, "combined/experimentalDesign.txt")) %>% 
  dplyr::rename(rawfile=Name, fraction=Fraction, msexperiment_mq=Experiment, is_ptm = PTM) %>%
  tidyr::extract("msexperiment_mq", c("mstype", "treatment", "timepoint", "replicate"),
                 "HFF_(FP)_([^_]+)_(\\d+)h_(\\d)?", remove = FALSE) %>% 
  dplyr::mutate(mstype = tolower(mstype),
                is_ptm = replace_na(is_ptm, FALSE),
                timepoint_num = as.integer(as.character(timepoint)),
                timepoint = factor(timepoint),
                treatment = str_replace(treatment, "MPX", "MPXV"),
                treatment = factor(treatment, c("mock", "MPXV")),
                condition = str_c(treatment, "_", timepoint),
                replicate = as.integer(replicate),
                msexperiment = str_c(condition, "_", replicate),
                msexperiment = factor(msexperiment, levels=unique(msexperiment)),
                is_skipped = is.na(mstype)
  ) %>% 
  dplyr::arrange(treatment, timepoint_num, replicate)

fasta.dfs <- list(
  human = read_innate_uniprot_fasta(file.path(data_path, msfolder, "fasta/uniprot-reviewed_yes+taxonomy_9606.fasta")),
  mpxv = read_innate_uniprot_fasta(file.path(data_path, msfolder, "fasta/MPXV_reformatted_20220811.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, msfolder, "fasta/contaminants_20200405.fasta"))
)

msdata.wide <- read.MaxQuant.ProteinGroups(file.path(msdata_path, 'combined/txt'), import_data = c(data_info$quant_type, "ident_type"))
msdata_colgroups <- attr(msdata.wide, "column_groups")

mqevidence <- read.MaxQuant.Evidence(file.path(msdata_path, 'combined', 'txt'),
                                     mschannel_annotate.f = function(mschans_df) {
                                       res <- dplyr::inner_join(dplyr::select(mschans_df, mstag, rawfile),
                                                                msruns.df,
                                                                by="rawfile")
                                       attr(res, "column_scopes") <- c(timepoint = "msexperiment",
                                                                       timepoint_num = "msexperiment",
                                                                       treatment = "msexperiment",
                                                                       condition = "msexperiment",
                                                                       replicate = "msexperiment") #include all the columns you want to keep later!
                                       return(res)
                                     })

mqevidence$peptides <- read.MaxQuant.Peptides(file.path(msdata_path, 'combined', 'txt'), file_name='peptides.txt',
                                              import_data='ident_type')
mqevidence$peaks <- NULL # exclude big data frame

strlist_label <- function(strs) {
  str_c(strs[[1]], if_else(n_distinct(strs) > 1, '...', ''))
}
strlist_label2 <- function(strs, delim=fixed(';')) {
  sapply(str_split(strs, delim), strlist_label)
}

# all ms data
msdata_full <- mqevidence[c('msexperiments', 'msruns')]

msdata_full <- append_protgroups_info(msdata_full, msdata.wide,
                                      proteins_info = dplyr::bind_rows(
                                        dplyr::mutate(fasta.dfs$contaminants, is_viral = FALSE, is_contaminant=TRUE),
                                        dplyr::mutate(fasta.dfs$mpxv, is_viral = TRUE, is_contaminant=FALSE),
                                        dplyr::mutate(fasta.dfs$human, is_viral = FALSE, is_contaminant = FALSE)) %>%
                                        dplyr::mutate(protein_ac_noiso = str_remove(protein_ac, "-\\d+(?:#.+)*$"),
                                                      protein_isoform_ix = replace_na(as.integer(str_match(protein_ac, "-(\\d+)$")[, 2]), 1L)),
                                      import_columns = c("is_viral", "is_contaminant", "organism"))

msdata_full$peptides <- mqevidence$peptides %>%
  dplyr::left_join(select(msdata_full$proteins, lead_razor_protein_ac = protein_ac)) %>% 
  mutate(peptide_rank = case_when(is_reverse ~ -1L,
                                  TRUE ~ 1L))

msdata_full$pepmods <- mqevidence$pepmods %>%
  dplyr::left_join(select(msdata_full$peptides, peptide_id)) %>% 
  mutate(peptide_rank = case_when(is_reverse ~ -1L,
                                  TRUE ~ 1L))
# pepmods and peptides could be ranked by their relevance to the experiment, e.g. APMS, compartment... see examples in adhoc
msdata_full$pepmodstates <- mqevidence$pepmodstates

# redefine protein groups (protregroups)
peptides.df <- dplyr::select(msdata_full$peptides, peptide_id, protgroup_ids, protein_acs, lead_razor_protein_ac, peptide_seq, is_reverse, peptide_rank)
proteins.df <- msdata_full$proteins
save(file = file.path(msdata_path, str_c(project_id, "_", msfolder, "_peptides.RData")),
     peptides.df, proteins.df)

# .. run protregroup.jl
msdata_full$protregroups <- read_tsv(file.path(msdata_path, 
                                               str_c(project_id, "_", msfolder, "_protregroups.txt")),
                                     col_types = list(protregroup_id = "i"))

msdata_full$protein2protregroup <- dplyr::select(msdata_full$protregroups, protregroup_id, protein_ac=majority_protein_acs) %>%
  separate_rows(protein_ac, sep=fixed(";"), convert=TRUE) %>%
  dplyr::mutate(is_majority = TRUE) %>%
  dplyr::group_by(protregroup_id) %>%
  dplyr::mutate(protein_ac_rank = row_number()) %>%
  dplyr::ungroup()

msdata_full$protregroup2peptide <- bind_rows(
  select(msdata_full$protregroups, protregroup_id, peptide_id=spec_peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = TRUE),
  select(msdata_full$protregroups, protregroup_id, peptide_id=peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = FALSE)) %>%
  dplyr::group_by(protregroup_id, peptide_id) %>%
  dplyr::summarise(is_specific = any(is_specific)) %>%
  dplyr::ungroup()
msdata_full$protregroup2pepmod <- dplyr::inner_join(msdata_full$protregroup2peptide,
                                                    dplyr::select(msdata_full$pepmods, pepmod_id, peptide_id)) %>%
  dplyr::select(-peptide_id)

msdata_full$protregroups <- msdata_full$protregroups %>%
  dplyr::mutate(gene_label = strlist_label2(gene_names),
                protein_label = strlist_label2(protein_names),
                protein_description = strlist_label2(protein_descriptions),
                is_viral = replace_na(str_detect(organism, "MPXV"), FALSE),
                protac_label = strlist_label2(majority_protein_acs),
                protregroup_label = case_when(is_viral ~ protein_label,
                  !is.na(gene_label) ~ gene_label,
                  !is.na(protac_label) ~ protac_label,
                  TRUE ~ str_c('#', protregroup_id))) %>%
  dplyr::left_join(dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmods) %>%
                     dplyr::group_by(protregroup_id) %>%
                     dplyr::summarise(npeptides = n_distinct(peptide_id),
                                      npepmods = n_distinct(pepmod_id),
                                      npeptides_unique = n_distinct(peptide_id[is_specific]),
                                      npepmods_unique = n_distinct(pepmod_id[is_specific])) %>%
                     dplyr::ungroup() %>%
                     dplyr::mutate(npeptides_unique_razor = npeptides_unique,
                                   npeptides_razor = 0L,
                                   npepmods_unique_razor = npepmods_unique,
                                   npepmods_razor = 0L))

msdata_full$pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::mutate(is_idented = str_detect(ident_type, "MSMS"))

# prepare protgroup intensities (wider format: all mstags in one row)
intensity_prespec_df <- tibble(.name = msdata_colgroups$LFQ) %>%
  extract(.name, c("mstag", "msexperiment_mq"), remove=FALSE,
          str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)")) %>%
  mutate(.value = str_c("intensity.", mstag)) %>%
  dplyr::inner_join(select(msruns.df, msexperiment_mq, rawfile))

# prepare protgroup intensities (longer format: each mstag on its own row)
protgroup_intensities_all.df <- tidyr::pivot_longer_spec(
  dplyr::select(msdata.wide, protgroup_id, !!msdata_colgroups$LFQ),
  mutate(intensity_prespec_df, .value = "intensity")) %>%
  select(-mstag)

msdata_full$protgroup_intensities <- protgroup_intensities_all.df %>%
  dplyr::semi_join(select(msdata_full$msruns, msexperiment_mq, rawfile)) %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::select(-rawfile)

msdata_full$protregroup_idents <- dplyr::inner_join(msdata_full$protregroup2pepmod,
                                                    msdata_full$pepmodstate_intensities) %>%
  dplyr::group_by(msexperiment, protregroup_id) %>%
  dplyr::summarise(npepmods_quanted = sum(!is.na(intensity)),
                   nspecpepmods_quanted = sum(!is.na(intensity) & is_specific),
                   npepmods_idented = sum(is_idented),
                   nspecpepmods_idented = sum(is_idented & is_specific)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ident_type = factor(if_else(nspecpepmods_idented > 0, "By MS/MS", "By matching")))

# condition = treatment X timepoint
conditions.df <- dplyr::select(msdata_full$msruns, condition, treatment, timepoint, timepoint_num) %>%
  dplyr::distinct() %>% 
  arrange(treatment, timepoint_num)
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
  left_join(dplyr::mutate(dplyr::select(conditions.df, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint_num)),
            dplyr::mutate(dplyr::select(conditions.df, metacondition_rhs = condition, treatment_rhs = treatment, timepoint_lhs = timepoint_num))) %>%
    filter((as.integer(treatment_lhs) > as.integer(treatment_rhs)) &
             !((treatment_lhs == "infected") & (treatment_rhs != "mock"))) %>%
    mutate(timepoint_rhs = timepoint_lhs,
           contrast_kind = "treatment_vs_treatment",
           contrast = str_c(treatment_lhs, "_vs_", treatment_rhs, "@", timepoint_rhs, "h")),
  # all timepoints of the same treatment
  left_join(dplyr::mutate(dplyr::select(conditions.df, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint_num)),
            dplyr::mutate(dplyr::select(conditions.df, metacondition_rhs = condition, treatment_lhs = treatment, timepoint_rhs = timepoint_num))) %>%
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

msdata <- import_msglm_data(msdata_full, msglm_def, object="protregroup", quantobject="pepmodstate", verbose = TRUE) 

##########
## normalization
msdata4norm.df <- msdata_full$pepmodstate_intensities %>% ungroup() %>%
  dplyr::semi_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id)) %>%
  dplyr::semi_join(dplyr::semi_join(msdata_full$pepmods,
                                    dplyr::filter(msdata_full$protregroups, !is_reverse & !is_contaminant & !is_viral) %>%
                                      dplyr::inner_join(msdata_full$protregroup2pepmod) %>% dplyr::select(pepmod_id)) %>%
                     dplyr::filter(!is_reverse & !is_contaminant) %>%
                     dplyr::select(pepmod_id))


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
