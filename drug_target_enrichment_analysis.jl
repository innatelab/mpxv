# Enrichment analysis for the drug targets from network diffusion (Valter)
#
# Author: yhuang
###############################################################################

#---
proj_info = (id = "mpxv",
             modelobj = "protein",
             analysis_ver = "20230306b")

const base_scripts_path = "/home/ge54heq/projects"
const base_analysis_path = "/pool/analysis/yhuang"

using Pkg
Pkg.activate(@__DIR__)

using Revise
using Distances, DataFrames, CategoricalArrays, RData, JSON, CodecZstd, StatsBase
using JLD2
using DelimitedFiles

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const reports_path = joinpath(analysis_path, "reports")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")
const party3rd_data_path = "/pool/pub3rdparty"

@info "Project '$(proj_info.id)' (analysis ver=$(proj_info.analysis_ver)))"

includet(joinpath(misc_scripts_path, "frame_utils.jl"));
includet(joinpath(misc_scripts_path, "msglm_utils.jl"));

objid_col = Symbol("protein_id")

#---
import_rdata = true
if import_rdata

@info "Loading $(joinpath(reports_path, "downstream_analysis", "drug_targets_from_Valter_long_$(proj_info.analysis_ver)")).txt..."
data, header = readdlm(joinpath(reports_path, "downstream_analysis", "drug_targets_from_Valter_long_$(proj_info.analysis_ver).txt"), header = true)
contrasts_df = DataFrame(data, vec(header))
contrasts_df.protein_id = Int64.(contrasts_df.protein_id)
contrasts_df.protein_ac = string.(contrasts_df.protein_ac)
contrasts_df.gene_name = string.(contrasts_df.gene_name)
contrasts_df.comparison = string.(contrasts_df.comparison)
contrasts_df.is_hit = contrasts_df.is_hit .== "TRUE"
obj2protac_df = select(contrasts_df, :protein_id, :gene_name, :protein_ac) |> unique!


Revise.includet(joinpath(misc_scripts_path, "gmt_reader.jl"));
#add https://github.com/alyst/OptEnrichedSetCover.jl
include(joinpath(misc_scripts_path, "optcover_utils.jl"));
include(joinpath(misc_scripts_path, "omics_collections.jl"));

@info "Loading Human annotations..."
# human mappings from http://download.baderlab.org/EM_Genesets/December_01_2021/Human/UniProt/
# FIXME using all evidence codes
genesets_df, genesets_coll = GMT.read(String,
        joinpath(party3rd_data_path, "Human_GO_AllPathways_with_GO_iea_December_01_2021_UniProt.gmt"),
        id_col = :term_id, src_col = :term_src);
# Note: there's a spillover of GOBP terms into GOMF, might affect the GOMF outcome. The Balder lab is gradually correcting it, so make sure to always use the latest genesets from them!
#=
genesets_df, genesets_coll = GMT.read(String,
        joinpath(base_analysis_path, "pub3rdparty", "Human_GO_AllPathways_with_GO_iea_February_08_2023_UniProt.gmt"),
        id_col = :term_id, src_col = :term_src);
duplicates_df = genesets_df[findall(nonunique(genesets_df, :name)), :]
duplicates_df2 = genesets_df[findall(nonunique(genesets_df, :term_id)), :]
=#

# strip Reactome version
genesets_df.term_id = [ifelse(r.term_src == "Reactome", replace(r.term_id, r"\.\d+$" => ""), r.term_id)
                       for r in eachrow(genesets_df)]
genesets_coll = Dict(ifelse(contains(k, r"^R-HSA-"), replace(k, r"\.\d+$" => ""), k) => v
                     for (k, v) in pairs(genesets_coll))

pcomplexes_df, pcomplex_iactors_df, pcomplex_iactor2ac_df =
    OmicsCollections.ppicollection(joinpath(party3rd_data_path, "complexes_20191217.RData"), seqdb=:uniprot);
pcomplexes_df[!, :coll_id] .= "protein_complexes";
# make complexes collections, keep complexes with at least 2 participants
uprot_pcomplex_coll = FrameUtils.frame2collection(innerjoin(pcomplex_iactors_df, pcomplex_iactor2ac_df,
            on=[:file, :entry_index, :interaction_id, :interactor_id]),
            set_col=:complex_id, obj_col=:protein_ac, min_size=2)

protac_sets = merge!(genesets_coll, uprot_pcomplex_coll)

terms_df = vcat(rename(genesets_df[!, [:term_src, :term_id, :name, :descr]],
                       :term_src => :coll_id, :name=>:term_name, :descr=>:term_descr),
                rename(pcomplexes_df[!, [:coll_id, :complex_id, :interaction_label, :interaction_name]],
                       :complex_id=>:term_id, :interaction_label=>:term_name, :interaction_name=>:term_descr));
protac2term_df = FrameUtils.collection2frame(protac_sets, terms_df,
                                             setid_col=:term_id, objid_col=:protein_ac)

# link protein group IDs to annots and create protgroup collections
obj2term_df = select!(innerjoin(obj2protac_df, protac2term_df, on = :protein_ac),
                      Not([:protein_ac])) |> unique!
protac_colls = FrameUtils.frame2collections(protac2term_df, obj_col=:protein_ac,
                                            set_col=:term_id, coll_col=:coll_id)
obj_colls = FrameUtils.frame2collections(obj2term_df, obj_col=objid_col,
                                         set_col=:term_id, coll_col=:coll_id)

@info "Preparing mosaics..."
observed_protacs = Set(string.(obj2protac_df.protein_ac)) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls, protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);

jlannots_path = joinpath(scratch_path, "$(proj_info.id)_$(proj_info.analysis_ver)_annotations.jld2")
@info "Serializing annotations to JLD2 format: $jlannots_path"
@save(jlannots_path,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      obj2protac_df, contrasts_df)

else

jlannots_path = joinpath(scratch_path, "$(proj_info.id)_$(proj_info.analysis_ver)_annotations.jld2")
@info "Deserializing annotations from JLD2 format: $jlannots_path"
@load(jlannots_path,
proj_info, protac_colls, obj_colls, obj_mosaics,
obj2term_df, terms_df,
obj2protac_df, contrasts_df)

end

ENV["MKL_NUM_THREADS"] = 1
using OptEnrichedSetCover
include(joinpath(misc_scripts_path, "omics_collections.jl"));
cover_params = CoverParams(setXset_factor=0.5,
                           uncovered_factor=0.0, covered_factor=0.0)
# GO enrichment for the drug targets
obj_drug_sets = FrameUtils.frame2collection(
    filter(r -> r.is_hit, contrasts_df),
 set_col=:comparison, obj_col=:protein_id)  

obj_drug_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by contract hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_drug_sets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

obj_drug_hit_mosaics_v = collect(pairs(obj_drug_hit_mosaics))
obj_drug_hit_covers_v = similar(obj_drug_hit_mosaics_v, Pair)


Threads.@threads for i in eachindex(obj_drug_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_drug_hit_mosaics_v[i]
    @info "Covering $mosaic_name by drug targets..."
    obj_drug_hit_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.1, 0.1], MaxSteps=2_000_000, WeightDigits=2, NWorkers=1),
            true)
end
obj_drug_hit_covers = Dict(k => v for (k, v) in obj_drug_hit_covers_v)

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_drug_hits_$(proj_info.analysis_ver)_covers.jld2")
@save(hit_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      obj2protac_df, contrasts_df, 
      obj_drug_sets, 
      obj_drug_hit_mosaics, 
      cover_params, obj_drug_hit_covers)
if !@isdefined(obj_drug_hit_covers)
using JLD2, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      obj2protac_df, contrasts_df, 
      obj_drug_sets, 
      obj_drug_hit_mosaics, 
      cover_params, obj_drug_hit_covers)
end

include(joinpath(misc_scripts_path, "optcover_utils.jl"));
#FIXME: the term_collection and coll_id should also match! Currently they mismatch when the term_id is the same, but in reality only term_collection correctly indicates the source of this term.
@info "Preparing protgroup->gene_name map..."
obj_id2name = Dict(r.protein_id => r.gene_name
                   for r in eachrow(obj2protac_df)) 

obj_drug_hit_covers_df = OptCoverUtils.covers_report(
    obj_drug_hit_covers, obj_drug_sets, 
    obj_colls, 
    obj_id2name, terms_df,
    cover_params = cover_params,
    experimentid_col=[:comparison], weightedset_col_prefix="hit");
obj_drug_hit_covers_df.intersect_genes = [join(unique(vcat(split.(split(genes, ' '), ';')...)), ' ') for genes in obj_drug_hit_covers_df.intersect_genes]

obj_drug_hit_covers_signif_df = combine(groupby(obj_drug_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    df = OptCoverUtils.filter_multicover(coll_df, set_cols=[:comparison],
                                         max_term_pvalue=1E-3, max_set_pvalue=1E-2, min_set_overlap=nothing)
    return select!(df, Not(:term_collection))
end

using CSV

CSV.write(joinpath(reports_path, "$(proj_info.id)_drug_hit_enrichment_$(proj_info.analysis_ver).txt"),
obj_drug_hit_covers_df[obj_drug_hit_covers_df.nmasked .> 0, :],
missingstring="", delim='\t');

CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_drug_hit_enrichment_signif_$(proj_info.analysis_ver).txt"),
          obj_drug_hit_covers_signif_df[obj_drug_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');
