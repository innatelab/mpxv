proj_info = (id = "mpxv",
             oesc_ver = "20230804")

datasets = Dict(
    :mpxv => (fit_ver="20221104", data_ver = "20221104",
           analysis=:msglm, ptm=nothing,
           datatype=:fp, entity=:protein,
           label="MPXV",
         color="#000000"),
    :vacv => (fit_ver=nothing, 
                 analysis=:literature, ptm=nothing,
                 datatype=:fp, entity=:protein,
                 label="VACV",
                 color="#000000"),
    :mva => (fit_ver=nothing,
                analysis=:literature, ptm=nothing,
                datatype=:fp, entity=:protein,
                label="MVA",
                color="#000000"),
)

using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, CSV, DataFrames, FastaIO
using JLD2
using StatsBase

@info "Project '$(proj_info.id)' analysis version=$(proj_info.oesc_ver)"

const base_scripts_path = "/home/ge54heq/projects"
const base_analysis_path = "/pool/analysis/yhuang"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")
const party3rd_data_path = "/pool/pub3rdparty"

includet(joinpath(misc_scripts_path, "frame_utils.jl"));
includet(joinpath(misc_scripts_path, "msglm_utils.jl"));
includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));
includet(joinpath(misc_scripts_path, "fasta_reader.jl"));
include(joinpath(misc_scripts_path, "ms_import.jl"));
includet(joinpath(misc_scripts_path, "protgroup_assembly.jl"));
includet(joinpath(misc_scripts_path, "protgroup_crossmatch.jl"));

# load data and fits
msglm_rdata = Dict(begin
    @info "Loading $dsname analysis (fit_ver=$(dsinfo.fit_ver))..."
    dsname => (input_data = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(dsinfo.datatype)_$(dsinfo.fit_ver).RData"), convert=true),
               full_data = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(dsinfo.datatype)_$(dsinfo.data_ver).RData"), convert=true),
               fit = load(joinpath(results_path, "$(proj_info.id)_msglm_fit_$(dsinfo.datatype)_$(dsinfo.fit_ver).RData"), convert=true))
end for (dsname, dsinfo) in pairs(datasets) if dsinfo.analysis == :msglm);

# fix some data differences
for (dsname, rdata) in pairs(msglm_rdata)
    @info "Fixing $dsname"
    msdata_full = rdata.full_data["msdata_full"]

    proteins_df = msdata_full["proteins"]
    proteins_df.organism = [ismissing(org) ? org : replace(org, r"\s+OX=\d+$" => "") for org in proteins_df.organism]
    
    hasproperty(proteins_df, :genename) && rename!(proteins_df, :genename => :gene_name)
    if hasproperty(proteins_df, :protein_name) && hasproperty(proteins_df, :protein_description)
        rename!(proteins_df, :protein_name => :protein_code, :protein_description => :protein_name)
    end

    peptides_df = msdata_full["peptides"]
    hasproperty(peptides_df, :peptide_seq) || rename!(peptides_df, :seq => :peptide_seq)

    pepmodstates_df = msdata_full["pepmodstates"]
    if haskey(msdata_full, "pepmods")
        pepmods_df = msdata_full["pepmods"]
        missing_pms_cols = intersect(setdiff([:peptide_id, :nselptms], propertynames(pepmodstates_df)), propertynames(pepmods_df))
        if !isempty(missing_pms_cols)
            pepmodstates_df = innerjoin(pepmodstates_df, pepmods_df[!, [[:pepmod_id]; missing_pms_cols]], on=:pepmod_id)|> unique! #sometimes contaminants might be counted twice in pepmods
            @assert nrow(pepmodstates_df) == nrow(msdata_full["pepmodstates"])
            msdata_full["pepmodstates"] = pepmodstates_df
        end
    end

    pms_intensities_df = msdata_full["pepmodstate_intensities"]
    hasproperty(pms_intensities_df, :psm_pvalue) || insertcols!(pms_intensities_df, :psm_pvalue =>1) #fp doesn't have psm_pvalue, use 1 instead. Later the peptides from fp won't be influcing the protgroup of ptms anyway
    if !hasproperty(pms_intensities_df, :peptide_id)
        pms_intensities_df = innerjoin(pms_intensities_df, pepmodstates_df[!, [:pepmodstate_id, :peptide_id]], on=:pepmodstate_id)
        @assert nrow(pms_intensities_df) == nrow(msdata_full["pepmodstate_intensities"])
        msdata_full["pepmodstate_intensities"] = pms_intensities_df
    end
end

peptide2protein_df = combine(groupby(reduce(vcat, begin
    @info "Processing peptides of $dsname"
    isptm = !isnothing(datasets[dsname].ptm)
    msdata = dsdata.full_data["msdata_full"]
    if isptm
    pep2prot_df = copy(msdata["peptide2protein"], copycols=false)
    else
        pep2prot_df = select(msdata["peptides"], [:peptide_id, :peptide_seq, :protein_acs])
        pep2prot_df = filter!(r -> !ismissing.(r.protein_acs), pep2prot_df)
        pep2prot_df[!, "protein_ac"] = split.(pep2prot_df.protein_acs, ";")
        pep2prot_df = flatten(pep2prot_df, "protein_ac")
    end
    if !hasproperty(pep2prot_df, :peptide_seq)
        pep2prot_df = innerjoin(msdata["peptide2protein"], msdata["peptides"],
                                on=:peptide_id)
    end
    if isptm # only use modified peptides for PTM datasets
        pepmodstates_df = filter(r -> r.nselptms > 0, msdata["pepmodstates"])
        # remove phosphosites from ubi data
        (dsname == :ubi) && filter!(r -> occursin("[GlyGly", r.pepmod_seq), pepmodstates_df)
        select!(pepmodstates_df, [:pepmodstate_id, :peptide_id])
    else
        pepmodstates_df = select(msdata["pepmodstates"], [:pepmodstate_id, :peptide_id])
    end
    pms_intensities_df = copy(msdata["pepmodstate_intensities"], copycols=false)
    if !hasproperty(pms_intensities_df, :peptide_seq)
        pms_intensities_df = innerjoin(semijoin(pms_intensities_df, pepmodstates_df, on=:pepmodstate_id),
                                       select(msdata["peptides"], [:peptide_id, :peptide_seq]), on=:peptide_id)
    end
    pepstats_df = combine(groupby(pms_intensities_df, :peptide_seq),
                          :psm_pvalue => (x -> minimum(coalesce(x,1))) => :psm_pvalue) #note that some peptides has nothing in psm_pvalue at all
    # set the rank of non-PTM peptides to -1, so that proteome-only observed peptides are ignored (was: still used to define
    # protein groups, but they would not split the protein group derived from PTM datasets)
    pepstats_df.peptide_rank = ifelse.(coalesce.(pepstats_df.psm_pvalue, 1.0) .<= 1E-3, ifelse(isptm, 1, 2), -1)
    pep2prot_df = innerjoin(pep2prot_df, pepstats_df, on=:peptide_seq)
    pep2prot_df.peptide_seq = replace.(pep2prot_df.peptide_seq, Ref('_' => ""))
    select!(pep2prot_df, [:protein_ac, :peptide_seq, :peptide_rank])
    unique!(pep2prot_df)
    pep2prot_df[!, :dataset] .= dsname
    @info "  $(nrow(pep2prot_df)) associations of $(length(unique(pep2prot_df.peptide_seq))) peptide(s) to $(length(unique(pep2prot_df.protein_ac))) protein(s)"
    pep2prot_df
end for (dsname, dsdata) in pairs(msglm_rdata)), [:protein_ac, :peptide_seq]),
    :peptide_rank => (x -> any(>(0), x) ? minimum(filter(>(0), x)) : -1) => :peptide_rank,
    :dataset => (x -> join(sort(x), ' ')) => :datasets)

peptide2proteins = Dict(df.peptide_seq[1] => (Set{String}(df.protein_ac), df.peptide_rank[1])
                        for df in groupby(peptide2protein_df, :peptide_seq))

proteins_df = reduce(vcat, begin
    select(rdata.full_data["msdata_full"]["proteins"], [:protein_ac, :gene_name, :protein_code, :protein_name, :protein_existence, :src_db, :is_contaminant, :is_viral, :organism], copycols=false)
end for (dsname, rdata) in msglm_rdata) |> unique!
proteins_df.protein_ac_noiso = Fasta.strip_uniprot_isoform.(proteins_df.protein_ac)
for noiso_group_df in groupby(proteins_df, :protein_ac_noiso) # extrapolate best PE across isoforms
    pe = noiso_group_df.protein_existence
    if any(ismissing, pe) && any(!ismissing, pe)
        noiso_group_df.protein_existence .= coalesce.(pe, minimum(skipmissing(pe)))
    elseif any(ismissing, pe) && any(noiso_group_df.is_viral)
        noiso_group_df.protein_existence .= 1
    end
end
# calculate AC ranks - the more canonical proptein is, the better
protein_ac_ranks = Dict(r.protein_ac => ProtgroupXMatch.rank_uniprot(r) for r in eachrow(proteins_df))

mpxv_protgroups_df = rename!(copy(msglm_rdata[:mpxv].full_data["msdata_full"]["protregroups"]),
                           :protregroup_id => :protgroup_id)

literature_rdata = load(joinpath(results_path, "proteome_literature_20230804.RData"), convert=true)
vacv_contrasts_df = copy(literature_rdata["soday4enrichment.df"], copycols=false)
mva_contrasts_df = copy(literature_rdata["albarnaz4enrichment.df"], copycols=false)
                            
vacv_protgroups_df = unique(select(vacv_contrasts_df, :majority_protein_acs, :organism, :gene_names))
mva_protgroups_df = unique(select(mva_contrasts_df, :majority_protein_acs, :organism, :gene_names))

vacv_protgroups_df[!, :is_contaminant] .= false
vacv_protgroups_df[!, :is_reverse] .= false
vacv_protgroups_df.protein_acs = vacv_protgroups_df.majority_protein_acs
vacv_protgroups_df.protgroup_id = 1:nrow(vacv_protgroups_df)
mva_protgroups_df[!, :is_contaminant] .= false
mva_protgroups_df[!, :is_reverse] .= false
mva_protgroups_df.protein_acs = mva_protgroups_df.majority_protein_acs
mva_protgroups_df.protgroup_id = 1:nrow(mva_protgroups_df)

protgroups_dfs = Dict(
    :mva => mva_protgroups_df,
    :vacv => vacv_protgroups_df,
    :mpxv => mpxv_protgroups_df,
)

pg_matches_df = ProtgroupXMatch.match_protgroups(collect(pairs(protgroups_dfs)), protein_ac_ranks);
objid_col = :protgroup_id_united;
pg_matches_df.protgroup_id_united = 1:nrow(pg_matches_df)

# add protgroup ids of each dataset and add protgroup_id_common to each dataset
for (dsname, pg_df) in pairs(protgroups_dfs)
    rowix_col = Symbol("rowix_", dsname)
    if !hasproperty(pg_matches_df, rowix_col)
        @warn "Dataset $dsname key does not exist, skipping"
        continue
    end
    pg_col = Symbol("protgroup_id_", dsname)
    pg_matches_df[!, pg_col] = missings(nonmissingtype(eltype(pg_df.protgroup_id)), nrow(pg_matches_df))
    pg_df[!, :protgroup_id_united] = missings(Int, nrow(pg_df))
    for (i, rowix) in enumerate(pg_matches_df[!, rowix_col])
        if !ismissing(rowix)
            pg_matches_df[i, pg_col] = pg_df[rowix, :protgroup_id]
            pg_df[rowix, :protgroup_id_united] = pg_matches_df[i, :protgroup_id_united]
        end
    end
    select!(pg_matches_df, Not(rowix_col))
end
pg_matches_long_df = FrameUtils.pivot_longer(pg_matches_df, [:protgroup_id_united, :protein_acs, :majority_protein_acs, :gene_names,
                                             :is_contaminant, :is_reverse],
                        measure_vars_regex=r"^(?<value>rowix|protgroup_id|pgrank|acrank)_(?<var>[^u].+)$",
                        var_col=:dataset)
pg_matches_long_df.dataset = Symbol.(pg_matches_long_df.dataset)
pg_matches_long_expanded_df = DelimDataUtils.expand_delim_column(
    pg_matches_long_df, key_col=[:dataset, :protgroup_id_united, :protgroup_id],
    list_col=:protein_acs, elem_col=:protein_ac)
#=
length(union(
    unique(filter(r -> occursin(r"SARS_.+_vs_mock", r.contrast) && r.std_type == "median" && r.change in ["+", "-"],
              innerjoin(msglm_rdata[:cov2ts_proteome].fit["object_contrasts.df"],
                        select(pg_matches_df, [:protgroup_id_united, :protgroup_id_cov2ts_proteome]),
                        on=[:protgroup_id => :protgroup_id_cov2ts_proteome])).protgroup_id_united),
    unique(filter(r -> occursin(r"SARS_.+_vs_mock", r.contrast) && r.std_type == "median" && r.change in ["+", "-"],
              innerjoin(msglm_rdata[:cov2el_proteome].fit["object_contrasts.df"],
                        select(pg_matches_df, [:protgroup_id_united, :protgroup_id_cov2el_proteome]),
                        on=[:protregroup_id => :protgroup_id_cov2el_proteome])).protgroup_id_united)
    ))
=#

obj2protac_df = unique!(dropmissing!(select(pg_matches_long_expanded_df,
                                            [objid_col, :protein_ac])))

nonunique_matches_dfs = [dsname => begin
    pg_col = Symbol("protgroup_id_", dsname)
    combine(groupby(filter(r -> !ismissing(r[pg_col]), pg_matches_df), pg_col)) do df
        return nrow(df) > 1 ? df : df[1:0, :]
    end
end for (dsname, _) in pairs(protgroups_dfs)]

contrast_cols = [:contrast, :timepoint_lhs, :timepoint_rhs, :treatment_lhs, :treatment_rhs]

obj_contrasts_msglm_dfs = Dict(ds => begin
    @info "Processing $ds results..."
    orig_objcontrasts_df = rdata.fit["object_contrasts.df"]
    sel_cols = [:ci_target; contrast_cols;
        [:object_id,
        :median, :p_value,
        :is_hit, :change, :is_hit_nomschecks, :is_signif]]
    for col in [:ptm_type, :ptmngroup_id, :ptmn_id, :protregroup_id]
        hasproperty(orig_objcontrasts_df, col) && push!(sel_cols, col)
    end
    objcontrasts_df = select(orig_objcontrasts_df, sel_cols)
    objcontrasts_df[!, :dataset] .= ds
    objcontrasts_df.change = ifelse.(objcontrasts_df.is_hit,
            ifelse.(objcontrasts_df.median .> 0, "+", "-"), ".")
    @show nrow(objcontrasts_df)
    isptm = !isnothing(datasets[ds].ptm)
    if isptm
        @info "  Associating PTM groups to protein groups..."
        objcontrasts_df = innerjoin(objcontrasts_df,
                                    ptmngroup2protgroup_dfs[ds], on=:ptmngroup_id)
        #@assert hasproperty(objcontrasts_df, :ptm_type)
        objcontrasts_df[!, :ptm_type] .= datasets[ds].ptm
    else
        if hasproperty(objcontrasts_df, :protregroup_id) &&
          !hasproperty(objcontrasts_df, :protgroup_id)
            rename!(objcontrasts_df, :protregroup_id => :protgroup_id)
        end
        objcontrasts_df.ptmgroup_id = missings(Int, nrow(objcontrasts_df))
        objcontrasts_df.ptmn_id = missings(Int, nrow(objcontrasts_df))
        objcontrasts_df.ptm_type = missings(String, nrow(objcontrasts_df))
    end
    @show nrow(objcontrasts_df)
    @info "  Associating dataset-specific protein groups to united protein groups..."
    objcontrasts_df = leftjoin(objcontrasts_df,
                               select(protgroups_dfs[isptm ? :ptm : ds], [:protgroup_id, :protgroup_id_united], copycols=false),
                               on=:protgroup_id)
    @show nrow(objcontrasts_df)
    objcontrasts_df
end for (ds, rdata) in pairs(msglm_rdata))

sel_ci_target = "average"
# in FP one protgroup could belong to a single +/-/. group
obj_contrasts_mpxv_agg_df = combine(groupby(filter(r -> (r.ci_target == sel_ci_target) && !ismissing(r.protgroup_id_united), obj_contrasts_msglm_dfs[:mpxv]),
   [[:ci_target, :dataset, :ptm_type]; contrast_cols; :protgroup_id_united])
) do df
    min_ix = findmin(df.p_value)[2]
    return df[min_ix:min_ix, :]
end

obj_contrasts_vacv_agg_df = combine(groupby(innerjoin(vacv_contrasts_df, protgroups_dfs[:vacv],
                                             on=[:majority_protein_acs, :gene_names, :organism]),
                        [contrast_cols; :protgroup_id_united])) do df
   max_ix = findmax(abs.(df.log2FC))[2]
   return df[max_ix:max_ix, :]
end
obj_contrasts_vacv_agg_df[!, :dataset] .= :vacv
obj_contrasts_vacv_agg_df.timepoint_lhs = string.(Int.(obj_contrasts_vacv_agg_df[!,:timepoint_lhs]))
obj_contrasts_vacv_agg_df.timepoint_rhs = string.(Int.(obj_contrasts_vacv_agg_df[!,:timepoint_rhs]))

obj_contrasts_mva_agg_df = combine(groupby(innerjoin(mva_contrasts_df, protgroups_dfs[:mva],
                                             on=[:majority_protein_acs, :gene_names, :organism]),
                        [contrast_cols; :protgroup_id_united])) do df
   max_ix = findmax(df.change)[2] #this study doesn't offer fold change, but there's only one protein that has both up and down so I just picked "up"...
   return df[max_ix:max_ix, :]
end
obj_contrasts_mva_agg_df[!, :dataset] .= :mva
obj_contrasts_mva_agg_df.timepoint_lhs = string.(Int.(obj_contrasts_mva_agg_df[!,:timepoint_lhs]))
obj_contrasts_mva_agg_df.timepoint_rhs = string.(Int.(obj_contrasts_mva_agg_df[!,:timepoint_rhs]))

united_cols = [:dataset; contrast_cols; :protgroup_id_united; :is_hit; :change]
obj_contrasts_united_df = vcat(select(obj_contrasts_mpxv_agg_df, united_cols, copycols=false),
                               select(obj_contrasts_vacv_agg_df, united_cols, copycols=false),
                               select(obj_contrasts_mva_agg_df, united_cols, copycols=false))
filter!(r -> r.treatment_rhs == "mock" && r.timepoint_lhs == r.timepoint_rhs,
        obj_contrasts_united_df)
countmap(collect(zip(obj_contrasts_united_df.is_hit, obj_contrasts_united_df.change)))

contrasts_df = unique!(select(obj_contrasts_united_df, [:dataset; contrast_cols; :change]))
contrasts_df.dataset = String.(contrasts_df.dataset)
contrasts_df.timepoint_lhs = Vector{Union{String}}(contrasts_df[!,:timepoint_lhs])
contrasts_df.timepoint_lhs = parse.(Int, contrasts_df.timepoint_lhs)
contrasts_df.timepoint_rhs = Vector{Union{String}}(contrasts_df[!,:timepoint_rhs])
contrasts_df.timepoint_rhs = parse.(Int, contrasts_df.timepoint_rhs)
contrasts_df.change_alt = getindex.(Ref(Dict("+" => "▲",  "-" => "▼", "." => ".")),
                                        contrasts_df.change)

include(joinpath(misc_scripts_path, "optcover_utils.jl"));
includet(joinpath(misc_scripts_path, "gmt_reader.jl"));
includet(joinpath(misc_scripts_path, "omics_collections.jl"));

@info "Loading Human annotations..."
# human mappings from http://download.baderlab.org/EM_Genesets/December_01_2021/Human/UniProt/
# FIXME using all evidence codes
genesets_df, genesets_coll = GMT.read(String,
        joinpath(party3rd_data_path, "Human_GO_AllPathways_with_GO_iea_December_01_2021_UniProt.gmt"),
        id_col = :term_id, src_col = :term_src);
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

@info "Preparing hit sets"
ObjectType = eltype(obj2protac_df[!, objid_col])
obj_hit_sets = Dict{Tuple{String, String, String}, Set{ObjectType}}()
for hits_df in groupby(filter(r -> coalesce(r.is_hit, false), obj_contrasts_united_df),
                       [:dataset, :contrast, :change])
    obj_hit_sets[(string(hits_df[1, :dataset]),
                  hits_df[1, :contrast], hits_df[1, :change])] =
        Set(skipmissing(hits_df[!, objid_col]))
end

# only relevant ones
obj_hit_selsets = obj_hit_sets

@info "Preparing mosaics..."
observed_protacs = Set(obj2protac_df.protein_ac) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls,
                                      protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);

# remove broad terms larger than 200 elements (too unspecific)
obj_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_hit_selsets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200, verbose=true)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.5,
                           uncovered_factor=0.0, covered_factor=0.0)#, covered_factor=0.002)

ENV["MKL_NUM_THREADS"] = 1

obj_hit_mosaics_v = collect(pairs(obj_hit_mosaics))
obj_hit_covers_v = similar(obj_hit_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_hit_mosaics_v[i]
    @info "Covering $mosaic_name by hits..."
    obj_hit_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2,
                                    NWorkers=1,#Threads.nthreads()-1,
                                    MaxRestarts=200),
            true)
end
obj_hit_covers = Dict(k => v for (k, v) in obj_hit_covers_v)

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_proteome_literature_hit_covers_$(proj_info.oesc_ver).jld2")
@save(hit_covers_filename,
      proj_info, datasets, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      #objects_df,
      protgroups_dfs, pg_matches_df,
      contrasts_df,
      obj_contrasts_united_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
if !@isdefined(obj_hit_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
      proj_info, datasets, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      #objects_df,
      protgroups_dfs, pg_matches_df,
      contrasts_df,
      obj_contrasts_united_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
end

include(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing protgroup?gene_name map..."
obj_id2name = Dict(r.protgroup_id_united => !ismissing(r.gene_names) ? DelimDataUtils.rejoin_unique_substrings([r.gene_names]) :
                    !ismissing(r.majority_protein_acs) ? DelimDataUtils.rejoin_unique_substrings([r.majority_protein_acs]) :
                    string("PGU_", r.protgroup_id_united)
                   for r in eachrow(filter(r -> !ismissing(r.protgroup_id_united), pg_matches_df)))

obj_hit_covers_df = innerjoin(
    OptCoverUtils.covers_report(
    obj_hit_covers, obj_hit_selsets, obj_colls, #obj_mosaics, 
    obj_id2name, terms_df,
    cover_params = cover_params,
    experimentid_col=[:dataset, :contrast, :change],
    weightedset_col_prefix="hit"),
    contrasts_df, on=[:dataset, :contrast, :change])
obj_hit_covers_df.intersect_genes = [join(unique(vcat(split.(split(genes, ' '), ';')...)), ' ') for genes in obj_hit_covers_df.intersect_genes]

# don't remove the sets since they are timecourses timepoints
obj_hit_covers_signif_df = combine(groupby(obj_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:dataset, :contrast, :change],
                        max_term_pvalue=1E-3, max_set_pvalue=nothing, min_set_overlap=nothing),
                    Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_proteome_literature_hit_covers_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_df[obj_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_proteome_literature_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_signif_df[obj_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_plots.jl"))
include(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

using PlotlyJS, TextWrap

heatmap_layout_attrs = Dict(
    ("SigDB_C2", true) => Dict(:margin_l => 500),
    ("SigDB_C2", false) => Dict(:margin_l => 500),
    #("Reactome", true) => Dict(:margin_l => 600),
    #("Reactome", false) => Dict(:margin_l => 600),
    #("GO_CC", true) => Dict(:margin_l => 200),
    #("GO_CC", false) => Dict(:margin_l => 200),
)

stylize_dataset(ds) =
    "<span style=\"font-weight: bold; color: $(datasets[Symbol(ds)].color);\">" * datasets[Symbol(ds)].label * "</span>"

stylize_contrast(str) = foldl(replace, [
    r"(MPXVAC?)_vs_mock@(\d+)h" => s"\1:<span style=\"font-weight: bold; color: black;\">\2</span>h",
    "infection:" => "<span style=\"font-wieght: bold; color: #F9CB40;\">infection</span> ",
    ],
    init = str)

function process_contrast_axis(contrast_df)
    contrast_df,
    stylize_dataset.(contrast_df.dataset) .* ":&nbsp;" .*
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change),
    stylize_dataset.(contrast_df.dataset) .* ":&nbsp;" .*
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change)#stylize_effect.(effect_df.effect)
end

dataset_order = Dict(:mpxv => 1, :vacv => 2, :mva => 3)

heatmaps_path = joinpath(plots_path, "proteome_literature_hits_oesc_$(sel_ci_target)_$(proj_info.oesc_ver)")
isdir(heatmaps_path) || mkdir(heatmaps_path)

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")hit heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? obj_hit_covers_signif_df : obj_hit_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end

    for outformat in ("html", "pdf", "svg")
        coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="protein",
            experiment_axis_title = "contrast",
            experiment_cols = [:dataset, :contrast, :treatment_lhs, :timepoint_lhs, :change, :nhit],
            process_experiment_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 200), #bottom margin, increase if the bottom is cut off
            transpose=false,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",#outformat in ["svg", "pdf"] ? "#000" : "#BBB",
            cell_width=40, cell_height=30, gridwidth=2,
            experiment_order=contrasts -> begin
                contrasts.dataset_order = [dataset_order[Symbol(r.dataset)]
                                            for r in eachrow(contrasts)]
                return sortperm(contrasts, [:change, :dataset_order, :timepoint_lhs, :treatment_lhs, :dataset])
            end)
        
        (coll_heatmap === nothing) && continue
        for (k, v) in [#:width=>800, :height=>400,
                       :margin_r=>80, #:margin_t=>20,
                       :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
            coll_heatmap.plot.layout[k] = v
        end
        
    for outformat in ["svg", "pdf", "html"]
            plot_fname = joinpath(heatmaps_path,
                              "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_contrast$(signif ? "_signif" : "")_heatmap.$(outformat)")
        try
                savefig(coll_heatmap, plot_fname, format=outformat, width=coll_heatmap.plot.layout[:width], height=coll_heatmap.plot.layout[:height]);
        catch e
            if e isa InterruptException
                rethrow(e)
            else
                @warn "$term_coll generation failed: $e"
            end
        end
    end
end
end