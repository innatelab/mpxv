proj_info = (id = "mpxv",
            modelobj = "protregroup",
            fit_ver_ivip = "20221130",
            fit_ver_mpxv = "20221201",
            data_ver_ivip = "20221129",
            data_ver_mpxv = "20221130",
            oesc_ver = "20221205"
)

const base_scripts_path = "/home/ge54heq/projects"
const base_analysis_path = "/pool/analysis/yhuang"

using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, DataFrames
using StatsBase

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const reports_path = joinpath(analysis_path, "reports")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")
const party3rd_data_path = "/pool/pub3rdparty"

includet(joinpath(misc_scripts_path, "frame_utils.jl"));
includet(joinpath(misc_scripts_path, "msglm_utils.jl"));
includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));

objid_col = Symbol(string(proj_info.modelobj, "_id"));

input_rdata_ivip = load(joinpath(scratch_path, "$(proj_info.id)_ivip_msglm_data_$(proj_info.data_ver_ivip).RData"), convert=true)
full_rdata_ivip = load(joinpath(scratch_path, "$(proj_info.id)_ivip_msdata_full_$(proj_info.data_ver_ivip).RData"), convert=true)
fit_rdata_ivip = load(joinpath(results_path, "$(proj_info.id)_ivip_msglm_fit_$(proj_info.fit_ver_ivip).RData"), convert=true)

input_rdata_mpxv = load(joinpath(scratch_path, "$(proj_info.id)_HeLa_msglm_data_$(proj_info.data_ver_mpxv).RData"), convert=true)
full_rdata_mpxv = load(joinpath(scratch_path, "$(proj_info.id)_HeLa_msdata_full_$(proj_info.data_ver_mpxv).RData"), convert=true)
fit_rdata_mpxv = load(joinpath(results_path, "$(proj_info.id)_HeLa_msglm_fit_$(proj_info.fit_ver_mpxv).RData"), convert=true)

obj_contrasts_ivip_df = fit_rdata_ivip["object_contrasts.df"]
obj_contrasts_mpxv_df = fit_rdata_mpxv["object_contrasts.df"]
obj_contrasts_df = vcat(obj_contrasts_ivip_df, obj_contrasts_mpxv_df)|> MSGLMUtils.fix_object_id!

contrast_cols = [:contrast, :timepoint_lhs, :treatment_lhs]
contrasts_df = unique!(select(obj_contrasts_df, [contrast_cols; :change]))
contrasts_df.contrast_alt = contrasts_df.contrast
contrasts_df.contrast = string.(contrasts_df.contrast_alt, contrasts_df.change) 

objects_ivip_df = copy(input_rdata_ivip["msdata"][string(proj_info.modelobj, "s")]) |> MSGLMUtils.fix_object_id!;
objects_mpxv_df = copy(input_rdata_mpxv["msdata"][string(proj_info.modelobj, "s")]) |> MSGLMUtils.fix_object_id!;
cols4objects = [:protregroup_id, :majority_protein_acs, :protein_acs, :protac_label, :protregroup_label, :object_id]
select!(objects_ivip_df, cols4objects)
select!(objects_mpxv_df, cols4objects)
objects_df = outerjoin(objects_ivip_df, objects_mpxv_df, on = cols4objects)|> MSGLMUtils.fix_object_id!;

obj2protac_df = select!(
         filter(r -> r.is_majority, full_rdata_ivip["msdata_full"][string("protein2", proj_info.modelobj)]),
         [objid_col, :protein_ac]) |> unique! |> MSGLMUtils.fix_object_id! #this dataframe is identical from both datasets as they were quanted together

includet(joinpath(misc_scripts_path, "gmt_reader.jl"));
include(joinpath(misc_scripts_path, "optcover_utils.jl"));
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

@info "Preparing mosaics..."
observed_protacs = Set(obj2protac_df.protein_ac) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls, 
                                    protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);
                                      
ENV["MKL_NUM_THREADS"] = 1
using OptEnrichedSetCover
cover_params = CoverParams(setXset_factor=0.5,
                           uncovered_factor=0.0, covered_factor=0.0)

# GO enrichment for contrasts, all vs mock comparisons
obj_contrast_hit_sets = OmicsCollections.effects_collection(
    filter(r -> r.ci_target == "average", obj_contrasts_df),
    obj_col=:object_id, change_col=:median, sel_col=:is_hit,
    group_cols=[:var], #don't skip the group_cols even though it's not always meaningful
    effect_col=:contrast)    

obj_contrast_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by contrast hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_contrast_hit_sets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

obj_contrast_hit_mosaics_v = collect(pairs(obj_contrast_hit_mosaics))
obj_contrast_hit_covers_v = similar(obj_contrast_hit_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_contrast_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_contrast_hit_mosaics_v[i]
    @info "Covering $mosaic_name by contrast hits..."
    obj_contrast_hit_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.1, 0.1], MaxSteps=2_000_000, WeightDigits=2, NWorkers=1),
            true)
end
obj_contrast_hit_covers = Dict(k => v for (k, v) in obj_contrast_hit_covers_v)

# GO enrichment for contrasts, only at 24h
obj_contrast_24h_hit_sets = OmicsCollections.effects_collection(
    filter(r -> r.ci_target == "average" && r.timepoint_lhs == "24", obj_contrasts_df),
    obj_col=:object_id, change_col=:median, sel_col=:is_hit,
    group_cols=[:var], 
    effect_col=:contrast)    

obj_contrast_24h_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by contrast hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_contrast_24h_hit_sets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

obj_contrast_24h_hit_mosaics_v = collect(pairs(obj_contrast_24h_hit_mosaics))
obj_contrast_24h_hit_covers_v = similar(obj_contrast_24h_hit_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_contrast_24h_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_contrast_24h_hit_mosaics_v[i]
    @info "Covering $mosaic_name by contrast hits..."
    obj_contrast_24h_hit_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.1, 0.1], MaxSteps=2_000_000, WeightDigits=2, NWorkers=1),
            true)
end
obj_contrast_24h_hit_covers = Dict(k => v for (k, v) in obj_contrast_24h_hit_covers_v)


using JLD2

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_HeLa_hit_$(proj_info.oesc_ver)_covers.jld2")
@save(hit_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      objects_df, obj_contrasts_df,
      contrasts_df,
      obj_contrast_hit_sets, obj_contrast_24h_hit_sets, obj_contrast_hit_mosaics, obj_contrast_24h_hit_mosaics,
      cover_params, obj_contrast_hit_covers, obj_contrast_24h_hit_covers)
if !@isdefined(obj_hit_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
        proj_info, protac_colls, obj_colls, obj_mosaics,
        obj2term_df, terms_df,
        objects_df, obj_contrasts_df,
        contrasts_df,
        obj_contrast_hit_sets, obj_contrast_24h_hit_sets, obj_contrast_hit_mosaics, obj_contrast_24h_hit_mosaics,
        cover_params, obj_contrast_hit_covers, obj_contrast_24h_hit_covers)
end

include(joinpath(misc_scripts_path, "optcover_utils.jl"));
@info "Preparing protgroup->gene_name map..."
obj_id2name = Dict(r.object_id => r[Symbol(proj_info.modelobj, "_label")]
                   for r in eachrow(objects_df))

obj_contrast_hit_covers_df = innerjoin(
    OptCoverUtils.covers_report(
    obj_contrast_hit_covers, obj_contrast_hit_sets, 
    obj_colls, 
    obj_id2name, terms_df,
    cover_params = cover_params,
    experimentid_col=[:var, :contrast], weightedset_col_prefix="hit"),
    contrasts_df, on = [:contrast]);
obj_contrast_hit_covers_df.intersect_genes = [join(unique(vcat(split.(split(genes, ' '), ';')...)), ' ') for genes in obj_contrast_hit_covers_df.intersect_genes]

# don't remove the sets since they are timecourses timepoints
obj_contrast_hit_covers_signif_df = combine(groupby(obj_contrast_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    df = OptCoverUtils.filter_multicover(coll_df, set_cols=[:contrast],
                                         max_term_pvalue=1E-3, max_set_pvalue=nothing, min_set_overlap=nothing)
    return select!(df, Not(:term_collection))
end

obj_contrast_24h_hit_covers_df = innerjoin(
    OptCoverUtils.covers_report(
    obj_contrast_24h_hit_covers, obj_contrast_24h_hit_sets, 
    obj_colls, 
    obj_id2name, terms_df,
    cover_params = cover_params,
    experimentid_col=[:var, :contrast], weightedset_col_prefix="hit"),
    contrasts_df, on = [:contrast]);
obj_contrast_24h_hit_covers_df.intersect_genes = [join(unique(vcat(split.(split(genes, ' '), ';')...)), ' ') for genes in obj_contrast_24h_hit_covers_df.intersect_genes]

# don't remove the sets since they are timecourses timepoints
obj_contrast_24h_hit_covers_signif_df = combine(groupby(obj_contrast_24h_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    df = OptCoverUtils.filter_multicover(coll_df, set_cols=[:contrast],
                                         max_term_pvalue=1E-3, max_set_pvalue=nothing, min_set_overlap=nothing)
    return select!(df, Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_HeLa_united_hit_covers_$(proj_info.oesc_ver).txt"),
          obj_contrast_hit_covers_df[obj_contrast_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_HeLa_united_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_contrast_hit_covers_signif_df[obj_contrast_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_HeLa24h_united_hit_covers_$(proj_info.oesc_ver).txt"),
          obj_contrast_24h_hit_covers_df[obj_contrast_24h_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_HeLa24h_united_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_contrast_24h_hit_covers_signif_df[obj_contrast_24h_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');


Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_plots.jl"))
include(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

oesc_plots_path = joinpath(plots_path, "oesc_HeLa_$(proj_info.oesc_ver)")
isdir(oesc_plots_path) || mkdir(oesc_plots_path)

using PlotlyJS, TextWrap#, ORCA

heatmap_layout_attrs = Dict(
    ("SigDB_C2", true) => Dict(:margin_l => 500),
    ("SigDB_C2", false) => Dict(:margin_l => 500),
    #("Reactome", true) => Dict(:margin_l => 500),
    #("Reactome", false) => Dict(:margin_l => 500),
    ("Reactome", true) => Dict(:margin_l => 600),
    ("Reactome", false) => Dict(:margin_l => 600),
    ("GO_CC", true) => Dict(:margin_l => 200),
    ("GO_CC", false) => Dict(:margin_l => 200),
)

treatment_info = Dict("mock" => (color="gray", label = "Mock"),
                      "MVA_F" => (color="#2DA8FF", label = "MVA"),
                      "MVA_dE9L" => (color="#2D3FFF", label = "MVA_dE3L"),
                      "CVA_152" => (color="#83e000", label = "CVA"),
                      "VACV_WR" => (color="#f74a0f", label = "VACV_WR"),
                      "MPXV" => (color="#F9CB40", label = "MPXV"),
                      )

function stylize_contrast_multi(contrast::AbstractString, condition_info::AbstractDict{<:AbstractString})
    contrast_match = match(r"^(?<lhs>.+)_(?<time>\d+h)_VS_(?<rhs>[^+-]+)_(\d+)h(?<sign>[+-])?$", contrast)
    if (isnothing(contrast_match))
        @warn "No format match for contrast '$(contrast)'"
        return contrast
    else
        res = OptCoverHeatmap.stylize_treatment(contrast_match[:lhs], condition_info) * " <b>vs</b> " *
        OptCoverHeatmap.stylize_treatment(contrast_match[:rhs], condition_info);
        if (!isnothing(contrast_match[:time]))
            res *= " <b>@</b> " * contrast_match[:time];
        end
        if (!isnothing(contrast_match[:sign]))
            res *= OptCoverHeatmap.stylize_change(contrast_match[:sign]);
        end    
        return res
    end
end

function process_contrast_axis(contrasts_df)
    contrast_labels = stylize_contrast_multi.(contrasts_df.contrast, Ref(treatment_info))

    contrasts_df,
    contrast_labels,
    string.("contrast: ", contrast_labels)
end

treatment_order = Dict(:MPXV => 1, :VACV_WR => 2, :CVA_152 => 3, :MVA_F => 4, :MVA_dE9L => 5)

for term_coll in unique(obj_contrast_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")contrast heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? obj_contrast_hit_covers_signif_df : obj_contrast_hit_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="protein",
            experiment_axis_title = "contrast",
            experiment_cols = [:contrast, :treatment_lhs, :timepoint_lhs, :change, :nhit],
            process_experiment_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 80),
            cell_width=40, cell_height=30,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",
            experiment_order=contrasts -> begin
                contrasts.treatment_order = [treatment_order[Symbol(r.treatment_lhs)]
                                            for r in eachrow(contrasts)]
                return sortperm(contrasts, [:change, :treatment_order, :timepoint_lhs])
            end)
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>80,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    for outformat in ["svg", "pdf","html"]
        plot_fname = joinpath(oesc_plots_path,
            "$(proj_info.id)_HeLa_$(proj_info.oesc_ver)_$(term_coll)_contrasts$(signif ? "_signif" : "")_heatmap.$(outformat)")
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


for term_coll in unique(obj_contrast_24h_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")contrast heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? obj_contrast_24h_hit_covers_signif_df : obj_contrast_24h_hit_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="protein",
            experiment_axis_title = "contrast",
            experiment_cols = [:contrast, :treatment_lhs, :timepoint_lhs, :change, :nhit],
            process_experiment_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 80),
            cell_width=40, cell_height=30,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",
            experiment_order=contrasts -> begin
                contrasts.treatment_order = [treatment_order[Symbol(r.treatment_lhs)]
                                            for r in eachrow(contrasts)]
                return sortperm(contrasts, [:change, :treatment_order, :timepoint_lhs])
            end)
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>80,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    for outformat in ["svg", "pdf","html"]
        plot_fname = joinpath(oesc_plots_path,
            "$(proj_info.id)_HeLa_$(proj_info.oesc_ver)_$(term_coll)_24h_contrasts$(signif ? "_signif" : "")_heatmap.$(outformat)")
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
