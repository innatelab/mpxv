proj_info = (id = "mpxv",
             data_ver = "20220812",
             fit_ver = "20220813",
             oesc_ver = "20220829",
             modelobj = "ptmngroup",
             countobj= "ptmngroup",
             mstype = "phospho",
             msfolder = "phospho_20220812")

using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, CSV, CodecZlib, DataFrames, FastaIO
using JLD2
using StatsBase

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const base_scripts_path = "/home/ge54heq/projects"
const base_analysis_path = "/pool/analysis/yhuang"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "msglm_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));

objid_col = Symbol(string(proj_info.modelobj, "_id"));

input_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(proj_info.mstype)_$(proj_info.fit_ver).RData"), convert=true)
full_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(proj_info.mstype)_$(proj_info.data_ver).RData"), convert=true)
fit_rdata = load(joinpath(results_path, "$(proj_info.id)_msglm_fit_$(proj_info.mstype)_$(proj_info.fit_ver).RData"), convert=true)
ptmns_df = copy(input_rdata["msdata"]["ptmngroups"])
ptm2group_df = select(ptmns_df, [:ptm_id, :ptmngroup_id]) |> unique!

obj_contrasts_df = copy(fit_rdata["object_contrasts.df"]);
contrast_cols = [:contrast, :timepoint_lhs, :timepoint_rhs, :treatment_lhs, :treatment_rhs]
contrasts_df = unique!(select(obj_contrasts_df, [contrast_cols; :change]))
contrasts_df.timepoint_lhs = Vector{Union{String}}(contrasts_df[!,:timepoint_lhs])
contrasts_df.timepoint_lhs = parse.(Int, contrasts_df.timepoint_lhs)
contrasts_df.timepoint_rhs = Vector{Union{String}}(contrasts_df[!,:timepoint_rhs])
contrasts_df.timepoint_rhs = parse.(Int, contrasts_df.timepoint_rhs)

objects_df = copy(input_rdata["msdata"][string(proj_info.modelobj, "s")]);

for df in [obj_contrasts_df]
    df |> MSGLMUtils.fix_quantile_columns!
end
# pick the most significant comparison/effect
#=
obj_contrasts_df.rowix = 1:nrow(obj_contrasts_df)
objcontr_selixs = Vector{Int}()
for df in groupby(obj_contrasts_df, [:contrast, :ci_target, objid_col])
    selix = findmin(df.p_value)[2]
    push!(objcontr_selixs, df.rowix[selix])
end
obj_contrasts_df = obj_contrasts_df[objcontr_selixs, :]
obj_contrasts_df.object_id = copy(obj_contrasts_df[!, objid_col])
=#
ptm2gene_df = copy(full_rdata["msdata_full"]["ptm2gene"])
ptm2gene_df.ptm_id = convert(Vector{Int}, ptm2gene_df.ptm_id)
ptm2gene_df.ptm_pos = convert(Vector{Int}, ptm2gene_df.ptm_pos)
ptm2gene_df.ptm_site = ["$(r.ptm_AA_seq)$(r.ptm_pos)" for r in eachrow(ptm2gene_df)]

include(joinpath(misc_scripts_path, "optcover_utils.jl"));

@load(joinpath(scratch_path, "phosphositeplus_annotations_20220808.jld2"), psitep_annot_dfs, psitep_annots_df)
for df in values(psitep_annot_dfs)
    df.ptm_id = df.protein_ac .* "_" .* uppercase.(df.ptm_AA) .* string.(df.ptm_pos)
end
psitep_annots_df = let kinsub_df = copy(psitep_annot_dfs[:KinaseSubstrate])
    filter!(r -> r.organism=="human", dropmissing!(kinsub_df, :kinase_gene_name))
    kinsub_df[!, :coll_id] .= :KinaseSubstrate
    kinsub_df.term_id = uppercase.(kinsub_df.kinase_gene_name)
    kinsub_df
end

#=
perseus_report_df = CSV.read(joinpath(data_path, proj_info.msfolder,
        "cov2_snaut_parsars_phospho_20201005_contrasts_report_20201012_long_psp.txt"),
        header=1, datarow=3, comment="#", delim='\t')
perseus_motifs_df = dropmissing!(DelimDataUtils.expand_delim_column(select(perseus_report_df, [:protein_ac, :ptm_site, :Motifs]) |> unique!,
                                 list_col=:Motifs, elem_col=:term_id, key_col=[:protein_ac, :ptm_site]),
                                 [:term_id, :protein_ac])
perseus_motifs_df.term_id = replace.(perseus_motifs_df.term_id, Ref(r"\s+(?:motif|sequence)$" => ""))
perseus_motifs_df.ptm_AA = [r.ptm_site[1:1] for r in eachrow(perseus_motifs_df)]
perseus_motifs_df.ptm_pos = [parse(Int, r.ptm_site[2:end]) for r in eachrow(perseus_motifs_df)]
perseus_motifs_df[!, :coll_id] .= :Motifs

motif_colls = FrameUtils.frame2collections(ptmgroup2motif_df, set_col=:term_id, obj_col=:ptmgroup_id, coll_col=:coll_id)

motifXregulators_collapsed_df = CSV.read(joinpath(data_path, proj_info.msfolder, "Perseus_motifsXregulators.txt"))
motifXregulators_df = DelimDataUtils.expand_delim_column(motifXregulators_collapsed_df,
                            list_col=:regulators, elem_col=:regulator, key_col=:motif, delim=",")

function motif2regulators(df::DataFrame; keycol::Union{Symbol, AbstractVector}=:ptmgroup_id)
    dfnew = rename!(select!(join(df, motifXregulators_df,
                         on=[:term_id=>:motif], kind=:left), Not(:term_id)),
            :regulator => :term_id) |> unique!
    dropmissing!(dfnew, [:term_id; keycol])
    dfnew[!, :coll_id] .= :Regulators
    return dfnew
end
perseus_regulators_df = motif2regulators(perseus_motifs_df, keycol=[:protein_ac, :ptm_AA, :ptm_pos])

# https://linkphinder.insight-centre.org/
linkphinder_df = CSV.read(joinpath(party3rd_data_path, "linkPhinder_data.csv"))
countmap(linkphinder_df.KinaseLabel)
linkphinder_annots_df = rename!(filter(r -> r.Score >= 0.75, linkphinder_df), :KinaseLabel => :term_id, :ProteinSubstrate_ID => :protein_ac)
linkphinder_annots_df.ptm_pos = [parse(Int, r.Site[2:end]) for r in eachrow(linkphinder_annots_df)]
linkphinder_annots_df.ptm_AA = getindex.(linkphinder_annots_df.Site, Ref(1:1))
linkphinder_annots_df.coll_id = :LinkPhinder
=#

ptm_annot_cols = [:coll_id, :term_id, :protein_ac, :ptm_AA, :ptm_pos]
ptm_annots_df = select(psitep_annots_df, ptm_annot_cols, copycols=false)

ptmgroup_annots_df = dropmissing!(select!(leftjoin(leftjoin(ptm_annots_df,
            rename!(select(ptm2gene_df, [:ptm_id, :protein_ac, :ptm_AA_seq, :ptm_pos], copycols=false), :ptm_AA_seq => :ptm_AA),
            on=[:protein_ac, :ptm_AA, :ptm_pos]),
            ptm2group_df, on=:ptm_id, matchmissing = :notequal),
            [:ptmngroup_id, :coll_id, :term_id]), [:ptmngroup_id, :term_id]) |> unique!
countmap(ptmgroup_annots_df.coll_id)

obj_colls = FrameUtils.frame2collections(ptmgroup_annots_df,
    set_col=:term_id, obj_col=:ptmngroup_id, coll_col=:coll_id)

terms_df = unique!(select(ptmgroup_annots_df, [:coll_id, :term_id]));
terms_df.term_name = copy(terms_df.term_id)
terms_df[!, :term_descr] .= missing

@info "Preparing sets of hits"
sel_ci_target = "average"
ObjectType = eltype(ptmgroup_annots_df.ptmngroup_id)
obj_hit_sets = Dict{Tuple{String, String}, Set{ObjectType}}()
for hit_df in groupby(filter(r -> coalesce(r.is_hit_nofp, false) && (r.ci_target == sel_ci_target) && !r.is_contaminant, obj_contrasts_df), [:contrast, :change])
    obj_hit_sets[(hit_df[1, :contrast], hit_df[1, :change])] = Set(skipmissing(hit_df.ptmngroup_id))
end
# only relevant ones
sel_contrasts_df = filter(r -> r.treatment_rhs == "mock", contrasts_df)
obj_hit_selsets = filter(kv -> kv[1][1] ∈ sel_contrasts_df.contrast, obj_hit_sets)

@info "Preparing mosaics..."
observed_ptms = Set(skipmissing(semijoin(obj_contrasts_df, ptmgroup_annots_df, on=:ptmngroup_id).ptmngroup_id))
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls, setXset_frac_extra_elms=0.05, verbose=true);

obj_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_hit_selsets,
                                          max_sets=2000, min_nmasked=1, max_setsize=2000, verbose=true)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.1,
                           uncovered_factor=0.0, covered_factor=0.0)#, covered_factor=0.002)

obj_hit_mosaics_v = collect(pairs(obj_hit_mosaics))
obj_hit_covers_v = similar(obj_hit_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_hit_mosaics_v[i]
    @info "Covering $mosaic_name by hits..."
    obj_hit_covers_v[i] = mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.02], MaxSteps=2_000_000, WeightDigits=2,
                                    NWorkers=1, #Threads.nthreads()-1,
                                    MaxRestarts=200),
            true)
end
obj_hit_covers = Dict(k => v for (k, v) in obj_hit_covers_v)

using JLD2

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_$(proj_info.oesc_ver).jld2")
@save(hit_covers_filename,
      proj_info, obj_colls, obj_mosaics,
      terms_df, ptm_annots_df, ptmgroup_annots_df, contrasts_df,
      ptmns_df, obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
if !@isdefined(obj_effect_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
      proj_info, obj_colls, obj_mosaics,
      terms_df, ptm_annots_df, ptmgroup_annots_df, contrasts_df,
      ptmns_df, obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
end

include(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing protgroup?gene_name map..."
obj_id2name = Dict(r.ptmngroup_id => replace(r.ptmngroup_label_no_ptm_type, r"_M\d$" => "")
                   for r in eachrow(ptmns_df))

obj_hit_covers_df = innerjoin(OptCoverUtils.covers_report(
    obj_hit_covers, obj_hit_selsets, obj_colls, #obj_mosaics, 
    obj_id2name, terms_df,
    experimentid_col=[:contrast, :change],
    weightedset_col_prefix="hit"),
    contrasts_df, on=[:contrast, :change])
filter!(r -> r.term_id != "CHEK2" && r.term_id != "CK1E", obj_hit_covers_df) # redundant with Chk2 and CSK1E

obj_hit_covers_signif_df = combine(groupby(obj_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:contrast, :change],
                                                   max_term_pvalue=1E-2, max_set_pvalue=nothing, min_set_overlap=nothing),
                   Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_$(proj_info.oesc_ver).txt"),
          filter(r -> r.nmasked > 0, obj_hit_covers_df),
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          filter(r -> r.nmasked > 0, obj_hit_covers_signif_df),
          missingstring="", delim='\t');

includet(joinpath(misc_scripts_path, "frame_utils.jl"))
includet(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

using PlotlyJS, TextWrap

heatmap_layout_attrs = Dict(
    (:Motifs, true) => Dict(:margin_l => 300, :margin_b => 200),
    (:Motifs, false) => Dict(:margin_l => 300, :margin_b => 200),
    (:Regulators, true) => Dict(:margin_l => 150, :margin_b => 200),
    (:Regulators, false) => Dict(:margin_l => 150, :margin_b => 200),
)

stylize_contrast(str) = foldl(replace, [
    r"(MPXV?)_vs_mock@(\d+)h" => s"\1:<span style=\"font-weight: bold; color: black;\">\2</span>h",
    "MPXV:" => "<span style=\"font-wieght: bold; color: #F9CB40;\">MPXV</span> ",
    ],
    init = str)

function process_contrast_axis(contrast_df)
    contrast_df,
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change),
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change)#stylize_effect.(effect_df.effect)
end

heatmaps_path = joinpath(plots_path, "$(proj_info.msfolder)_$(proj_info.fit_ver)", "hits_oesc_$(sel_ci_target)_$(proj_info.oesc_ver)")
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
            elements_label="PTM",
            experiment_axis_title = "contrast",
            experiment_cols = [:contrast, :treatment_lhs, :timepoint_lhs, :change, :nhit],
            process_experiment_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 320),
            margin_b=get(layout_attrs, :margin_b, 400),
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",#outformat in ["svg", "pdf"] ? "#000" : "#BBB",
            zmin=-10, cell_width=50, cell_height=30, gridwidth=2,
            experiment_order=contrasts -> begin
                return reverse(sortperm(contrasts, [order(:change, rev=true), :timepoint_lhs, :treatment_lhs]))
            end,
            transpose=true)
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>60,
                   #:margin_t=>20,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    plot_fname = joinpath(heatmaps_path,
                          "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_X_contrast$(signif ? "_signif" : "")_heatmap.$(outformat)")
    for outformat in ["svg", "pdf", "html"]
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
