# Map the modified peptides back to the correct positions in protein(groups)
#
# Author: Yiqi Huang
###############################################################################

proj_info = (id = "mpxv",
             data_ver = "20221105",
             msfolder = "phospho_20221105",
             ptm_locprob_min = 0.75,
             )
using Pkg
Pkg.activate(@__DIR__)
using Revise
using StatsBase, DataFrames, CSV, FastaIO, JLD2, CodecZlib, Base.Filesystem, CategoricalArrays

const base_scripts_path = "/home/ge54heq/projects"
const base_analysis_path = "/pool/analysis/yhuang"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl")
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data", proj_info.msfolder)
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")
const party3rd_data_path = "/pool/pub3rdparty"

includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));
includet(joinpath(misc_scripts_path, "frame_utils.jl"));
include(joinpath(misc_scripts_path, "ms_import.jl")); # includet doesn't work so well with another include
includet(joinpath(misc_scripts_path, "fasta_reader.jl"));
includet(joinpath(misc_scripts_path, "ptm_extractor.jl"));
includet(joinpath(misc_scripts_path, "phosphositeplus_reader.jl"));

@info "Add more info to msruns..."
msruns_df = CSV.read(joinpath(data_path, "combined", "experimentalDesign.txt"), DataFrame) #check that the file is UTF8 encoded!
rename!(msruns_df, :Name=>:rawfile, :Experiment=>:msexperiment, :PTM => :is_ptm, :Fraction => :msfraction)
msrun_matches = match.(Ref(r"(.*)HFF_(\d+)h_(.*)_(\d)$"), msruns_df.msexperiment)
msruns_df.treatment = getindex.(msrun_matches, 3);
msruns_df.treatment = replace.(msruns_df.treatment, "mpx" => "MPXV") #instead of MPX, use the full name MPXV!
msruns_df.timepoint = getindex.(msrun_matches, 2);
msruns_df.condition = string.(msruns_df.treatment, "_", msruns_df.timepoint);
msruns_df.replicate = parse.(Int, getindex.(msrun_matches, 4));
msruns_df.msexperiment = string.(msruns_df.condition, "_", msruns_df.replicate);
msruns_df.is_skipped = contains.((getindex.(msrun_matches, 1)), "MBR")

# read direct MaxQuant output
@info "Reading phospho report..."
pmsreport = let
    report_df, colgroups = MSImport.MaxQuant.read_evidence(joinpath(data_path, "combined", "txt"), #limit_rows=1000,
            msruns=msruns_df,
            import_data=[:peptide, :pepmodstate, :pepmod_ptms], delim='\t', verbose=true)
    (data = report_df, colgroups = colgroups)
end #if it doesn't work, check if you have CSV and DataFrame in your Manifest.toml! Add them if not.

msruns_df_new = unique( pmsreport.data[!, last.(pmsreport.colgroups[:mschannel])])

# fix contaminant names (strip "CON__<SRC>:")
for col in [:lead_protein_acs, :lead_razor_protein_ac, :peptide_protein_acs]
    df = pmsreport.data
    mask = .!ismissing.(df[!, col])
    df[mask, col] .= replace.(replace.(replace.(df[mask, col],
                              Ref(r"CON__(?:[^:;]+):CON__" => "CON__")),
                              Ref(r"CON__(?:[^:;]+):" => "CON__")),
                              Ref(r"CON__Streptavidin" => "CON__P22629"))
end

msfasta_path = joinpath(data_path, "fasta")
proteins_df = let
    human_df = Fasta.read_uniprot(joinpath(msfasta_path, "uniprot-reviewed_yes+taxonomy_9606.fasta"))
    human_df[!, :is_contaminant] .= false
    human_df[!, :is_viral] .= false

    mpxv_df = Fasta.read_uniprot(joinpath(msfasta_path, "MPXV_reformatted_20220811.fasta"))
    mpxv_df.genename = "MPXV_" .* mpxv_df.genename
    mpxv_df[!, :is_contaminant] .= false
    mpxv_df[!, :is_viral] .= true

    contaminants_df = Fasta.read_contaminants(joinpath(msfasta_path, "contaminants_20200405.fasta"))
    contaminants_df[!, :is_contaminant] .= true
    contaminants_df[!, :is_viral] .= false
    contaminants_df[!, :protein_code] .= missing
    contaminants_df[!, :genename] .= missing
    contaminants_df[!, :protein_existence] .= missing
    contaminants_df[!, :seq_version] .= missing

    res = vcat(human_df, mpxv_df, contaminants_df)
    res.orig_genename = res.genename
    res.genename = coalesce.(res.orig_genename, res.protein_ac)
    res.protein_ac_isoform = [(isnothing(m) ? 1 : parse(Int, m[1])) for m in match.(Ref(r"-(\d+)(?:#.+)?$"), res.protein_ac)]
    res
end
protein_seqs = Dict(r.protein_ac => r.seq for r in eachrow(proteins_df))

# find all modified_peptide -to- sequence matches
peptides_df = sort!(unique(pmsreport.data[!, last.(pmsreport.colgroups[:peptide])]), :peptide_id)
#peptides_df.peptide_id = FrameUtils.indexunique(peptide2protgroups_df.peptide_seq)
# global pepmods
pepmodstates_df = sort!(unique(pmsreport.data[!, last.(pmsreport.colgroups[:pepmodstate])]), #"peptide_id" is already in the pepmodstate colgroup
                        [:pepmodstate_id, :peptide_id, :pepmod_id, :charge]) #remove "msfraction" because there's only one
@assert pepmodstates_df.pepmodstate_id == 0:nrow(pepmodstates_df)-1

peptide2protein_df = PTMExtractor.match_peptides(filter(r -> !ismissing(r.peptide_protein_acs), peptides_df), protein_seqs)

ptm2pms_df, pms_ptms_stats_df = PTMExtractor.extract_ptms(pepmodstates_df, format=:maxquant, objid_col=:pepmodstate_id, selptms=r"Phospho|GlyGly")
filter!(r -> !in(r.ptm_AAs, ["Protein N-term", "Protein C-term"]), ptm2pms_df)
ptm2pms_df = innerjoin(unique!(ptm2pms_df), select(pepmodstates_df, [:pepmodstate_id, :peptide_id], copycols=false),
                       on=:pepmodstate_id) |> unique!
sort!(select!(ptm2pms_df, Not(:peptide_seq)), [:pepmodstate_id, :ptm_offset])
pepmodstates_df = innerjoin(pepmodstates_df, pms_ptms_stats_df, on=:pepmodstate_id)
ptm2pms2protein_df = innerjoin(ptm2pms_df, peptide2protein_df, on=:peptide_id)
ptm2pms2protein_df.ptm_pos = ptm2pms2protein_df.peptide_pos .+ ptm2pms2protein_df.ptm_offset
ptm2protein_df = leftjoin(unique!(select(ptm2pms2protein_df, [:ptm_type, :ptm_AAs, :ptm_AA_seq, :ptm_pos, :protein_ac])),
                          select(proteins_df, [:protein_ac, :genename, :organism, :is_viral, :is_contaminant], copycols=false),
                          on=:protein_ac)
ptm2gene_df = PTMExtractor.group_aaobjs(ptm2protein_df, proteins_df,
                                        seqgroup_col=[:genename, :is_viral, :is_contaminant],
                                        seqid_col=:protein_ac, seqrank_col=:protein_ac_isoform,
                                        obj_prefix=:ptm_, objid_col=[:ptm_type, :ptm_pos, :ptm_AA_seq],
                                        force_refseqs=true, verbose=true) #Note: up till now the BioAlginments still hasn't included Alexey's update (seq2aln), please add the package with the git URL https://github.com/BioJulia/BioAlignments.jl.git instead for this part to work!
ptm2gene_df.flanking_15AAs = PTMExtractor.flanking_sequence.(getindex.(Ref(protein_seqs), ptm2gene_df.protein_ac), ptm2gene_df.ptm_pos, flanklen=15)
ptm2protein_df2 = leftjoin(ptm2protein_df, unique!(select(ptm2gene_df, [:ptm_id, :protein_ac, :ptm_pos, :ptm_AA_seq])),
                           on=[:protein_ac, :ptm_pos, :ptm_AA_seq])
@assert nrow(ptm2protein_df)==nrow(ptm2protein_df2)
ptm2protein_df = ptm2protein_df2

filter(r -> r.ptm_is_reference && coalesce(r.is_viral, false), ptm2gene_df) |> print
countmap(filter(r -> r.ptm_is_reference, ptm2gene_df).ptm_type)

ptmn2pms_df = innerjoin(select!(innerjoin(
        select(ptm2gene_df, [:ptm_id, :ptm_label, :protein_ac, :ptm_type, :ptm_AA_seq, :ptm_pos], copycols=false),
        ptm2pms2protein_df, on=[:protein_ac, :ptm_type, :ptm_AA_seq, :ptm_pos]),
        [:ptm_id, :ptm_label, :ptm_type, :peptide_id, :pepmodstate_id, :ptm_offset]) |> unique!,
        select(pepmodstates_df, [:pepmodstate_id, :nselptms], copycols=false), on=:pepmodstate_id)
ptmns_df = unique!(select(ptmn2pms_df, [:ptm_id, :ptm_label, :ptm_type, :nselptms]))
sort!(ptmns_df, [:ptm_type, :ptm_id, :nselptms])
ptmns_df.ptmn_id = 1:nrow(ptmns_df)
ptmns_df.ptmn_label = categorical(string.(ptmns_df.ptm_label) .* "_M" .* string.(ptmns_df.nselptms))
ptmn2pms_df = innerjoin(ptmn2pms_df, select(ptmns_df, Not(:ptm_label), copycols=false), on=[:ptm_id, :ptm_type, :nselptms])
sort!(ptmn2pms_df, [:ptmn_id, :pepmodstate_id])

ptm2ptmgroup_df, ptmn2ptmngroup_df, ptmgroup2ptmngroup_df = PTMExtractor.ptmgroups(ptmn2pms_df, verbose=true)
ptmngroups_df = select(ptmn2ptmngroup_df, [:ptmngroup_id, :ptm_type, :nselptms])|> unique!
ptmngroup2pms_df = leftjoin(ptmn2ptmngroup_df, select(ptmn2pms_df, [:pepmodstate_id, :nselptms, :ptmn_id, :ptmn_label]),on = [:ptmn_id, :nselptms])

print(countmap(collect(values(countmap(ptmn2pms_df.pepmodstate_id)))))
print(countmap(collect(values(countmap(filter(r -> r.nselptms == 1, ptmns_df).ptmn_id)))))

pms_intensities_df = sort!(select(pmsreport.data, [:evidence_id, :pepmodstate_id, :msrun, :intensity, :ident_type, :psm_pvalue]), :evidence_id)
filter!(r -> !ismissing(r.intensity), pms_intensities_df)

pms_locprobs_wide_df = select(pmsreport.data, [["pepmodstate_id", "pepmod_seq", "peptide_id", "msrun", "evidence_id"]; last.(pmsreport.colgroups[:pepmod_ptms])])
pms_locprobs_df = MSImport.MaxQuant.pivot_ptms_longer(pms_locprobs_wide_df, [:evidence_id, :pepmodstate_id, :peptide_id, :pepmod_seq, :msrun])
filter!(r -> !ismissing(r.ptm_locprob_seq), pms_locprobs_df)
disallowmissing!(pms_locprobs_df, [:ptm_locprob_seq, :ptm_scorediff_seq])

ptm2pms_locprobs_df = PTMExtractor.extract_ptm_locprobs(pms_locprobs_df, format=:maxquant, modseq_col=:pepmod_seq, msrun_col=:msrun)
ptm2pms_locprobs_df.evidence_id = pms_locprobs_df.evidence_id[ptm2pms_locprobs_df.report_ix]
ptm2pms_locprobs_df.peptide_id = pms_locprobs_df.peptide_id[ptm2pms_locprobs_df.report_ix]
FrameUtils.matchcategoricals!(ptm2pms_locprobs_df, ptmns_df)
ptmn_locprobs_df = innerjoin(ptm2pms_locprobs_df, ptmn2pms_df,
                             on=intersect(names(ptm2pms_locprobs_df), names(ptmn2pms_df)))
unique!(select!(ptmn_locprobs_df, Not([:report_ix, :ptm_offset])))

countmap(collect(coalesce.(ptmn_locprobs_df.ptm_locprob, 0.0) .>= proj_info.ptm_locprob_min))

# save the files
output_path = joinpath(data_path, "ptm_extractor_$(proj_info.data_ver)")
isdir(output_path) || mkdir(output_path)
for (fname, df) in ["rawfiles_info.txt" => msruns_df_new,
                    "proteins.txt.gz" => proteins_df,
                    "ptmn_locprobs.txt.gz" => ptmn_locprobs_df,
                    "pms_intensities.txt.gz" => pms_intensities_df,
                    #"protgroups.txt.gz" => protgroups_df,
                    #"peptide_to_protgroup.txt.gz" => peptide2protgroups_df,
                    "peptides.txt.gz" => peptides_df,
                    "pepmodstates.txt.gz" => pepmodstates_df,
                    #"protein_to_protgroup.txt.gz" => protein2protgroup_df,
                    "peptide_to_protein.txt.gz" => peptide2protein_df,
                    "ptm_to_protein.txt.gz" => ptm2protein_df,
                    "ptm_to_gene.txt.gz" => ptm2gene_df,
                    "ptmns.txt.gz" => ptmns_df,
                    "ptmn_to_pepmodstate.txt.gz" => ptmn2pms_df,
                    "ptmn_to_ptmngroup.txt.gz" => ptmn2ptmngroup_df,
                    "ptmngroups.txt.gz" => ptmngroups_df,
                    "ptmngroup2pms.txt.gz" => ptmngroup2pms_df]
    @info "Saving $fname..."
    if endswith(fname, ".gz")
        open(GzipCompressorStream, joinpath(output_path, fname), "w") do io
            CSV.write(io, df, delim='\t')
        end
    else
        CSV.write(joinpath(output_path, fname), df, delim='\t')
    end
end

# read PhosphoSitePlus annotations and sequences
psitep_path = joinpath(party3rd_data_path, "PhosphoSitePlus", "20220808")
psitep_annot_dfs = PhosphoSitePlus.read_annotations(psitep_path, verbose=true)
psitep_annots_df = PhosphoSitePlus.combine_annotations(psitep_annot_dfs)

# read reference sequences of PhosphoSitePlus
psitepseqs_df = Fasta.read_phosphositeplus(joinpath(psitep_path, "Phosphosite_PTM_seq_no_header.fasta"))
psitepseqs_df.is_ref_isoform = .!occursin.(Ref(r"-\d+$"), psitepseqs_df.protein_ac);

@save(joinpath(scratch_path, "phosphositeplus_annotations_20220808.jld2"), psitep_annot_dfs, psitep_annots_df)

# map PTMs to PhosphoSitePlus sequences
ptm2psitep_df = PTMExtractor.map_aapos(ptm2protein_df, proteins_df, psitepseqs_df,
                                       destmap_prefix=:psitep_, obj_prefix=:ptm_, objid_col=[:ptm_type, :ptm_pos, :ptm_AA_seq], verbose=true)
filter!(r -> !ismissing(r.destseq_ix), ptm2psitep_df)
ptm2psitep_df.pos_match = coalesce.(ptm2psitep_df.ptm_pos .== ptm2psitep_df.psitep_ptm_pos, false)
ptm2psitep_df.AA_match = coalesce.(ptm2psitep_df.ptm_AA_seq .== uppercase.(coalesce.(ptm2psitep_df.psitep_ptm_AA, '?')), false)
ptm2psitep_df.protein_ac_isoform = [(isnothing(m) ? 1 : parse(Int, m[1])) for m in match.(Ref(r"-(\d+)(?:#.+)?$"), ptm2psitep_df.protein_ac)]
sort!(ptm2psitep_df, [:protein_ac_isoform, :ptm_pos, order(:pos_match, rev=true), order(:AA_match, rev=true)])

ptm_annots_df = innerjoin(select!(filter(r -> r.pos_match && r.AA_match, ptm2psitep_df),
                                  [:ptm_type, :ptm_pos, :ptm_AA_seq, :protein_ac]),
                          psitep_annots_df, on=[:protein_ac, :ptm_pos])
filter(r -> !ismissing(r.kinase_gene_names), ptm_annots_df)
CSV.write(joinpath(output_path, "ptm_annots_$(proj_info.data_ver).txt"), ptm_annots_df, delim='\t')