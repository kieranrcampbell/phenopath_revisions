

Ns = [200, 500]
Gs = [500]
prop_interactions = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
reps = 40
algorithms_no_pp = ["dpt", "monocle2", "tscan"]
noises = ["low", "high"]

sig_str = expand("N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}",
                N = Ns, G = Gs, p = prop_interactions, rep = list(range(1, reps + 1)),
                noise = noises)

pseudotime_str = expand("N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}",
                N = Ns, G = Gs, p = prop_interactions, 
                rep = list(range(1, reps + 1)),
                noise = noises,
                alg = algorithms_no_pp)

phenopath_str = expand("N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_phenopath",
                N = Ns, G = Gs, p = prop_interactions, 
                rep = list(range(1, reps + 1)),
                noise = noises)



sim_data_dir = "data/simulations/"
scesets = [sim_data_dir + "scesets/sceset_" + s + ".rds" for s in sig_str]

pseudotimes_no_pp = [sim_data_dir + "pseudotimes/pseudofit_" + s + ".csv" for s in pseudotime_str]
pseudotimes_phenopath = [sim_data_dir + "pseudotimes/phenopathfit_" + s + ".csv" for s in phenopath_str]

dex_qvals = [sim_data_dir + "qvals/qvals_" + s + ".csv" for s in pseudotime_str]
dex_qvals_deseq2 = [sim_data_dir + "deseq2_qvals/qvals_" + s + ".csv" for s in pseudotime_str]
dex_qvals_mast = [sim_data_dir + "mast_qvals/qvals_" + s + ".csv" for s in pseudotime_str]
dex_qvals_monocle = [sim_data_dir + "monocle_qvals/qvals_" + s + ".csv" for s in pseudotime_str]

phenopath_fdata = [sim_data_dir + "phenopath_fdata/fdata_" + s + ".csv" for s in phenopath_str]


datasets = ["chu", "dulken", "hsc", "trapnell"]

lin_scesets = expand("data/scesets/{dataset}-sce.rds", dataset = datasets)
linear_psts = expand("data/linpst/{dataset}_pseudotimes.csv", dataset = datasets)
linear_coefs = expand("data/lincoef/{dataset}.csv", dataset = datasets)

hvg_datasets = ["coad", "brca", "shalek"]
hvgs = [100, 200, 500, 1000, 2000, 4000]
hvg_algorithms = ["phenopath", "monocle"]

paper_scesets = expand("data/paper-scesets/sce_{hvg_dset}_clvm.rds", hvg_dset = hvg_datasets)
hvg_pseudotimes = expand("data/hvg/pseudotime_{hvg_dset}_{hvg}_{hvg_algorithm}.csv",
                        hvg_dset = hvg_datasets, hvg = hvgs, hvg_algorithm = hvg_algorithms)

hvgs_shalek = ["500", "1000", "2000", "4000", "6000", "all"]
shalek_algs = ["dpt", "monocle2", "tscan", "phenopath"]

shalek_pseudotimes = expand("data/shalek_cor/pseudotime_{hvg_shalek}_{hvg_shalek_algorithm}.csv",
                            hvg_shalek = hvgs_shalek, 
                            hvg_shalek_algorithm = shalek_algs)


rule all:
    input:
        # phenopath_fdata,
        # "data/simulations/all_pseudotime_correlations.csv",
        # dex_qvals,
        # "data/simulations/roc.csv",
        # "data/simulations/roc_deseq.csv",
        # "data/simulations/roc_phenopath.csv",
        # "data/simulations/roc_mast.csv",
        # "data/simulations/roc_monocle.csv",
        # linear_psts,
        # linear_coefs,
        # "figs/mast.png"
        # hvg_pseudotimes,
        # "figs/hvg.png",
        shalek_pseudotimes


# Shalek correlation stuff -----------

rule fit_shalek_pseudotimes:
    input:
        "data/paper-scesets/sce_shalek_clvm.rds"
    output:
        "data/shalek_cor/pseudotime_{hvg_shalek}_{hvg_shalek_algorithm}.csv"
    shell:
        "Rscript analysis/shalek_cor/fit_shalek_pseudotime.R --input_sceset {input} --algorithm {wildcards.hvg_shalek_algorithm} --hvg {wildcards.hvg_shalek} --output_csv {output}"


# HVG stuff ------------------

rule fit_hvg_pseudotimes:
    input:
        "data/paper-scesets/sce_{hvg_dset}_clvm.rds"
    output:
        "data/hvg/pseudotime_{hvg_dset}_{hvg}_{hvg_algorithm}.csv"
    shell:
        "Rscript analysis/hvg/fit_hvg_pseudotime.R --input_sceset {input} --algorithm {wildcards.hvg_algorithm} --hvg {wildcards.hvg} --dataset {wildcards.hvg_dset} --output_csv {output}"

rule hvg_figure:
    input:
        hvg_pseudotimes
    output:
        "figs/hvg.png"
    shell:
        "Rscript analysis/hvg/collate_pseudotimes.R"

# Simulations ----------------

rule make_scesets:
    output:
        scesets
    shell:
        "Rscript analysis/simulations/simulate.R"


rule fit_pseudotimes_no_pp:
    input:
        "data/simulations/scesets/sceset_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}.rds"
    output:
        "data/simulations/pseudotimes/pseudofit_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv"
    shell:
        "Rscript analysis/simulations/pseudotime_inference.R --algorithm {wildcards.alg} --input_file {input} --output_file {output}"

rule fit_pseudotimes_phenopath:
    input:
        "data/simulations/scesets/sceset_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}.rds"
    output:
        pseudotime="data/simulations/pseudotimes/phenopathfit_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_phenopath.csv",
        fdata="data/simulations/phenopath_fdata/fdata_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_phenopath.csv"
    shell:
        "Rscript analysis/simulations/pseudotime_inference.R --algorithm phenopath --input_file {input} --output_file {output.pseudotime} --phenopath_fdata_file {output.fdata}"

rule parse_pseudotime_results:
    input:
        scesets,
        pseudotimes_no_pp,
        pseudotimes_phenopath
    output:
        "data/simulations/all_pseudotime_correlations.csv"
    shell:
        "Rscript analysis/simulations/compare_pseudotimes.R --output_file {output}"

rule differential_expression:
    input:
        pseudotime="data/simulations/pseudotimes/pseudofit_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv",
        sceset="data/simulations/scesets/sceset_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}.rds"
    output:
        "data/simulations/qvals/qvals_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv"
    shell:
        "Rscript analysis/simulations/differential_expression.R --input_sceset {input.sceset} --pseudotime_file {input.pseudotime} --output_file {output}"

rule differential_expression_deseq:
    input:
        pseudotime="data/simulations/pseudotimes/pseudofit_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv",
        sceset="data/simulations/scesets/sceset_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}.rds"
    output:
        "data/simulations/deseq2_qvals/qvals_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv"
    shell:
        "Rscript analysis/simulations/differential_expression_deseq2.R --input_sceset {input.sceset} --pseudotime_file {input.pseudotime} --output_file {output}"

rule differential_expression_mast:
    input:
        pseudotime="data/simulations/pseudotimes/pseudofit_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv",
        sceset="data/simulations/scesets/sceset_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}.rds"
    output:
        "data/simulations/mast_qvals/qvals_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv"
    shell:
        "Rscript analysis/simulations/differential_expression_mast.R --input_sceset {input.sceset} --pseudotime_file {input.pseudotime} --output_file {output}"

rule differential_expression_monocle:
    input:
        pseudotime="data/simulations/pseudotimes/pseudofit_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv",
        sceset="data/simulations/scesets/sceset_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}.rds"
    output:
        "data/simulations/monocle_qvals/qvals_N_{N}_G_{G}_p_{p}_rep_{rep}_noise_{noise}_alg_{alg}.csv"
    shell:
        "Rscript analysis/simulations/differential_expression_monocle.R --input_sceset {input.sceset} --pseudotime_file {input.pseudotime} --output_file {output}"



rule roc:
    input:
        scesets,
        dex_qvals
    output:
        "data/simulations/roc.csv"
    shell:
        "Rscript analysis/simulations/calculate_auc.R --qval_dir qvals --output_file {output}"

rule roc_deseq2:
    input:
        scesets,
        dex_qvals_deseq2
    output:
        "data/simulations/roc_deseq.csv"
    shell:
        "Rscript analysis/simulations/calculate_auc.R --qval_dir deseq2_qvals --output_file {output}"

rule roc_mast:
    input:
        scesets,
        dex_qvals_mast
    output:
        "data/simulations/roc_mast.csv"
    shell:
        "Rscript analysis/simulations/calculate_auc.R --qval_dir mast_qvals --output_file {output}"

rule roc_monocle:
    input:
        scesets,
        dex_qvals_monocle
    output:
        "data/simulations/roc_monocle.csv"
    shell:
        "Rscript analysis/simulations/calculate_auc.R --qval_dir monocle_qvals --output_file {output}"


rule roc_phenopath:
    input:
        scesets,
        phenopath_fdata
    output:
        "data/simulations/roc_phenopath.csv"
    shell:
        "Rscript analysis/simulations/calculate_auc_phenopath.R --output_file {output}"

rule sim_example_figure:
    output:
        "figs/simulation_example.rds"
    shell:
        "Rscript analysis/simulations/graph_sims.R"

rule demonstrate_mast:
    input:
        "data/simulations/scesets/sceset_N_200_G_500_p_0.5_rep_8.rds"
    output:
        "figs/mast.png"
    shell:
        "Rscript analysis/simulations/demonstrate_mast.R"

# Testing for linearity of pseudotime


rule fit_linear:
    input:
        "data/scesets/{dataset}-sce.rds"
    output:
        pst="data/linpst/{dataset}_pseudotimes.csv",
        coefs="data/lincoef/{dataset}.csv"
    shell:
        "Rscript analysis/linear/pseudotime_inference.R --input_file {input} --output_file {output.pst} --output_linear_model {output.coefs} --dataset {wildcards.dataset}"



