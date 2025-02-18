# lmn_bdnf
RNA-seq and proteomic analyses on i3 LMN

Most of the analyses rely on 'new_ratio_bayesian_p_de' which is table S1 in the manuscript

Scripts - in no particular order

`new_ratio_functions` - contains the helper functions used in `full_merging_gene_estimates`

`full_merging_gene_estimates` - combines the outputs from GRAND-SLAM and estimates the condition specific new_total_rna functions

`comparing_other_bdnf_paper` - contains the comparisons on log2FC in the GABAergic neurons as well as the 
gene set enrichment analyses

`making_sample_betas` - Helper functions to draw beta distributions something like the supplementary figure S2B

`slam_go_comparison_timepoints` - Runs the GO analyses for genes by the differential labelling and differential expression and time

`category_plots_plus_log2_ierg` - Makes the categorization bar chart, plots the distribution of the NTR/total log2fold change by time and control/BNDF, 

`GO_total_comparisons` - runs the GO analyses for genes by just the total RNA differential expression by time 

`run_phos_de_analysis` - run limma to perform differential phosphosite analyses, also pulls back the sequence of the protein from UniProt based on the uniprot ID inorder to get the full region flanking the phosphosite for the analyses on kinase prediction

`proteome_analysis` - performs the GO analyses on the significant phosphosite genes, does plotting of phospho-volcanoes

`splicing_bdnf_effect` - reads in splicing analyses from MAJIQ, annotates as CDS, 3 or 5'utr, makes plots - use table S2 for the input here

`write_fasta_for_motif_hunting_logo` - reads in the results of PhosophoSitePlus kinase prediction, makes plots, also the sequence logo

`plot_tf_activity` - makes the plots on TF activity


Data

`estimate_list_full` - a long format of all the BDNF and control samples outputs from GRAND-SLAM

`enrichment-analysis-result-table_onehour_pvalue` - output of the PhosphoSitePlus kinase prediction tool after giving the surrounding sequence, log2FC, and pvalue of limma as input to the PhosphoSitePlus server - 1 hr

`enrichment-analysis-result-table_onehour_pvalue` - output of the PhosphoSitePlus kinase prediction tool after giving the surrounding sequence, log2FC, and pvalue of limma as input to the PhosphoSitePlus server - 6 hr

`new_ratio_bayesian_p_de` - basically table S1 in the manscript, contains all the differential labelling, expression, etc stuff

These are the Rdata objects containing the results of DESEQ2 on the counts by featureCounts on the aligned BAM files
`hour_one_featurecounts`
`hour_two_featurecounts`
`hour_six_featurecounts`