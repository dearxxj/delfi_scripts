# DELFI (DNA Evaluation of Fragments for early Interception)
Detection of abnormalities in cell-free DNA (cfDNA) by genome-wide analysis of fragmentation patterns based on low-coverage Whole Genome Sequencing of isolated cfDNA.
For more details, please refer to Cristiano et al. 2019 [ref](https://www.nature.com/articles/s41586-019-1272-6)
Forked from [https://github.com/cancer-genomics/delfi_scripts](https://github.com/cancer-genomics/delfi_scripts)

-------------------------------------------------------------------------------------


## Step 1. Obtain genome-wide fragment profile using bin size of 100kb
```
Usage: Rscript --vanilla delfi_WGS.R --bam <bam> --out_prefix <output_prefix> --out_dir <output_directory>
Params:
    --bam            A WGS bam file of 1-2X coverage
    --out_prefix     Output prefix
    --out_dir        Output directory

Output:
    *_bin_100kb.rds
```
This script generate a genome-wide fragment profile at 100kb bins from a bam file. The fragment profile is GC-corrected and outlier-filtered.

Here, I down-sampled a WGS for K562 (Leukemia) cell line to 2X coverage, and generated 100kb-binned fragment profile as cancer sample. Fragment profile for healthy control was obtained from a WGS for NA12878 (lymphocyte) downsampled to 2X coverage.

WGS for K562 was downloaded from [https://www.encodeproject.org/experiments/ENCSR045NDZ/](https://www.encodeproject.org/experiments/ENCSR045NDZ/)
WGS for NA12878 was downloaded from SRA (accession number: ERR3239334)


## Step 2. (Testing only) Permutation of the fragment profile from Step 1 to obtain the same number of fragment profiles used in the paper
Since the original dataset used in the paper is access-controlled, a set of sample-size matched fragment profiles were generated by random sampling from the fragment profile for K562 (Cancer) and NA12878 (Healthy).
```
Rscript --vanilla random_sample_frag.R

Output:
    A set of fragment profiles for 423 samples under ../random_100kb folder
```


## Step3. Aggregate the fragment profile for all samples and predict healthy/cancer individuals with stochastic gradient boosting model.
This script combined all individual fragment profile and merged 100kb bins to 5mb bins, calculated summary statistics, and performed the GBM training and prediction.
```
Usage: Rscript --vanilla summarize_data_gbm.R --bindir <100kb bin directory> --sample <sample_reference.csv> --outdir <output directory>

Params:
    --bindir         Directory for fragment profiles at 100kb bins
    --sample         Sample reference csv file
    --outdir        Output directory

Output:
    bins_100kbcompartments.rds
    bins_5mbcompartments.rds
    summary_tibble.rds
    models_list.rds
    predictions_gbm.csv
```
