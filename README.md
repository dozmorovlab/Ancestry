# Continental Ancestry Pipeine

Main script directory: `/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev`
Figures: https://docs.google.com/presentation/d/1mwSDvPSpa1S0DQtcua2zZLOZdVX32bXQS2J6XfdCpn4/edit?slide=id.g34552a1d370_0_23#slide=id.g34552a1d370_0_23
Repository: https://github.com/VCUWrightCenter/U4HELPP_PDXProgram

## Notes
* `00_*` are helper scripts, used when needed
* BAM files must have @RG tags to work properly with GATK
    * See https://gatk.broadinstitute.org/hc/en*us/articles/360035890671*Read*groups
    * `02b_Add_RG.sh`
* RNA-Seq BAM files start with `01_MarkDuplicates.sh`
* WGS BAM files were started from `04_HaplotypeCaller_WGS.sh`

## Scripts

### `00_Ancestry_install.sh`
* **Summary:** This script sets up the conda environment, "Ancestry" used for almost all of the following (AEon and gnomAD have separate environments).
* **Inputs:** None.
* **Outputs:** A conda environment.

### `00_combine_WGS.sh`
* **Summary:** This script consolidates whole-genome sequencing (WGS) data by creating symbolic links from multiple source directories to a single target directory.
* **Inputs:** Files located in `$src_dirs`.
* **Outputs:** Symbolic links to the input files, organized in the `/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS` directory.

### `00_Process_BAMs.sh`
* **Summary:** This script pre-processes BAM files. It sorts and indexes the files, filters them to keep only standard chromosomes, and renames the chromosomes in the header based on a mapping file (`human_rename.txt`).
* **Inputs:** BAM files from the directory specified in `00_Settings.sh` and the `human_rename.txt` file.
* **Outputs:** Processed, sorted, and re-headered BAM files in a new directory named `Processed_BAMs`.

### `00_Settings.sh`
* **Summary:** This is a central configuration file. It doesn't perform any data processing itself but instead defines and exports environment variables that are used by the other processing scripts. These variables include paths to data and reference files, computational resource allocations (CPUs and memory), and parameters for various tools. It also creates the directory structure for the analysis results.
* **Inputs:** None.
* **Outputs:** A set of environment variables and an output directory structure.

### `01_MarkDuplicates.sh`
* **Summary:** This script is the first step in the GATK pipeline. It uses the Picard tool `MarkDuplicates` to identify and flag duplicate reads in BAM files, which can arise from PCR artifacts during sequencing library preparation.
* **Inputs:** BAM files that have been sorted by coordinate.
* **Outputs:** New BAM files with duplicate reads marked, along with a metrics file detailing the number of duplicates found.

### `02_SplitNCigarReads.sh`
* **Summary:** This script is specific to RNA-seq data. It uses the GATK tool `SplitNCigarReads` to split reads that span introns (represented as 'N' in the CIGAR string). This makes the data more suitable for variant calling because it prevents variant callers from misinterpreting the large gaps as deletions. The script also creates a sequence dictionary from the reference genome if one doesn't already exist.
* **Inputs:** BAM files from the previous `MarkDuplicates` step and a reference genome in FASTA format.
* **Outputs:** BAM files with reads split across introns, ready for base quality score recalibration.

### `02b_Add_RG.sh`
* **Summary:** This script uses the Picard tool `AddOrReplaceReadGroups` to add or modify read group information in the BAM files. Read groups are essential for downstream GATK tools as they allow for tracking of sample-specific information.
* **Inputs:** BAM files from the `SplitNCigarReads` step.
* **Outputs:** BAM files with updated read group information.

### `03_BQSR.sh`
* **Summary:** This script performs Base Quality Score Recalibration (BQSR) using GATK's `BaseRecalibrator` and `ApplyBQSR` tools. This process corrects for systematic errors in the base quality scores assigned by the sequencing machine, leading to more accurate variant calls. (The script also generates diagnostic plots with `AnalyzeCovariates`).
* **Inputs:** BAM files from the previous step and a set of known variant sites (e.g., from dbSNP) to help differentiate true variation from sequencing errors.
* **Outputs:** Recalibrated BAM files with more accurate base quality scores, a recalibration table, and PDF plots for analysis.

### `04_HaplotypeCaller.sh`
* **Summary:** This is the variant calling step. The script uses GATK's `HaplotypeCaller` to identify SNPs (Single Nucleotide Polymorphisms) and indels (insertions and deletions) in the recalibrated BAM files.
* **Inputs:** The recalibrated BAM files from the BQSR step and the reference genome.
* **Outputs:** VCF (Variant Call Format) files containing the identified genetic variants for each sample.

### `04_HaplotypeCaller_gVCF.sh`
* **Summary:** This script uses GATK's `HaplotypeCaller` to call germline SNPs and indels on BAM files, but it operates in a special mode (`--emit-ref-confidence GVCF`) to produce a genomic VCF (gVCF) file. This gVCF format contains information about all sites in the genome, not just the variant ones, which is necessary for joint genotyping of multiple samples later on.
* **Inputs:** Recalibrated BAM files from the BQSR step (`03_BQSR.sh`).
* **Outputs:** A gVCF file for each input BAM file, stored in the `HaplotypeCaller_gvcf_data_dir`.

### `04_HaplotypeCaller_WGS.sh`
* **Summary:** This script is tailored for running `HaplotypeCaller` on Whole Genome Sequencing (WGS) data (assuming BQSR was already done). While functionally similar to the other `HaplotypeCaller` scripts, it points to specific input and output directories.
* **Inputs:** BAM files from a WGS processing pipeline, located in `/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS`.
* **Outputs:** Standard VCF files containing variant calls for each WGS sample, stored in `/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/WGS_VCF`.

### `04a_process_gVCFs.sh`
* **Summary:** This script handles the joint genotyping process. It first uses GATK's `CombineGVCFs` to merge all the individual gVCF files generated by `04_HaplotypeCaller_gVCF.sh` into a single, multi-sample gVCF. Then, it runs `GenotypeGVCFs` on the combined file to produce a final VCF with genotypes for all samples at all variant sites.
* **Inputs:** Multiple single-sample gVCF files.
* **Outputs:** A combined gVCF and a final, joint-called VCF file containing genotypes for all samples.

### `04b_filter_VCFs.sh`
* **Summary:** This script is for quality control. It applies a series of stringent filters to the VCF files using `bcftools` to remove low-quality or likely false-positive variant calls. It defines multiple filter criteria (e.g., based on quality score, read depth, strand bias) and applies them sequentially, logging the number of SNPs remaining after each step.
* **Inputs:** VCF files from the HaplotypeCaller or GenotypeGVCFs steps.
* **Outputs:** Filtered VCF files (`.vcf.gz`) with low-quality variants removed, along with a TSV file logging the effect of each filter.

### `04c_index_VCFs.sh`
* **Summary:** A simple utility script that uses `bcftools index` to create an index file (`.tbi`) for each compressed VCF file (`.vcf.gz`) in a specified directory. This index is crucial as it allows for rapid random access to the VCF, which is required by most downstream analysis tools.
* **Inputs:** Compressed VCF files.
* **Outputs:** An index file (`.vcf.gz.tbi`) for each input VCF.

### `05_Admixture.sh`
* **Summary:** This script performs an ancestry analysis using the ADMIXTURE program. It first uses `plink` and `bcftools` to prepare the data by merging the user's VCF files with a reference panel (like the 1000 Genomes Project data), filtering for common SNPs, and converting the format. It then runs ADMIXTURE for different numbers of ancestral populations (K) and uses an R script to plot the results.
* **Inputs:** Filtered VCF files and a reference VCF from a population genetics resource.
* **Outputs:** PLINK format files (`.bed`, `.bim`, `.fam`), ADMIXTURE results (`.Q` and `.P` files), and plots visualizing the ancestry components.

### `05_Aeon.sh`
* **Summary:** This script runs Aeon, a tool for predicting genetic ancestry from VCF files. It iterates through the input VCFs, running the `aeon` tool on each to generate ancestry predictions.
* **Inputs:** Filtered VCF files.
* **Outputs:** Ancestry prediction results from the Aeon tool for each sample.

### `05_EthSeq.sh`
* **Summary:** This script runs `EthSEQ` through a custom R script wrapper to generate ancestry predictions.
* **Inputs:** Filtered VCF files.
* **Outputs:** Reports detailing the predicted continental ancestries for each sample.

### `05_gnomAD.sh`
* **Summary:** This script runs a pipeline described by `gnomAD` using their databases through a custom python script wrapper to generate ancestry predictions.
* **Inputs:** Filtered VCF files.
* **Outputs:** Reports detailing the predicted continental ancestries for each sample.

### `05_RAIDS.sh`
* **Summary:** This script runs `RAIDS` through a custom R script wrapper to generate ancestry predictions.
* **Inputs:** Filtered VCF files.
* **Outputs:** Reports detailing the predicted continental ancestries for each sample.

### `05a_JAX.sh`
* **Summary:** This script runs the JAX ancestry pipeline using `nextflow` and generates ancestry predictions. (JAX uses SNPWeights internally here)
* **Inputs:** BAM files, a sample sheet in CSV format, and various reference files specific to the JAX pipeline.
* **Outputs:** A directory of results generated by the Nextflow pipeline, likely including ancestry predictions and quality control metrics.

### `05a2_JAX_standalone.sh`
* **Summary:** A work in progress - similar to `05a_JAX.sh`, but starting from VCF files instead of BAM files
* **Inputs:** VCF files (uncompressed), a sample sheet in CSV format, and various reference files specific to the JAX pipeline.
* **Outputs:** Similar to `05a_JAX.sh`, producing ancestry analysis results.

### `05b_compress_results.sh`
* **Summary:** This script compresses and indexes VCF files that are uncompressed. This is an important step for downstream analysis as `05a_JAX.sh` outputs uncompressed VCF files.
* **Inputs (predicted):** Uncompressed VCF or other text-based result files.
* **Outputs (predicted):** Compressed versions of the input files.

### `06a1_make_pred_table.R`
* **Summary:** This R script aggregates the ancestry predictions from the different tools (Aeon, EthSEQ, gnomAD, JAX SNPWeights, RAIDS, and Admixture) for both RNA-seq and WGS data. It reads the output files from each tool, cleans up sample names, and combines them into two master summary tables (one for RNA, one for WGS). Finally, it saves these summary tables as both CSV and Excel files.
* **Inputs:** Various output files from the ancestry prediction tools (e.g., CSVs, TSVs).
* **Outputs:** Two CSV files (`RNA_results.csv`, `WGS_results.csv`) and two Excel files (`RNA_results.xlsx`, `WGS_results.xlsx`) containing the collated ancestry predictions.

### `06a2_Count_reads.sh`
* **Summary:** This script uses `samtools` to count the total number of reads in each BAM file from either the WGS or the RNA-seq pipeline.
* **Inputs:** A directory of BAM files.
* **Outputs:** A TSV file with two columns: the BAM file name and the total read count.

### `06a3_GATK_vs_JAX.sh`
* **Summary:** This script compares the variant calls generated by two different pipelines: the GATK pipeline and the JAX pipeline. For each sample, it uses `bcftools isec` to find the variants that are unique to each pipeline and those that are shared between them. It then counts the number of SNPs and indels in each of these categories and compiles the results into a single summary CSV file.
* **Inputs:** Two directories of VCF files, one from the GATK pipeline and one from the JAX pipeline.
* **Outputs:** A CSV file summarizing the variant overlap between the two pipelines for each sample.

### `06a4_overlap_aims.sh`
* **Summary:** This script overlaps samples (VCF files) with each internal model of each tool 
* **Inputs:** VCF files for each sample and one or more BED files containing the models for tools used.
* **Outputs:** TSV files for each sample listing the SNPs that overlap.

### `06a4a_papers_to_bed.R`
* **Summary:** THis script converts the AIMs found in papers to GRCh38, and creates BED files representing those SNPs / AIMs.
* **Inputs:** Text or Excel files downloaded from scientific publications.
* **Outputs:** BED files, where each file represents a set of AIMs from a specific paper.

### `06b0_sources.R`
* **Summary:** This is a central R script for all the libraries, directories, and shared functions for downstream scripts generating figures and tables for the final analysis. It loads many R libraries and defines file paths to the summary tables created by `06a1_make_pred_table.R`. It also contains helper functions and aesthetic definitions (like official tool names and color palettes) that will be used.
* **Inputs:** Data sources, color choices, etc.
* **Outputs:** None.

### `06b1_f1_acc.R`
* **Summary:** This R script is for evaluating the accuracy of the different ancestry prediction tools. It reads the prediction results for both RNA-seq and WGS data, as well as SIRE information per sample. It then calculates the F1 score (a measure of accuracy) for each tool by comparing their predictions to the ground truth. Finally, it generates bar plots to visualize these F1 scores, allowing for a direct comparison of the tools' performance.
* **Inputs:** `06b0_sources.R`
* **Outputs:** Bar plots / tables showing confusion matricies and f1 / accracy scores for each ancestry prediction method on both RNA-seq and WGS data.

### `06b2_GATK_vs_JAX.R`
* **Summary:** This R script is dedicated to visualizing the comparison between the GATK and JAX variant calling pipelines. It reads the summary CSV file generated by `06a3_GATK_vs_JAX.sh` and creates a series of boxplots and bar plots. These plots illustrate the number of SNPs that are unique to each pipeline versus those that are shared, providing a clear visual representation of the concordance between the two methods for both RNA-seq and WGS data.
* **Inputs:** `06b0_sources.R`
* **Outputs:** A variety of plots (boxplots, bar plots) / tables comparing the SNP counts between the GATK and JAX pipelines.

### `06b3_visualize_overlap.R`
* **Summary:** This R script visualizes the number of Ancestry Informative Markers (AIMs) / SNPs that were successfully identified in each sample. It reads the BED files produced by `06a4_overlap_aims.sh` and generates boxplots showing the distribution of AIMs found per sample for different sets of markers (e.g., those defined by specific tools or papers). The script also calculates and plots the percentage of total AIMs / SNPs that were recovered before / after filtering.
* **Inputs:** `06b0_sources.R`
* **Outputs:** Boxplots / tables that show both the raw counts and percentages of AIMs detected in the samples.

### `06b4_visualize_snps_reads.R`
* **Summary:** This script generates plots / tables to explore the relationship between the number of sequencing reads and the number of SNPs called for each sample.
* **Inputs (predicted):** `06b0_sources.R`.
* **Outputs (predicted):** Plots showing the relationship between read counts and SNP counts for RNA-seq and WGS data.

### `06b5_visualize_models.R`
* **Summary:** This script analyzes and visualizes the overlap between the different sets of Ancestry Informative Markers (AIMs) used by the various tools and papers. It reads the BED files defining each set of markers, counts the number of markers in each, and generates plots to show how many markers are unique to each set versus how many are shared between them.
* **Inputs:** `06b0_sources.R`
* **Outputs:** Plots that visualize the size of and overlap between the different AIMs / SNPs panels.

### `06b6_resources_plots.R`
* **Summary:** This script is for performance analysis. It reads a file containing computational resource usage data (CPU time, memory) for each of the main ancestry analysis tools. It then creates bar plots to compare the resource consumption of the different tools, separated by data type (RNA-seq vs. WGS). This is useful for understanding the computational cost of each method.
* **Inputs:** `06b0_sources.R`
* **Outputs:** Plots comparing the maximum memory usage and total CPU time for each ancestry tool.

### `06b7_get_pca.R`
* **Summary:** This script visualizes Principal Component Analysis (PCA) on the sample variants, then create plots.
* **Inputs:** `06b0_sources.R`
* **Outputs:** PCA and other plots.

### `06b8_combine_figs.R`
* **Summary:** This script acts as a master assembler for creating the final, multi-panel figures for a publication. It doesn't generate new plots itself, but instead loads the plots created by the other `06b*` scripts (e.g., `06b1_f1_acc.R`, `06b2_GATK_vs_JAX.R`) and arranges them into composite figures using the `patchwork` R library. It defines the layout and adds labels (A, B, C, etc.) to each panel.
* **Inputs:** `06b*` scripts.
* **Outputs:** The final, publication-ready figures saved in formats like PNG and SVG.

### `06c1_merge_samples.sh`
* **Summary:** This script combines VCF files from all individual samples into a single, multi-sample VCF file.
* **Inputs:** A directory of single-sample VCF files.
* **Outputs:** A single merged VCF file containing the variant information for all samples.

### `06c2_overlap_1k.sh`
* **Summary:** This script is be used to find samples variants (SNPs) that are also present in the 1000 Genomes Project reference geneomes.
* **Inputs:** The merged VCF file of all user samples and a merged VCF (all chromosomes) from the 1000 Genomes Project.
* **Outputs:** A VCF file per sample containing only the overlapping variant sites.

### `06d_randomize_and_separate.sh`
* **Summary:** This script is a key part of a simulation study. It takes a VCF file that contains variants from many individuals from the 1000 Genomes Project present in the merged VCF of samples and randomly samples a specified number of SNPs (e.g., 10, 100, 1,000, etc.). It repeats this process multiple times with different random seeds to create many small VCF files, each containing a random subset of SNPs.
* **Inputs:** A VCF file containing variants that overlap between the user's samples and the 1000 Genomes Project.
* **Outputs:** A directory structure containing many small VCF files, organized by the number of SNPs and the random seed used.

### `06e_predict.sh`
* **Summary:** This script automates the process of running multiple ancestry prediction tools on the small, simulated VCF files created by `06d_randomize_and_separate.sh`. It iterates through the directories of simulated data and generates and submits SLURM batch jobs to run the ancestry prediction pipelines on each set of files.
* **Inputs:** The directories of small VCF files generated by the previous script.
* **Outputs:** SLURM job scripts for running the ancestry analyses, and ultimately, the prediction results for each simulated dataset.

### `06f_count_1k_snps.sh`
* **Summary:** This is a quality control and data summary script. It recursively searches through the directories of simulated VCF files and uses `bcftools` to count the number of SNPs in each file. It then compiles these counts into a single TSV file.
* **Inputs:** The directory structure containing the many small, simulated VCF files.
* **Outputs:** A single TSV file with the path to each VCF and the number of SNPs it contains.

### `06g_aggregate_1k_results.R`
* **Summary:** This R script is the final analysis step for the simulation study (1k for 1000 genomes project). It reads all the ancestry prediction results generated by the `06e_predict*.sh` scripts, as well as the true ancestry information for the 1000 Genomes samples. It then calculates the accuracy (specifically, the F1 score) of each tool for different numbers of input SNPs. Finally, it generates plots to visualize how the accuracy of each tool changes as the number of available SNPs increases or decreases.
* **Inputs:**
    * The directories containing the ancestry prediction results from the simulation runs.
    * The SNP count summary file from `06f_count_1k_snps.sh`.
    * The metadata for the 1000 Genomes Project samples.
* **Outputs:**
    * Plots showing the F1 scores for each tool at different SNP counts.
    * A CSV file that aggregates all the prediction results into a single table.

### `07a_Admixture_aims.sh`
* **Summary:** Extracts and creates BED file from internal tool's model.
* **Inputs:** Tool's internal files.
* **Outputs:** A BED file listing the genomic coordinates of the tool's AIMs / SNPs.

### `07a_Aeon_aims.R`
* **Summary:** Extracts and creates BED file from internal tool's model.
* **Inputs:** Tool's internal files.
* **Outputs:** A BED file listing the genomic coordinates of the tool's AIMs / SNPs.

### `07a_EthSEQ_aims.R`
* **Summary:** Extracts and creates BED file from internal tool's model.
* **Inputs:** Tool's internal files.
* **Outputs:** A BED file listing the genomic coordinates of the tool's AIMs / SNPs.

### `07a_gnomAD_aims.py`
* **Summary:** Extracts and creates BED file from internal tool's model.
* **Inputs:** Tool's internal files.
* **Outputs:** A BED file listing the genomic coordinates of the tool's AIMs / SNPs.

### `07a_RAIDS_aims.R`
* **Summary:** Extracts and creates BED file from internal tool's model.
* **Inputs:** Tool's internal files.
* **Outputs:** A BED file listing the genomic coordinates of the tool's AIMs / SNPs.

### `07a_snpweights_aims.R`
* **Summary:** Extracts and creates BED file from internal tool's model.
* **Inputs:** Tool's internal files.
* **Outputs:** A BED file listing the genomic coordinates of the tool's AIMs / SNPs.
