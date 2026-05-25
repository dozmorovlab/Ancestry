# Various settings for later scripts

# Activate conda environment
# 
# Ex:
# conda create --name Ancestry --file ~/../harrell_lab/bulkRNASeq/Ancestry/installed_packages.txt
# conda create --name Ancestry --file ~/../harrell_lab/bulkRNASeq/Ancestry/installed_packages.txt --override-channels -c conda-forge -c defaults
# srun -c 22 --nodelist=apollo17 --mem=494G -J AncestryWGS --pty bash
# ca
# conda activate /lustre/home/wallbp/mambaforge/envs/Ancestry
 
# Input BAM directory
BAM_DIR="/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/RNA"
BAM_DIR=$(realpath $BAM_DIR)

# REGEX for selecting BAMs within BAM_DIR
BAM_REGEX="*.genome.sorted.bam"

# Output data directory
OUT_DIR=./Ancestry_Results_GATK
OUT_DIR=$(realpath $OUT_DIR)

# Resouce directory
RESOURCE_DIR=~/Ancestry

# Computational resources
MAX_CPUS=22
MAX_MEM_GB=494

# Human reference
# https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA,RNA_FASTA,GENOME_GTF
HUMAN_GENOME_REF="./References/GRCh38/genome.fa"
HUMAN_GENOME_GTF=${RESOURCE_DIR}/References/Human/ncbi_dataset/data/GCF_000001405.40/genomic.gtf
HUMAN_RNA_REF=${RESOURCE_DIR}/References/Human/ncbi_dataset/data/GCF_000001405.40/rna.fna
HUMAN_SEQ_REPORT=${RESOURCE_DIR}/References/Human/GCF_000001405.40_sequence_report.jsonl

# Mouse reference
# https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.27/download?include_annotation_type=GENOME_FASTA,RNA_FASTA,GENOME_GTF
#MOUSE_GENOME_REF=${RESOURCE_DIR}/References/Mouse/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna
#MOUSE_GENOME_GTF=${RESOURCE_DIR}/References/Mouse/ncbi_dataset/data/GCF_000001635.27/genomic.gtf
#MOUSE_RNA_REF=${RESOURCE_DIR}/References/Mouse/ncbi_dataset/data/GCF_000001635.27/rna.fna

# Other files
chromosome_map=${RESOURCE_DIR}/References/Human/chromosome_chr.map

# https://labshare.cshl.edu/shares/krasnitzlab/aicsPaper/matGeno1000g.gds
matGeno1000g=${RESOURCE_DIR}/References/Population/matGeno1000g.gds

# https://labshare.cshl.edu/shares/krasnitzlab/aicsPaper/matAnnot1000g.gds
matAnnot1000g=${RESOURCE_DIR}/References/Population/matAnnot1000g.gds

# https://labshare.cshl.edu/shares/krasnitzlab/aicsPaper/snvSel0.01.vcf.gz
snvSel0=${RESOURCE_DIR}/References/Population/snvSel0.01.vcf.gz

# For gnomAD
#https://gnomad.broadinstitute.org/downloads
#gsutil -m cp -r -n gs://gcp-public-data--gnomad/release/4.0/pca ~/gnomad/
loadings_path=~/gnomad/pca/gnomad.v4.0.pca_loadings.ht
rf_model_path=~/gnomad/pca/gnomad.v4.0.RF_fit.onnx

#gsutil -m cp -r -n gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt ~/gnomad/
hgdp_1kg_mt_path=~/gnomad/gnomad.genomes.v3.1.2.hgdp_1kg_subset_dense.mt

# Platform for BAMs
PL=ILLUMINA

# RAIDS settings
# integer representing the number of samples that will be selected for
# each subcontinental population present in the 1KG GDS file
number_of_profiles=50

# https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard
# The maximum offset between two duplicate clusters
# in order to consider them optical duplicates. The
# default is appropriate for unpatterned versions of
# the Illumina platform. For the patterned flowcell
# models, 2500 is moreappropriate. For other
# platforms and models, users should experiment to
# find what works best.
OPTICAL_DUPLICATE_PIXEL_DISTANCE=100

# Known sites for variant calling
# https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz
#known_sites="${RESOURCE_DIR}/References/Human/snps/GCF_000001405.40.gz"
known_sites=(
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr1.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr2.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr3.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr4.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr5.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr6.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr7.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr8.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr9.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr10.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr11.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr12.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr13.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr14.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr15.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr16.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr17.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr18.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr19.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr20.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr21.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr22.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chrX.vcf.bgz"
    #"/lustre/home/wallbp/DNA/gnomad_4.1/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chrY.vcf.bgz"
    "/lustre/home/wallbp/DNA/gnomad.genomes.v4.1.sites.reduced.merged.vcf.gz"
)
# Create Directory Structure
# ---------------------------------------------------

mkdir -p $OUT_DIR

# Data
data_dir=$OUT_DIR/Data
mkdir -p $data_dir

MarkDuplicates_bam_dir="$data_dir/MarkDuplicates"
mkdir -p $MarkDuplicates_bam_dir
SplitNCigarReads_bam_dir="$data_dir/SplitNCigarReads"
mkdir -p $SplitNCigarReads_bam_dir
BaseRecalibrator_data_dir="$data_dir/BaseRecalibrator"
mkdir -p $BaseRecalibrator_data_dir
HaplotypeCaller_data_dir="$data_dir/HaplotypeCaller"
mkdir -p $HaplotypeCaller_data_dir
HaplotypeCaller_gvcf_data_dir="$data_dir/HaplotypeCaller_gVCF"
mkdir -p $HaplotypeCaller_gvcf_data_dir
RAIDS_data_dir="$data_dir/RAIDS"
mkdir -p $RAIDS_data_dir
SNVStory_singularity="$data_dir/singularity/snvstory"
mkdir -p $SNVStory_singularity
#SNVStory_data_dir="$data_dir/SNVStory"
#mkdir -p $SNVStory_data_dir
EthSEQ_data_dir="$data_dir/EthSEQ"
mkdir -p $EthSEQ_data_dir
renamed_vcf_data_dir="$data_dir/renamed_VCFs"
mkdir -p $renamed_vcf_data_dir
gnomAD_data_dir="$data_dir/gnomAD"
mkdir -p $gnomAD_data_dir
GenotypeGVCFs_data_dir="$data_dir/GenotypeGVCFs"
mkdir -p $GenotypeGVCFs_data_dir

# Reports
report_dir=$OUT_DIR/Reports
mkdir -p $report_dir

MarkDuplicates_report="$report_dir/MarkDuplicates"
mkdir -p $MarkDuplicates_report
BaseRecalibrator_report="$report_dir/BaseRecalibrator"
mkdir -p $BaseRecalibrator_report
analysis_report="$report_dir/Analysis"
mkdir -p $analysis_report
