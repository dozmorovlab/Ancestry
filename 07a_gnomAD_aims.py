import pandas as pd
import hail as hl

loadings_path = '/lustre/home/wallbp/gnomad/pca/gnomad.v4.0.pca_loadings.ht'

loadings = hl.read_table(loadings_path)

bed_table = loadings.select(
    chrom = loadings.locus.contig,
    start = loadings.locus.position - 1,
    end = loadings.locus.position
)

bed_pd = bed_table.to_pandas()
bed_pd = bed_pd.drop_duplicates()
bed_pd = bed_pd[['chrom', 'start', 'end']]
bed_pd.to_csv(
    "/lustre/home/harrell_lab/bulkRNASeq/Ancestry_dev/aims_tools/models/gnomAD.bed",
    sep="\t",
    index=False,
    header=False
)
