from modules.python.ActiveRegionFinder import ActiveRegionFinder

bam_file = '/data/users/common/giab_flowcells/bam/fl1_chr19.bam'
fasta_file = '/data/users/common/giab_flowcells/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa'
contig = "chr19"

# find active regions
ar = ActiveRegionFinder(bam_file, fasta_file, contig, 51937245, 51938245)
active_regions = ar.find_active_region()
print(active_regions)
