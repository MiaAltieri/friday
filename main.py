from modules.python.LocalRealignment import LocalAssembler

bam_file = '/data/users/common/giab_flowcells/bam/fl1_chr19.bam'
fasta_file = '/data/users/common/giab_flowcells/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa'
contig = "chr19"

# find active regions
ar = LocalAssembler(bam_file, fasta_file, contig, 8926678, 8927210)
all_reads = ar.perform_local_assembly()
print("Realignment done")