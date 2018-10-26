from modules.python.LocalRealignment import LocalAssembler
from modules.python.CandidateFinder import CandidateFinder
from modules.python.CandidateLabler import CandidateLabeler

bam_file = '/data/users/common/giab_flowcells/bam/fl1_chr19.bam'
vcf_file = '/data/users/common/giab_flowcells/vcf/NA12878_GRCh38.vcf.gz'
fasta_file = '/data/users/common/giab_flowcells/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa'
contig = "chr19"
pos_start = 8051441
pos_end = 8052441
# 8926678, 8927210
# find active regions
ar = LocalAssembler(bam_file, fasta_file, contig, pos_start, pos_end)
all_reads = ar.perform_local_assembly()

candidate_finder = CandidateFinder(fasta_file, contig, pos_start, pos_end)
candidates = candidate_finder.find_candidates(all_reads)
for candidate in candidates:
    candidate.print()
candidate_labler = CandidateLabeler(fasta_file, vcf_file)
labeled_candidates = candidate_labler.get_labeled_candidates(contig, pos_start, pos_end, candidates)

print("ALL DONE")