# build the image
docker build -t friday ./Dockerfile/

docker run friday python3 friday/generate_images.py

time sudo docker run \
  -v "/home/${USER}":"/home/${USER}" \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/call_variants \
  --outfile=${OUTPUT_DIR}/HG002.cvo.tfrecord-00000-of-00064.gz \
  --examples=${OUTPUT_DIR}/HG002.examples.tfrecord-00000-of-00064.gz \
  --checkpoint="${MODEL}"


BAM=/data/users/common/giab_flowcells/bam/fl1_chr19.bam
FASTA=/data/users/common/giab_flowcells/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
VCF=/data/users/common/giab_flowcells/vcf/NA12878_GRCh38.vcf.gz
BED=/data/users/common/giab_flowcells/bed/GRCH38_confident.bed
OUTPUT_DIR=/data/users/kishwar/docker_output
CHROMOSOME_NAME=chr19:6840300-6840310
THREADS=1

time seq 0 0 | parallel -k --lb docker run \
               -v /data/users/common:/data/users/common \
               -v /data/users/kishwar:/data/users/kishwar \
               friday python3 friday/generate_images.py \
               --bam ${BAM} \
               --fasta ${FASTA} \
               --vcf ${VCF} \
               --bed ${BED} \
               --output_dir ${OUTPUT_DIR} \
               --train_mode true \
               --chromosome_name ${CHROMOSOME_NAME} \
               --threads ${THREADS} \
               --thread_id {}
