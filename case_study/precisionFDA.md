# GATHER HG001 and HG002 fastq files

# GRCH37 download
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai

# install bwa
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
tar bwa-0.7.17.tar.bz2
cd bwa-0.7.17/

sudo apt-get install make gcc g++ zlib1g-dev
git clone https://github.com/lh3/bwa.git
cd bwa; make
sudo ln bwa /usr/local/bin/

# align using bwa
REF_FILE=/data/users/common/precisionFDA/GRCH37/human_g1k_v37.fasta
HG001_FA1=/data/users/common/precisionFDA/HG001/fastq/HG001-NA12878-50x_1.fastq
HG001_FA2=/data/users/common/precisionFDA/HG001/fastq/HG001-NA12878-50x_2.fastq
echo ${REF_FILE} ${HG001_FA1} ${HG001_FA2}

cd ${REF_FILE}
bwa index human_g1k_v37.fasta

cd /data/users/common/precisionFDA/HG001/bam/
time bwa mem -M -t 64 ${REF_FILE} ${HG001_FA1} ${HG001_FA2} | samtools sort -@64 -o HG001_GRCh37_pfda.bam -
real    512m41.277s
user    15204m49.407s
sys     53m29.651s

samtools index -@ 40 -b HG001_GRCh37_pfda.bam

### DeepVariant evaluation
### FROM HERE ### https://github.com/google/deepvariant/blob/master/docs/deepvariant-quick-start.md
# DOWNLOAD THE MODEL
BIN_VERSION="0.7.2"
MODEL_VERSION="0.7.2"


cd /data/users/common/dv_model/
MODEL_NAME="DeepVariant-inception_v3-${MODEL_VERSION}+data-wgs_standard"
MODEL_HTTP_DIR="https://storage.googleapis.com/deepvariant/models/DeepVariant/${MODEL_VERSION}/${MODEL_NAME}"
mkdir -p ${MODEL_NAME}
wget -P ${MODEL_NAME} ${MODEL_HTTP_DIR}/model.ckpt.data-00000-of-00001
wget -P ${MODEL_NAME} ${MODEL_HTTP_DIR}/model.ckpt.index
wget -P ${MODEL_NAME} ${MODEL_HTTP_DIR}/model.ckpt.meta


# RUN DEEPVARIANT
OUTPUT_DIR=/data/users/kishwar/deepvariant_outputs/examples/pfda_HG001_grch37
VCF_OUTPUT_DIR=/data/users/kishwar/deepvariant_outputs/vcf_outputs
FINAL_OUTPUT_VCF="${VCF_OUTPUT_DIR}/pfda_HG001_GRCh37.vcf.gz"

REF=/data/users/common/ref/GRCh37/human_g1k_v37.fasta
BAM=/data/users/common/bam/precision_fda/pfda_HG001_GRCh37.bam
MODEL=/data/users/common/dv_model/DeepVariant-inception_v3-0.7.2+data-wgs_standard/model.ckpt

# create images
LOGDIR=/data/users/kishwar/deepvariant_outputs/logs

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOGDIR}"
N_SHARDS=32

sudo docker pull gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}"

mkdir -p "${LOGDIR}"
time seq 0 $((N_SHARDS-1)) | \
  parallel --eta --halt 2 --joblog "${LOGDIR}/log" --res "${LOGDIR}" \
  sudo docker run \
    -v /data/:/data/ \
    -v /home/${USER}:/home/${USER} \
    gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
    /opt/deepvariant/bin/make_examples \
    --mode calling \
    --ref "${REF}" \
    --reads "${BAM}" \
    --sample_name "NA12878" \
    --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" \
    --task {}

CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/call_variants_output.tfrecord.gz"


sudo docker run \
  -v /home/${USER}:/home/${USER} \
  -v /data/:/data/ \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/call_variants \
 --outfile "${CALL_VARIANTS_OUTPUT}" \
 --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" \
 --checkpoint "${MODEL}"



sudo docker run \
   -v /data/:/data/ \
   -v /home/${USER}:/home/${USER} \
   gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
   /opt/deepvariant/bin/postprocess_variants \
   --ref "${REF}" \
   --infile "${CALL_VARIANTS_OUTPUT}" \
   --outfile "${FINAL_OUTPUT_VCF}"

### TO HERE ###


# GENERATE FRIDAY IMAGES "branch: friday_v2"
time seq 0 31 | parallel --ungroup python3 generate_images.py \
--bam /data/users/common/precisionFDA/HG001/bam/HG001_GRCh37_pfda.bam \
--fasta /data/users/common/precisionFDA/GRCH37/human_g1k_v37.fasta \
--threads 32 \
--chromosome_name 1-22 \
--output_dir /data/users/kishwar/friday/image_output/pfda_HG001/ \
--thread_id {}
real    278m28.937s
user    7639m45.944s
sys     29m1.152s