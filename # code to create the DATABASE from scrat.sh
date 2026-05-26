# code to create the DATABASE from scratch, starting with an empty trs.db and raw data files

# 1. create schema
python 01_create_schema.py --db data/trs.db

# 2. fetch STRipy (already done, skip if you have data/raw/stripy_loci.json)
python 02_fetch_stripy.py --out data/raw/stripy_loci.json

# 3. ingest TRExplorer
python 03_ingest_trexplorer.py \
    --tsv data/raw/TR_catalog.5599658_loci.20260123_034640.tsv \
    --db data/trs.db \
    --stripy data/raw/stripy_loci.json

# 4. export RExPRT input (all 1M loci)
python 04_prepare_rexprt_input.py --db data/trs.db --out data/processed/rexprt_input_all.tsv

# 5. run RExPRT overnight
cd ~/test_TREX/TREX
bash rexptr_standalone.sh /home/rachele/TRpsyDB/scripts/data/processed/rexprt_input_all.tsv


cd ~/db_rex
for chunk in /home/rachele/TRpsyDB/scripts/data/processed/rexprt_chunks/chunk_*.tsv; do
    echo "=== Processing $(basename $chunk) ==="
    bash rexptr_standalone.sh $chunk
done


# 6. once RExPRT done — clean output (back in scripts/)

cd ~/db_rex/results/rexprt

head -1 chunk_000_rex_input_TRsAnnotated_RExPRTscores.txt \
    > ~/TRpsyDB/scripts/data/processed/rexprt_all_merged.txt

for i in $(seq -w 000 021); do
    tail -n +2 chunk_${i}_rex_input_TRsAnnotated_RExPRTscores.txt \
        >> ~/TRpsyDB/scripts/data/processed/rexprt_all_merged.txt
done

wc -l ~/TRpsyDB/scripts/data/processed/rexprt_all_merged.txt

cd ~/TRpsyDB/scripts

python 05_clean_rexprt_output.py \
    --rexprt ~/TRpsyDB/scripts/data/processed/rexprt_all_merged.txt \
    --input  data/processed/rexprt_input_all.tsv \
    --out    data/processed/rexprt_cleaned.tsv


# 7. ingest RExPRT scores
python 06_ingest_rexprt.py \
    --input data/processed/rexprt_cleaned.tsv \
    --db    data/trs.db

    
# 8. download EPD promoters
mkdir -p ~/TRpsyDB/data/raw/regulatory
cd ~/TRpsyDB/data/raw/regulatory
wget https://epd.expasy.org/ftp/epdnew/H_sapiens/current/Hs_EPDnew.bed

# 9. download ENCODE cCREs
wget "https://downloads.wenglab.org/V3/GRCh38-cCREs.bed" -O cCREs_hg38.bed

# 10. download SEdb3 sample info + bulk SE file
# → go to http://www.licpathway.net/sedb/download.php in browser
# → download: Sample information file → save as SEdb_sample_info.txt
# → download: Bulk SE file (hg38, human) → save as SEdb_SE_hg38.txt
# → download: Bulk TE file (hg38, human) → save as SEdb_TE_hg38.txt

# 11. download TRxQTL brain files
mkdir -p ~/TRpsyDB/data/raw/trxqtl
cd ~/TRpsyDB/data/raw/trxqtl
wget https://wlcb.oit.uci.edu/TRxQTL/download/ROSMAP.DLPFC.TR-eQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/ROSMAP.DLPFC.TR-sQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/ROSMAP.DLPFC.TR-3aQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/ROSMAP.DLPFC.TR-hQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/ROSMAP.DLPFC.TR-mQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cortex_Frontal.TR-eQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cortex_Frontal.TR-sQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cortex_Frontal.TR-3aQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cortex_Temporal.TR-eQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cortex_Temporal.TR-sQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cortex_Motor.TR-eQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cortex_Motor.TR-sQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cerebellum.TR-eQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/NYGCALS.Cerebellum.TR-sQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/AnswerALS.iPSC-MN.TR-eQTL.txt.gz
wget https://wlcb.oit.uci.edu/TRxQTL/download/AnswerALS.iPSC-MN.TR-caQTL.txt.gz




# GTEx v8 — all 49 tissues TR-eQTL
BASE="https://wlcb.oit.uci.edu/TRxQTL/download"

wget $BASE/Adipose_Subcutaneous.TR-eQTL.txt.gz
wget $BASE/Adipose_Visceral_Omentum.TR-eQTL.txt.gz
wget $BASE/Adrenal_Gland.TR-eQTL.txt.gz
wget $BASE/Artery_Aorta.TR-eQTL.txt.gz
wget $BASE/Artery_Coronary.TR-eQTL.txt.gz
wget $BASE/Artery_Tibial.TR-eQTL.txt.gz
wget $BASE/Brain_Amygdala.TR-eQTL.txt.gz
wget $BASE/Brain_Anterior_cingulate_cortex_BA24.TR-eQTL.txt.gz
wget $BASE/Brain_Caudate_basal_ganglia.TR-eQTL.txt.gz
wget $BASE/Brain_Cerebellar_Hemisphere.TR-eQTL.txt.gz
wget $BASE/Brain_Cerebellum.TR-eQTL.txt.gz
wget $BASE/Brain_Cortex.TR-eQTL.txt.gz
wget $BASE/Brain_Cortex_Frontal_BA9.TR-eQTL.txt.gz
wget $BASE/Brain_Hippocampus.TR-eQTL.txt.gz
wget $BASE/Brain_Hypothalamus.TR-eQTL.txt.gz
wget $BASE/Brain_Nucleus_accumbens_basal_ganglia.TR-eQTL.txt.gz
wget $BASE/Brain_Putamen_basal_ganglia.TR-eQTL.txt.gz
wget $BASE/Brain_Spinal_cord_cervical_c-1.TR-eQTL.txt.gz
wget $BASE/Breast_Mammary_Tissue.TR-eQTL.txt.gz
wget $BASE/Cells_EBV-transformed_lymphocytes.TR-eQTL.txt.gz
wget $BASE/Cells_Cultured_fibroblasts.TR-eQTL.txt.gz
wget $BASE/Colon_Sigmoid.TR-eQTL.txt.gz
wget $BASE/Colon_Transverse.TR-eQTL.txt.gz
wget $BASE/Esophagus_Gastroesophageal_Junction.TR-eQTL.txt.gz
wget $BASE/Esophagus_Mucosa.TR-eQTL.txt.gz
wget $BASE/Esophagus_Muscularis.TR-eQTL.txt.gz
wget $BASE/Heart_Atrial_Appendage.TR-eQTL.txt.gz
wget $BASE/Heart_Left_Ventricle.TR-eQTL.txt.gz
wget $BASE/Liver.TR-eQTL.txt.gz
wget $BASE/Lung.TR-eQTL.txt.gz
wget $BASE/Muscle_Skeletal.TR-eQTL.txt.gz
wget $BASE/Nerve_Tibial.TR-eQTL.txt.gz
wget $BASE/Ovary.TR-eQTL.txt.gz
wget $BASE/Pancreas.TR-eQTL.txt.gz
wget $BASE/Pituitary.TR-eQTL.txt.gz
wget $BASE/Prostate.TR-eQTL.txt.gz
wget $BASE/Skin_Not_Sun_Exposed.gz
wget $BASE/Skin_Sun_Exposed.TR-eQTL.txt.gz
wget $BASE/Small_Intestine_Terminal_Ileum.TR-eQTL.txt.gz
wget $BASE/Spleen.TR-eQTL.txt.gz
wget $BASE/Stomach.TR-eQTL.txt.gz
wget $BASE/Testis.TR-eQTL.txt.gz
wget $BASE/Thyroid.TR-eQTL.txt.gz
wget $BASE/Uterus.TR-eQTL.txt.gz
wget $BASE/Vagina.TR-eQTL.txt.gz
wget $BASE/Whole_Blood.TR-eQTL.txt.gz

# 12. download RefFlat for TSS/SJ distances
mkdir -p ~/TRpsyDB/data/raw/refflat
cd ~/TRpsyDB/data/raw/refflat
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz

# 13. download APPRIS principal isoforms
wget "https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/homo_sapiens/GRCh38/appris_data.principal.txt" -O appris_data.principal.hg38.txt

# 14. download hg38 intron + exon files (needed for SJ analysis, same as your pipeline)
# these come from UCSC refFlat — generate them with R exactly as your pipeline does:
Rscript - << 'EOF'
library(dplyr)
refflat <- read.delim("refFlat.txt", header=FALSE,
    col.names=c("geneName","name","chrom","strand","txStart","txEnd",
                "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds"))

# generate exon BED
exons <- refflat %>%
    rowwise() %>%
    mutate(
        starts = strsplit(exonStarts, ",")[[1]],
        ends   = strsplit(exonEnds,   ",")[[1]]
    ) %>%
    unnest(c(starts, ends)) %>%
    transmute(chrom, start=as.integer(starts), end=as.integer(ends),
              geneName, name, strand)
write.table(exons, "hg38_exon_refFlat.txt", sep="\t", row.names=FALSE, quote=FALSE)

# generate intron BED
introns <- refflat %>%
    rowwise() %>%
    mutate(
        starts = strsplit(exonStarts, ",")[[1]],
        ends   = strsplit(exonEnds,   ",")[[1]]
    ) %>%
    unnest(c(starts, ends)) %>%
    group_by(name) %>%
    mutate(
        intron_start = as.integer(ends),
        intron_end   = lead(as.integer(starts))
    ) %>%
    filter(!is.na(intron_end)) %>%
    transmute(chrom, start=intron_start, end=intron_end,
              isoform=name, geneName, strand)
write.table(introns, "hg38_intron_refFlat.txt", sep="\t", row.names=FALSE, quote=FALSE)

cat("Done — exon and intron files written\n")
EOF

# 15. install PyRanges
pip install pyranges

# 16. annotate regulatory (EPD + ENCODE cCREs + SEdb + TSS/SJ)

cd ~/TRpsyDB/scripts
python 07_annotate_regulatory.py \
    --db        data/trs.db \
    --epd       /home/rachele/TRpsyDB/data/raw/regulatory/Hs_EPDnew.bed \
    --ccre      /home/rachele/TRpsyDB/data/raw/regulatory/cCREs_hg38.bed \
    --sedb_se   /home/rachele/TRpsyDB/data/raw/regulatory/SE.bed \
    --sedb_te   /home/rachele/TRpsyDB/data/raw/regulatory/SE_te.bed \
    --sedb_meta /home/rachele/TRpsyDB/data/raw/regulatory/Human_sample_information_sedb3.txt \
    --refflat   /home/rachele/TRpsyDB/data/raw/refflat/refFlat.txt \
    --introns   /home/rachele/TRpsyDB/data/raw/refflat/hg38_intron_refFlat.txt \
    --exons     /home/rachele/TRpsyDB/data/raw/refflat/hg38_exon_refFlat.txt \
    --appris108 /home/rachele/TRpsyDB/data/raw/refflat/appris_data.principal.refseq108.hg38.txt \
    --appris109 /home/rachele/TRpsyDB/data/raw/refflat/appris_data.principal.refseq109.hg38.txt




# 17. ingest TRxQTL
python 08_ingest_trxqtl.py \
    --input /home/rachele/TRpsyDB/data/raw/trxqtl/ \
    --db    data/trs.db

    
# 18. ingest Tanudisastro sc-eTRs
python 09_ingest_tanudisastro.py \
    --supp1 data/raw/tanudisastro/supp_table1.xlsx \
    --supp3 data/raw/tanudisastro/supp_table3.xlsx \
    --supp5 data/raw/tanudisastro/supp_table5.xlsx \
    --db    data/trs.db

# 19. ingest Manigbas + Mukamel PheWAS
python 10_ingest_phenotype_assoc.py \
    --manigbas data/raw/manigbas2024_supp.xlsx \
    --mukamel  data/raw/mukamel2021_supp.xlsx \
    --db       data/trs.db

# 20. create turso database
turso db show trpsydb --url
libsql://trpsydb-rachelerr.aws-eu-west-1.turso.io
turso db tokens create trpsydb
eyJhbGciOiJFZERTQSIsInR5cCI6IkpXVCJ9.eyJpYXQiOjE3Nzk0NTczODAsImlkIjoiMDE5ZTRmZWItODUwMS03MDMyLWJjYTctYzk5ODUyNzk3MGIwIiwicmlkIjoiYTI5ZDBlZDAtMmI3ZS00NGJiLWE5MjgtY2ZhYmQ3MGFmZTkzIn0.36gF-_GAIpAl0cOqz97lq36ppHFwoQyL9qlg5fbyrNpjE_eboFdN_y2lHIpC-W9RatRibkUFSxkvBFT3r9L8DQ



