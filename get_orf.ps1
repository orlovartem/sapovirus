$PROJECT_DIR = "C:\Users\orlov\Documents\term_project\"
$GENOMES_FILE = "sapovirus_genomes.gb"
$UTR5_FILE = "sapovirus_genomesUTR5.fasta"
$ORF1_FILE = "sapovirus_genomesORF1.fasta"
$ORF1_SLICES_FILE = "sapovirus_genomesORF1_slices.fasta"
$ORF1_LESSAMB_FILE = "sapovirus_genomesORF1_less_amb.fasta"
$ORF2_FILE = "sapovirus_genomesORF2.fasta"
$ORF12_FILE = "sapovirus_genomesORF12.fasta"
$UTR3_FILE = "sapovirus_genomesUTR3.fasta"

python .\parcer_orf.py -input ($PROJECT_DIR + $GENOMES_FILE) -r

python .\resolve_ambiguous.py -input ($PROJECT_DIR + $ORF1_FILE)

rm ($PROJECT_DIR + $ORF1_SLICES_FILE)
rm ($PROJECT_DIR + "local_db*")
rm ($PROJECT_DIR + "blast.out")

cat ($PROJECT_DIR + $ORF1_LESSAMB_FILE) | out-file -encoding ascii ($PROJECT_DIR + $ORF1_FILE)

rm ($PROJECT_DIR + $ORF1_LESSAMB_FILE)