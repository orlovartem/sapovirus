$PROJECT_DIR = "C:\Users\orlov\Documents\term_project\"
$GENOMES_FILE = "sapovirus_genomes.gb"
$UTR5 = "sapovirus_genomesUTR5"
$ORF1 = "sapovirus_genomesORF1"
$ORF2 = "sapovirus_genomesORF2"
$ORF12 = "sapovirus_genomesORF12"
$UTR3 = "sapovirus_genomesUTR3"
$LIST = $UTR5, $ORF1, $ORF2, $UTR3

python .\parcer_orf.py -input ($PROJECT_DIR + $GENOMES_FILE) -r

foreach ($ITEM in $LIST)
{
    if (($ITEM -eq $ORF1 -or $ITEM -eq $ORF2))
    {
        python .\resolve_ambiguous.py -input ($PROJECT_DIR + $ITEM + ".fasta")
        cat ($PROJECT_DIR + $ITEM + "_less_amb.fasta") | out-file -encoding ascii ($PROJECT_DIR + $ITEM + ".fasta")
    }
    if ($ITEM -eq $UTR5)
    {
        mafft ($PROJECT_DIR + $ITEM + ".fasta") | out-file -encoding ascii ($PROJECT_DIR + "alignments\align_UTR5.fasta")
    }
    if ($ITEM -eq $UTR3)
    {
        mafft ($PROJECT_DIR + $ITEM + ".fasta") | out-file -encoding ascii ($PROJECT_DIR + "alignments\align_UTR3.fasta")
    }

}



rm ($PROJECT_DIR + "*slices.fasta")
rm ($PROJECT_DIR + "local_db*")
rm ($PROJECT_DIR + "blast.out")
rm ($PROJECT_DIR + "*less_amb.fasta")