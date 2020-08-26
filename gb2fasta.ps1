$IN_GB_FILE = 'C:\Users\orlov\Desktop\norovirus\norovirus.gb'
$MIN_PARCER = 0 #nt
$MAX_PARCER = 8000 #nt
$MIN_REMOVE_SIMILAR = 0.5 #percentage e.g. '0.5'=0.5%
$MAX_REMOVE_SIMILAR = 100 #percentage e.g. '100'=100%
$PATH_GEN_ALIGNMENT = 'C:\Users\orlov\Documents\GitHub\GenAlignment\'
$PATH_RESOLVE_AMBIGUOUS = 'C:\Users\orlov\Documents\GitHub\resolve_ambiguous\'
$COORD_FILE = 'C:\Users\orlov\Desktop\norovirus\norovirus_orf.txt'

$path_dir = (($IN_GB_FILE.split('\'))[0..(($IN_GB_FILE.split('\')).GetUpperBound(0)-1)] -join '\') + '\'

python ($PATH_GEN_ALIGNMENT + 'parser_gb.py') -input $IN_GB_FILE -min $MIN_PARCER -max $MAX_PARCER -f 'country,host,collection_date'
python ($PATH_GEN_ALIGNMENT + 'remove_similar.py') -input ($IN_GB_FILE.split('.')[0] + '.fasta') -min $MIN_REMOVE_SIMILAR -max $MAX_REMOVE_SIMILAR
python ($PATH_RESOLVE_AMBIGUOUS + 'resolve_ambiguous.py') -input ($IN_GB_FILE.split('.')[0] + '_' + $MIN_REMOVE_SIMILAR.tostring() + '.fasta')
python ($PATH_GEN_ALIGNMENT + 'split_genome.py') -input ($IN_GB_FILE.split('.')[0] + '_' + $MIN_REMOVE_SIMILAR.tostring() + '_less_amb.fasta') -coord $COORD_FILE
python ($PATH_GEN_ALIGNMENT + 'trans_alignment.py') -input 
#дописать input (файлы по рамкам)
python ($PATH_GEN_ALIGNMENT + 'join_al.py') -input_list -output_name 
#дописать input_list и output_name
trimal -gt 0.8 -in -out
#дописать in и out 