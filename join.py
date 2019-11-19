dict_seq={}
with open("align_UTR5.fasta", "r")as in_f:
    iden = 'tmp'
    seq = ''
    for line in in_f:
        if line.startswith('>'):
            dict_seq[iden]=seq
            iden = line
            seq = ''
        else:
            seq += line.strip()
    dict_seq[iden]=seq
    #print(dict_seq['>CS790749_Mc114_unknown-country_unknown-host_Mar-2000\n'])

#for el in ['ORF1','ORF2','UTR3']:
for el in ['ORF1_less_amb','ORF2','UTR3']:
    with open("align_"+el+".fasta", "r")as in_f:
        iden = 'tmp'
        seq = ''
        for line in in_f:
            if line.startswith('>'):
                dict_seq[iden]+=seq
                iden = line
                seq = ''
            else:
                seq += line.strip()

#with open("align_genomes.fasta", "w") as out_f:
with open("align_genomes_less_amb.fasta", "w") as out_f:
    del dict_seq['tmp']
    for iden in dict_seq:
        out_f.write(iden+dict_seq[iden]+'\n')