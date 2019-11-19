import re

not_amb_line = re.compile(r"^[atgc\s]+$")

input_list = ["sapovirus_genomesORF1.fasta", "sapovirus_genomesORF1_less_amb.fasta"]

for input_file in input_list:
    with open(input_file, "r") as in_f:
        print("Filename: "+input_file+'\n')
        print("Sequence IDs with ambiguous nucleotides:"+'\n')
        count = 0
        file_lines = enumerate(in_f)
        for line_num, line in file_lines:
            if line.startswith('>'):
                iden = line
                amb = False
            else:
                if not not_amb_line.match(line):
                    if amb == True:
                        continue
                    amb = True
                    count += 1
                    print(iden+"Line: "+str(line_num+1))
    print("Total: "+str(count)+'\n' +'----------------------------------'+'\n')
