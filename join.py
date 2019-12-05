from Bio import SeqIO
import argparse


def join(path_to_dir):
    dict_out = {}
    for el in ['UTR5', 'ORF1', 'ORF2', 'UTR3']:
        fasta_file = path_to_dir + "align_" + el + ".fasta"
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if seq_record.name in dict_out:
                dict_out[seq_record.name] += seq_record
            else:
                dict_out[seq_record.name] = seq_record
    with open("align_orfs_utrs.fasta", "w") as out_file:
        for name in dict_out:
            out_file.write(">" + name + "\n" + str(dict_out[name].seq) + "\n")

    dict_out = {}
    for el in ['ORF1', 'ORF2',]:
        fasta_file = path_to_dir + "align_" + el + ".fasta"
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if seq_record.name in dict_out:
                dict_out[seq_record.name] += seq_record
            else:
                dict_out[seq_record.name] = seq_record
    with open("align_orfs.fasta", "w") as out_file:
        for name in dict_out:
            out_file.write(">" + name + "\n" + str(dict_out[name].seq) + "\n")


    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-path", "--path_to_dir", type=str,
                        help="Path to directory", required=True)
    args = parser.parse_args()
    join(args.path_to_dir)