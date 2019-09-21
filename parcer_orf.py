import argparse
import csv
import os
import re
import sys
from textwrap import fill

def parse_gb(input_file, min_length, max_length):
    '''
    J:
    Description of function.
    Input:
    variable - type - description
    Output:
    variable - type - description
    '''
    OUTPUT_FILE_UTR5 = '.'.join(input_file.split('.')[:-1]) + 'UTR5.fasta'
    OUTPUT_FILE_ORF1 = '.'.join(input_file.split('.')[:-1]) + 'ORF1.fasta'
    OUTPUT_FILE_ORF2 = '.'.join(input_file.split('.')[:-1]) + 'ORF2.fasta'
    OUTPUT_FILE_ORF12 = '.'.join(input_file.split('.')[:-1]) + 'ORF12.fasta'
    OUTPUT_FILE_UTR3 = '.'.join(input_file.split('.')[:-1]) + 'UTR3.fasta'

    MIN_ORIGIN_SIZE = min_length
    MAX_ORIGIN_SIZE = max_length

    accession = re.compile(r"^ACCESSION\s+([a-z_A-Z0-9]+)")  # accession number
    cds = re.compile(r"^\s{5}CDS[0-9\.\s<>]+")
    utr5 = re.compile(r"^\s{5}5\'UTR[0-9\.\s<>]+")
    utr3 = re.compile(r"^\s{5}3\'UTR[0-9\.\s<>]+")
    '''
    J:
    Artem, please read about the ambiguous nucleotides:
    https://genomevolution.org/wiki/index.php/Ambiguous_nucleotide
    They are very common. Your RegExp might loose the nucleotides in 
    unexpected positions. You'd better use [atgcnrykmswbdhvu\s]+ here.
    '''
    origin = re.compile(r"\s+\d+\s+([a,t,g,c ]+)$") 
    orf1 = re.compile(r".*(major|orf1|ORF1|VP1|polyprotein).*")
    orf2 = re.compile(r".*(minor|orf2|ORF2|VP2).*")


    def take_location(line_in):
        '''
        J:
        Description of the function.
        line_in - type - description
        [start_loc, end_loc] - type - descr
        '''
        loc = ((line_in.split())[1]).split("..")
        start_loc = int((re.findall(r"\w+", loc[0]))[0])
        end_loc = int((re.findall(r"\w+", loc[1]))[0])
        return [start_loc, end_loc]

    try:
        utr5_location = [0, 0]
        cds_location = [0, 0]
        orf1_location = [0, 0]
        orf2_location = [0, 0]
        utr3_location = [0, 0]
        out_utr5 = open(OUTPUT_FILE_UTR5, 'w')
        out_orf1 = open(OUTPUT_FILE_ORF1, 'w')
        out_orf2 = open(OUTPUT_FILE_ORF2, 'w')
        out_orf12 = open(OUTPUT_FILE_ORF12, 'w')
        out_utr3 = open(OUTPUT_FILE_UTR3, 'w')
        with open(input_file, "r") as in_f:
            cds_field = False  # outside cds field
            note_field = False
            product_field = False
            gene_field = False
            test_origin = ''
            for line in in_f:

                m = accession.match(line)  # finds ACCESSION field using RegExp
                if (m):
                    """
                    J:
                    min and max size of sequences can not be used since you write
                    every sequence into output file
                    """
                    test_accession = m.group(1)  # accession number
                    out_utr5.write('>'+test_accession+'\n')
                    out_orf1.write('>'+test_accession+'\n')
                    out_orf2.write('>'+test_accession+'\n')
                    out_orf12.write('>'+test_accession+'\n')
                    out_utr3.write('>'+test_accession+'\n')

                if re.match(r"^\s{5}[a-z_A-Z]+[\s0-9<>\.]+", line) or re.match(r"^[a-zA-z]+", line):
                    cds_field = False  # outside cds field

                # if re.match(r"\s{5}source[0-9\.\s<>]+", line):
                #    source_field = True #inside source field

                if utr5.match(line):
                    utr5_location = take_location(line)

                if cds.match(line):
                    cds_field = True  # inside cds field
                    cds_location = take_location(line)


                if utr3.match(line):
                    utr3_location = take_location(line)

                m = origin.match(line)

                if origin.match(line):
                    test_origin += re.sub("\s+","",m.group(1))

                if cds_field is True:
                    if re.match(r"^\s+/[a-zA-Z]+.*", line):
                        """
                        J:
                        Why not making note_field,product_field,product_field False 
                        when you reach the end of entry? line 130
                        """
                        note_field = False
                        product_field = False
                        product_field = False
                    if re.match(r"^\s+/note=.+", line) :
                        note_field = True  # inside /note
                    if re.match(r"^\s+/product=.+", line) :
                        product_field = True  # inside /product
                    if re.match(r"^\s+/gene=.+", line) :
                        gene_field = True  # inside /gene
                    if (note_field or product_field or gene_field) is True and orf1.match(line):
                        orf1_location = cds_location
                    if (note_field or product_field or gene_field) is True and orf2.match(line):
                        orf2_location = cds_location

                if re.match(r"^//$", line):  # end of record
                    """
                    J:
                    You might check the same for UTR3
                    """
                    if utr5_location == [0, 0]:
                        utr5_location = [0, orf1_location[0]-1]
                    print (test_accession, utr5_location, orf1_location, orf2_location, utr3_location)
                    out_utr5.write(fill(test_origin[utr5_location[0]-1:utr5_location[1]], 60)+'\n')
                    out_orf1.write(fill(test_origin[orf1_location[0]-1:orf1_location[1]], 60)+'\n')
                    out_orf2.write(fill(test_origin[orf2_location[0]-1:orf2_location[1]], 60)+'\n')
                    out_orf12.write(fill(test_origin[orf1_location[0]-1:orf2_location[1]], 60)+'\n')
                    out_utr3.write(fill(test_origin[utr3_location[0]-1:utr3_location[1]], 60)+'\n')
                    test_origin = ''
                    utr5_location = [0, 0]
                    cds_location = [0, 0]
                    orf1_location = [0, 0]
                    orf2_location = [0, 0]
                    utr3_location = [0, 0]

        out_utr5.close()
        out_orf1.close()
        out_orf2.close()
        out_orf12.close()
        out_utr3.close()

    except:
        print('Error: can\'t read genbank file')

    return 'Parced'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    """
    J:
    If the argument is obligatory, you can add  
    required=True to add_argument() function. My solution with checking
    the number of arguments was lame.
    """
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-min", "--min_length", type=int,
                        help="Minimal length of sequence.\
                        Sequences shorter than min length \
                        will not be included in the final dataset")
    parser.add_argument("-max", "--max_length", type=int,
                        help="Maximal length of sequence. \
                        Sequences longer than max length \
                        will not be included in the final dataset")
    args = parser.parse_args()


    parse_gb(args.input_file, args.min_length, args.max_length)
