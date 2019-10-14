import argparse
import csv
import os
import re
import sys
from textwrap import fill


def parse_gb(input_file):
    '''
    The function parses genbank files and formes files with open reading frames
    and untranslated regions.

    Input:
    input_file - string - the name of genbank file in the folder of this script

    Output:
    without variable
    5 output files with ORF and UTR
    '''

    OUTPUT_FILE_UTR5 = '.'.join(input_file.split('.')[:-1]) + 'UTR5.fasta'
    OUTPUT_FILE_ORF1 = '.'.join(input_file.split('.')[:-1]) + 'ORF1.fasta'
    OUTPUT_FILE_ORF2 = '.'.join(input_file.split('.')[:-1]) + 'ORF2.fasta'
    OUTPUT_FILE_ORF12 = '.'.join(input_file.split('.')[:-1]) + 'ORF12.fasta'
    OUTPUT_FILE_UTR3 = '.'.join(input_file.split('.')[:-1]) + 'UTR3.fasta'

    accession = re.compile(r"^ACCESSION\s+([a-z_A-Z0-9]+)")
    # accession number
    cds = re.compile(r"^\s{5}CDS[0-9\.\s<>]+")
    utr5 = re.compile(r"^\s{5}5\'UTR[0-9\.\s<>]+")
    utr3 = re.compile(r"^\s{5}3\'UTR[0-9\.\s<>]+")
    origin = re.compile(r"\s+\d+\s+([atgcnrykmswbdhvu\s]+)$")
    orf1 = re.compile(r".*(major|orf1|ORF1|VP1|polyprotein).*")
    orf1a = re.compile(r".*ORF1a.*")
    orf1b = re.compile(r".*ORF1b.*")
    orf2 = re.compile(r".*(minor|orf2|ORF2|VP2).*")
    source = re.compile(r"^\s{5}source[0-9\.\s<>]+")

    def take_from_commas(line):
        '''
        The function takes string value between commas in line.

        Input:
        line - string - an input line

        Output:
        result - string - symbols between two commas in line
        '''
        count = 0
        i1 = 0
        i2 = 0
        for i in line:
            count += 1
            if i == '\"':
                if not i1 == 0:
                    i2 = count
                else:
                    i1 = count
        result = line[i1:i2-1]
        return result

    def csv_reader(csv_file_input, input_str):
        '''
        Reads csv-file and outputs the short name of a variable.

        Input:
            csv_file_input - string - name of csv-file in the directory of the
            script input_str - string - a string value from feature in genbank file

        Output:
            output_str - string - short designation from csv-file
        '''
        with open(csv_file_input) as csv_file:
            output_str = ''
            reader = csv.DictReader(csv_file, delimiter=",", fieldnames=["base", "new"])
            for line in reader:
                if (line["base"] or line["new"]) in input_str:
                    output_str = line["new"]
            return output_str

    def take_location(line_in):
        '''
        The function takes location of ORF or UTR.

        Input:
        line_in - string - a line that contains location of ORF or UTR

        Output:
        [start_loc, end_loc] - integer (list) - location of ORF or UTR
        (example: [23, 1035])

        '''
        loc = ((line_in.split())[1]).split("..")
        start_loc = int((re.findall(r"\w+", loc[0]))[0])
        end_loc = int((re.findall(r"\w+", loc[1]))[0])
        return [start_loc, end_loc]
    try:
        count_total_entries = 0
        count_notread = 0
        source_location = [0, 0]
        utr5_location = [0, 0]
        cds_location = [0, 0]
        orf1_location = [1, 0]
        orf2_location = [1, 0]
        orf1a_location = [0, 0]
        orf1b_location = [0, 0]
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
                    test_accession = m.group(1)  # accession number

                if re.match(r"^\s{5}[a-z_A-Z]+[\s0-9<>\.]+", line) or (
                        re.match(r"^[a-zA-z]+", line)):
                    cds_field = False  # outside cds field

                if source.match(line):
                    source_location = take_location(line)

                if utr5.match(line):
                    utr5_location = take_location(line)

                if cds.match(line):
                    cds_field = True  # inside cds field
                    cds_location = take_location(line)

                if utr3.match(line):
                    utr3_location = take_location(line)

                m = origin.match(line)

                if origin.match(line):
                    test_origin += re.sub("\s+", "", m.group(1))

                if re.match(r"^\s+/strain=.+", line):
                    strain = take_from_commas(line)

                if re.match(r"^\s+/country=.+", line):
                    country = take_from_commas(line)
                    if not country == '':
                        country = csv_reader('country_map.csv', country)

                if re.match(r"^\s+/host=.+", line):
                    host = take_from_commas(line)
                    if not host == '':
                        host = csv_reader('host_map.csv', host)

                if re.match(r"^\s+/collection_date=.+", line):
                    collection_date = take_from_commas(line)

                if cds_field is True:
                    if re.match(r"^\s+/[a-zA-Z]+.*", line):
                        note_field = False
                        product_field = False
                    if re.match(r"^\s+/note=.+", line):
                        note_field = True  # inside /note
                    if re.match(r"^\s+/product=.+", line):
                        product_field = True  # inside /product
                    if re.match(r"^\s+/gene=.+", line):
                        gene_field = True  # inside /gene
                    if (note_field or product_field or gene_field) is True and orf1a.match(line):
                        orf1a_location = cds_location
                    if (note_field or product_field or gene_field) is True and orf1b.match(line):
                        orf1b_location = cds_location
                    if (note_field or product_field or gene_field) is True and orf1.match(line):
                        orf1_location = cds_location

                    if (note_field or product_field or gene_field) is True and orf2.match(line):
                        orf2_location = cds_location
                m = accession.match(line)  # finds ACCESSION field using RegExp
                if (m):
                    test_accession = m.group(1)  # accession number

                if re.match(r"^//$", line):  # end of record
                    out_utr5.write('>'+test_accession+'_'+strain+'_'+country+'_'+host+'_'+collection_date+'\n')
                    out_orf1.write('>'+test_accession+'_'+strain+'_'+country+'_'+host+'_'+collection_date+'\n')
                    out_orf2.write('>'+test_accession+'_'+strain+'_'+country+'_'+host+'_'+collection_date+'\n')
                    out_orf12.write('>'+test_accession+'_'+strain+'_'+country+'_'+host+'_'+collection_date+'\n')
                    out_utr3.write('>'+test_accession+'_'+strain+'_'+country+'_'+host+'_'+collection_date+'\n')

                    count_total_entries += 1
                    if not orf1a_location == [0, 0]:
                        orf1_location = [orf1a_location[0], orf1b_location[1]]
                    if utr3_location == [0, 0]:
                        if not orf2_location[1] == 0:
                            utr3_location = [orf2_location[1], source_location[1]]
                        else:
                            utr3_location = [orf1_location[1], source_location[1]]
                    if utr5_location == [0, 0]:
                        utr5_location = [1, orf1_location[0]-1]
                    if orf1_location == [1, 0] and orf2_location == [1, 0]:
                        count_notread += 1
                        utr3_location = [1, 0]
                    # print(test_accession, utr5_location, orf1_location, orf2_location, utr3_location)
                    out_utr5.write(fill(test_origin[utr5_location[0]-1:utr5_location[1]], 60)+'\n')
                    out_orf1.write(fill(test_origin[orf1_location[0]-1:orf1_location[1]], 60)+'\n')
                    out_orf2.write(fill(test_origin[orf2_location[0]-1:orf2_location[1]], 60)+'\n')
                    out_orf12.write(fill(test_origin[orf1_location[0]-1:orf2_location[1]], 60)+'\n')
                    out_utr3.write(fill(test_origin[utr3_location[0]-1:utr3_location[1]], 60)+'\n')
                    test_origin = ''
                    utr5_location = [0, 0]
                    cds_location = [0, 0]
                    orf1a_location = [0, 0]
                    orf1b_location = [0, 0]
                    orf1_location = [1, 0]
                    orf2_location = [1, 0]
                    utr3_location = [0, 0]
                    cds_field = False  # outside cds field
                    note_field = False
                    product_field = False
                    gene_field = False
                    country = 'unknown-country'
                    host = 'unknown-host'
        print('Total entries:', count_total_entries, '\nNot readable entries:', count_notread)
        out_utr5.close()
        out_orf1.close()
        out_orf2.close()
        out_orf12.close()
        out_utr3.close()

    except Exception:
        print('Error: can\'t read genbank file')

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    args = parser.parse_args()

    parse_gb(args.input_file)
