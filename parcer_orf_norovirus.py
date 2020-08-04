from Bio import SeqIO
import argparse
import copy

def parse_gb(input_file, remove_exceptions):
    '''
    

    Parameters
    ----------
    input_file : string
        The name of file in genbank format.
    remove_exceptions : boolean
        When true, remove records which added in "exceptions_file" (see below).
        Without changes if false.

    Returns
    -------
    Three files in fasta format with sequences of CDS.
    If amount of CDS is not equal to 3 the record will be skipped during parcing.

    '''
    exceptions_file = 'norovirus_exceptions.csv'
    records_dict = {'orf1':[], 'orf2':[], 'orf3':[]}
    for record in SeqIO.parse(input_file, "genbank"):
        cds_dict={'orf1':False, 'orf2':False, 'orf3':False}
        
        full_seq = record.seq
        is_exception = False
        if remove_exceptions == True:
            with open(exceptions_file, 'r') as exceptions_f:
                for line in exceptions_f:
                    if line.strip() == record.name:
                        print('Record', record.name, 'was skipped')
                        is_exception = True
        if is_exception == True:
            continue
        for feature in record.features:
            if feature.type == 'CDS':
                for key in cds_dict:
                    if cds_dict[key] == False:
                        cds_dict[key] = True
                        seq = feature.location.extract(full_seq)
                        new_record = copy.deepcopy(record)
                        setattr(new_record, 'seq', seq)
                        records_dict[key].append(new_record)
                        break
    for key in records_dict:
        out_file_name = 'norovirus_'+key+'.fasta'
        out_file = open(out_file_name, 'w+')
        records_list = records_dict[key]
        SeqIO.write(records_list, out_file_name, 'fasta')
        out_file.close()
    
    print('\nDone')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-r", "--remove_exceptions",
                        help="Remove exceptions", action="store_true")
    args = parser.parse_args()
    parse_gb(args.input_file, args.remove_exceptions)