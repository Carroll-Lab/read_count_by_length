#!/usr/bin/env python

import csv
import argparse


def multi_get_seq(sample_names, out_file, min_read_len, max_read_len):
    sRNA_sizes = max_read_len + 1 - min_read_len
    single_seq_sRNA_count = []
    all_seq_sRNA_count = []
    
    samples = parse_file_names(sample_names)

    for i in range(len(samples)):
        single_seq_sRNA_count=[]
        for j in range(sRNA_sizes):
            single_seq_sRNA_count.append(0) #make an empty list fith sRNA_sizes * 0 counts
        all_seq_sRNA_count.append(single_seq_sRNA_count)


    seq_file_count = 0
    

    sRNA_len_count = {} #Dictionary of len:count
    read_count = 0
    file_count = 0

    for sample_name in samples:
        loaded_seq = open(sample_name, 'rU')
        file_count = 0
        sRNA_len_count = 0
        sRNA_len_count = {}
        for i in loaded_seq:

            a=i.strip()
            if a[0] == '>':
                b=a.split('-')
                read_count = int(b[1])
            
            else:
                sRNA_len = len(a)
                if min_read_len <= sRNA_len <= max_read_len and sRNA_len in sRNA_len_count:
                    sRNA_len_count[sRNA_len]+=read_count
                elif min_read_len <= sRNA_len <= max_read_len:
                    sRNA_len_count[sRNA_len] = read_count
                else:
                    pass

        for length, count in sRNA_len_count.viewitems():
            if min_read_len <= length <= max_read_len: 
                file_count += count
        for length, count in sRNA_len_count.viewitems():

            all_seq_sRNA_count[seq_file_count][length - min_read_len] = count*(1000000/float(file_count))

        seq_file_count+=1

    mycsv = csv.writer(open(out_file, 'wb'))
    
    headers = ['Length']
    for i in samples:
        headers.append(i)


    mycsv.writerow(headers)
    for i in range(sRNA_sizes):
        row=[i+min_read_len]
        for j in all_seq_sRNA_count: 
            row.append(j[i])
    
        mycsv.writerow(row)


def parse_file_names(in_list):
    sample_names=[]
    with open(in_list, 'rU') as loaded_list:

        for line in loaded_list:

            sample_names.append(line.strip())
    return sample_names


def comline():

    """
    Command line parser - add agruments here if needed
    """
    parser = argparse.ArgumentParser("read_count_by_length.py")
    parser.add_argument('in_list', type = str)
    parser.add_argument('out_file', type = str)
    parser.add_argument('-min','--min_read_len', type = int, default = 18)
    parser.add_argument('-max','--max_read_len', type=int, default = 32)
    args = parser.parse_args()
    return args

def main():
    """
    Run the analysis.  
    """

    a = comline()
    
    multi_get_seq(a.in_list, a.out_file, a.min_read_len,a.max_read_len)

if __name__ == "__main__":
    main()  
