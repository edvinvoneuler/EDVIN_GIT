#!/usr/bin/python
from sys import argv
from itertools import islice
from collections import defaultdict
import gzip

def remove_duplicates_paired_end(fwd_read, rev_read, basename):

    output_fwd = open(str(basename)+"_dedup_1.fastq", "w")
    output_rev = open(str(basename)+"_dedup_2.fastq", "w")
    forward_dict = defaultdict(list)

    with gzip.open(fwd_read, 'r') as fwd:
        with gzip.open(rev_read, 'r') as rev:
            while True:
                read_F = list(islice(fwd,4))
                read_R = list(islice(rev,4))
                if not (read_F and read_R):
                    break
                if not (read_F[1].strip() in forward_dict and read_R[1].strip() in forward_dict[read_F[1].strip()]):
                    forward_dict[read_F[1].strip()].append(read_R[1].strip())
                    output_fwd.write(read_F[0] + read_F[1] + read_F[2] + read_F[3])
                    output_rev.write(read_R[0] + read_R[1] + read_R[2] + read_R[3])
                    
            output_fwd.close()
            output_rev.close()
if __name__ == "__main__":
    fwd_read = argv[1]
    rev_read = argv[2]
    basename = argv[3]
    remove_duplicates_paired_end(fwd_read, rev_read, basename)
