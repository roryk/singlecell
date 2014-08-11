from __future__ import print_function
from utils import file_transaction, safe_makedir, file_exists
import os
import itertools
import gzip

def format_fastq(buf):
    name, seq, qual = buf
    return "\n".join(["@" + name, seq, "+", qual])

def mask(seq, qual, min_qual=10):
    return "".join((b if (ord(q) - 33) >= min_qual else "N")
                   for b, q in itertools.izip(seq, qual))

def prep_r2_with_barcode(fq1, fq2, out_file):

    safe_makedir(os.path.dirname(out_file))
    if file_exists(out_file):
        print ("%s and %s have already been barcode-prepped, skipping."
               % (fq1, fq2))
        return out_file

    with open_fastq(fq1) as r1_file, open_fastq(fq2) as r2_file:
        with file_transaction(out_file) as tx_out_file:
            out_handle = open(tx_out_file, "w")
            read_count = 0
            buf = list()
            r1_r2 = itertools.izip(r1_file, r2_file)
            for header1, header2 in r1_r2:
                seq1, seq2 = r1_r2.next()
                plus1, plus2 = r1_r2.next()
                qual1, qual2 = r1_r2.next()

                read_name1, read_name2 = header1.split()[0][1:], header2.split()[0][1:]
                assert read_name1 == read_name2, "FASTQ files may be out of order."
                seq2, qual2 = seq2.rstrip(), qual2.rstrip()
                barcode, seq, qual = mask(seq1[0:6], qual1[0:6], min_qual=10) + \
                                     mask(seq1[6:], qual1[6:]), seq2, qual2
                barcoded_name = ":".join([read_name2, barcode])

                print(format_fastq([barcoded_name, seq, qual]), file=out_handle)
            out_handle.close()
    return out_file

def open_fastq(in_file):
    """ open a fastq file, using gzip if it is gzipped
    """
    _, ext = os.path.splitext(in_file)
    if ext == ".gz":
        return gzip.open(in_file, 'rb')
    if ext in [".fastq", ".fq"]:
        return open(in_file, 'r')

