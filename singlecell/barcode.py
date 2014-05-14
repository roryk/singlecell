from __future__ import print_function
from bcbio.distributed.transaction import file_transaction
from bcbio.bam import fastq
from bcbio.utils import file_exists, safe_makedir
import os

# accepted_barcode_pattern = re.compile(r"[ACGT]+[ACG]$")
# polyA_tail_pattern = re.compile(r"A{20,}$")
# max_edit_dist = 1
# max_best = 10

# def prep_barcodes(bamfile, out_dir, sample_id, ercc_dict, refseq_dict):
#     "read an alignment file and calculate a set of the UMI mapping to reads"
#     with bam.open_samfile(bamfile) as in_handle:
#         for read in in_handle:


# def skip_read(read):
#     barcode =
#     read.is_unmapped or re.match(accepted_barcode_pattern

def format_fastq(buf):
    name, seq, qual = buf
    return "\n".join(["@" + name, seq, "+", qual])

def mask(seq, qual, min_qual=10):
    return "".join((b if (ord(q) - 33) >= min_qual else "N")
                   for b, q in itertools.izip(seq, qual))

def prep_r2_with_barcode(fq1, fq2, out_file):

    safe_makedir(os.path.dirname(out_file))
    if file_exists(out_file):
        return out_file

    with bam.fastq.open_fastq(fq1) as r1_file, bam.fastq.open_fastq(fq2) as r2_file:
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
                                     mask(seq1[6:r1_length], qual1[6:r1_length]), seq2, qual2
                barcoded_name = ":".join([read_name2, barcode])

                print(format_fastq([barcoded_name, seq, qual]), file=out_handle)
            out_handle.close()
    return out_file
