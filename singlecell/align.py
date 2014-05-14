import re
from bcbio.provenance import do
from bcbio.utils import file_exists

accepted_barcode_pattern = re.compile(r"[ACGT]+[ACG]$")
polyA_tail_pattern = re.compile(r"A{20,}$")
max_edit_dist = 1
max_best = 10

def bwa_align(fastq_path, reference_prefix, alignment_dir):
    out_file = fastq_path + ".sam"
    if file_exists(out_file):
        return out_file

    cmd = ("bwa aln -l 24 {reference_prefix} {fastq_path} | "
           "bwa samse {reference_prefix} - {fastq_path} "
           "> {fastq_path}.sam").format(**locals())
    do.run(cmd, "Aligning %s to %s with bwa." % (fastq_path, reference_prefix),
           None)
    return out_file
