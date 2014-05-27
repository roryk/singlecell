import re
from utils import file_transaction, file_exists
import pysam
import os
import subprocess

accepted_barcode_pattern = re.compile(r"[ACGT]+[ACG]$")
polyA_tail_pattern = re.compile(r"A{20,}$")
MAX_EDIT_DISTANCE = 1
MAX_BEST = 10


def star_align(fastq_path, reference_prefix, out_prefix, cores=1):
    max_best = MAX_BEST
    out_file = out_prefix + "Aligned.out.sam"
    if file_exists(out_file):
        print ("%s has already been aligned, skipping." % (fastq_path))
        return out_file

    cmd = ("STAR --genomeDir {reference_prefix} --readFilesIn {fastq_path} "
           "--runThreadN {cores} --outFileNamePrefix {out_prefix} "
           "--outFilterMultimapNmax {max_best} "
           "--outSAMattributes NH HI NM MD AS "
           "--outSAMstrandField intronMotif").format(**locals())
    subprocess.check_call(cmd, shell=True)
    return out_file

def bwa_align(fastq_path, reference_prefix, out_file, cores=1):
    edit_distance = MAX_EDIT_DISTANCE
    if file_exists(out_file):
        print ("%s has already been aligned, skipping." % (fastq_path))
        return out_file

    with file_transaction(out_file) as tx_out_file:

        cmd = ("bwa aln -n {edit_distance} -l 24 {reference_prefix} "
               "{fastq_path} -t {cores} | bwa samse {reference_prefix} - {fastq_path} "
               "> {tx_out_file}").format(**locals())
        subprocess.check_call(cmd, shell=True)
    return out_file

def clean_align(align_file, out_file):
    seen = {}
    duped = {}
    if file_exists(out_file):
        print ("%s has already been UMI deduped, skipping." % (align_file))
        return out_file

    count_total_reads = 0
    count_assigned_reads = 0
    count_assigned_aligned_reads = 0
    with pysam.Samfile(align_file, "r") as in_handle, file_transaction(out_file) as tx_out_file:
        out_handle = pysam.Samfile(tx_out_file, "wh", template=in_handle)
        for read in in_handle:
            count_total_reads += 1
            count_assigned_reads += 1
            if poorly_mapped_read(read):
                continue
            count_assigned_aligned_reads += 1
            out_handle.write(read)
        out_handle.close

    return out_file

def unassigned_read(read, barcode_to_well):
    barcode = read.qname.split(":")[-1]
    if not re.match(accepted_barcode_pattern, barcode):
        return True
    well = barcode_to_well.get(barcode[0:6], None)
    if not well:
        return True
    return False

def poorly_mapped_read(read):
    if read.is_unmapped:
        return True
    nm = get_tag(read, "NM")
    if nm and (nm > MAX_EDIT_DISTANCE):
        return True
    x0 = get_tag(read, "X0")
    if x0 and (x0 > MAX_BEST):
        return True
    return False

def get_tag(read, tag):
    matched_tag = [x for x in read.tags if tag == x[0]]
    return matched_tag[0][1] if matched_tag else None
