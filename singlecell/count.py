from __future__ import print_function
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.utils import file_exists
from singlecell.align import get_tag
from collections import defaultdict, Counter
import itertools
import os
import pysam
import HTSeq



def htseq_count(sam_file, gtf_file, multimappers=False):
    base, _ = os.path.splitext(sam_file)
    out_file = base + ".counted.sam"
    count_file = sam_file + ".counts"
    htseq_file = out_file + ".tmp"
    if file_exists(out_file):
        return out_file
    if multimappers:
        mm_flag = " --multimappers=yes "
    else:
        mm_flag = ""
    with file_transaction([count_file, htseq_file]) as tx_files:
        tx_count_file = tx_files[0]
        tx_htseq_file = tx_files[1]
        cmd = ("htseq-count --stranded=no --format=sam --samout={tx_htseq_file} "
               " --minaqual=0 --mode=intersection-nonempty "
               " {mm_flag} {sam_file} {gtf_file} > {tx_count_file}")
        message = "Count reads in {sam_file} mapping to {gtf_file}."
        do.run(cmd.format(**locals()), message.format(**locals()))
    fixed_file = fix_header(sam_file, htseq_file)
    os.rename(fixed_file, out_file)
    return out_file

def fix_header(orig_file, broken_file):
    "htseq-count returns a SAM file without a valid header on top, this fixes that"
    base, _ = os.path.splitext(broken_file)
    fixed_file = base + ".fixed" + ".sam"
    with file_transaction(fixed_file) as tx_fixed_file:
        with open(orig_file) as orig_handle, open(tx_fixed_file, "w") as fixed_handle:
            for line in orig_handle:
                if not line.startswith("@"):
                    break
                fixed_handle.write(line)
            with open(broken_file) as broken_handle:
                for line in broken_handle:
                    fixed_handle.write(line)
    return fixed_file

def count_reads(sam_file, feature_names, barcode_to_well):
    wells = sorted(barcode_to_well.values())
    base, _ = os.path.splitext(sam_file)
    out_file = base + ".counts"
    if file_exists(out_file):
        return out_file

    seen_umi = defaultdict(set)
    with pysam.Samfile(sam_file, "r", check_sq=False, check_header=False) as in_handle:
        for read in in_handle:
            aligned_id = get_tag(read, "XF")
            if discard_read(aligned_id):
                continue
            barcode, umi = get_barcode_and_umi(read)
            if barcode not in barcode_to_well:
                continue
            seen_umi[(aligned_id, barcode_to_well[barcode])].add(umi)
    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            print("\t".join(["feature"] + wells), file=out_handle)
            for feature in feature_names:
                counts = [len(seen_umi[(feature, well)]) for well in wells]
                print("\t".join([feature] + map(str, counts)), file=out_handle)

    return out_file

def count_umi(sam_file, gtf_file, barcode_to_well, multimappers=False):
    """
    stripped down implementation of the HTSeq algorithm for counting
    """
    base, _ = os.path.splitext(sam_file)
    out_file = base + ".counts"
    if file_exists(out_file):
        return out_file

    wells = sorted(barcode_to_well.values())
    seen_umi = defaultdict(set)
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    gtf_handle = HTSeq.GFF_Reader(gtf_file)
    for feature in gtf_handle:
        if feature.type == "exon":
            exons[feature.iv] += feature.attr["gene_id"]

    sam_handle = HTSeq.SAM_Reader(sam_file)
    for read in sam_handle:
        if not read.aligned:
            continue
        if not multimappers:
            try:
                if read.optional_field("NH") > 1:
                    continue
            except KeyError:
                pass
        iv_seq = (co.ref_iv for co in read.cigar if co.type == "M" and co.size > 0)
        fs = set()
        for iv in iv_seq:
            for iv2, fs2 in exons[iv].steps():
                if not fs:
                    fs = fs2.copy()
                else:
                    fs = fs.intersection(fs2)
        if len(fs) == 1:
            barcode, umi = get_barcode_and_umi(read)
            if barcode not in barcode_to_well:
                continue
            seen_umi[(list(fs)[0], barcode_to_well[barcode])].add(umi)

    with file_transaction(out_file) as tx_out_file:
        with open(tx_out_file, "w") as out_handle:
            print("\t".join(["feature"] + wells), file=out_handle)
            for feature in get_feature_names(gtf_file):
                counts = [len(seen_umi[(feature, well)]) for well in wells]
                print("\t".join([feature] + map(str, counts)), file=out_handle)

def get_feature_names(gtf_file):
    features = set()
    with open(gtf_file) as in_handle:
        for line in in_handle:
            info = line.split("\t")[8]
            info_fields = info.split(";")
            feature_field = filter(lambda x: x.startswith("gene_id"), info_fields)[0]
            feature = feature_field.split()[1].strip()
            features.add(feature.replace("\"", ""))
    return sorted(list(features))

def get_barcode_and_umi(read):
    identifier = read.read.name.split(":")[-1]
    return identifier[0:6], identifier[6:]

# def get_barcode_and_umi(read):
#     identifier = read.qname.split(":")[-1]
#     return identifier[0:6], identifier[6:]

def discard_read(aligned_id):
    if not aligned_id:
        return True
    # htseq-count reports unassignable reads with this prefix in the XF field
    if aligned_id.startswith("__"):
        return True
    return False

