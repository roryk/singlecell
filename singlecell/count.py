from __future__ import print_function
from utils import file_transaction, file_exists
from collections import defaultdict, Counter
import os
import HTSeq

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
    """
    from the paper the identifier is
    BARCODE[6]-UMI[10]-[T or A 10]
    """
    identifier = read.read.name.split(":")[-1]
    return identifier[0:6], identifier[6:16]

