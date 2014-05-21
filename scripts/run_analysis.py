from argparse import ArgumentParser
from singlecell import barcode, align, count
import os
import glob
from collections import OrderedDict

def get_sample(line, sample_map_filename):
    keys = ["sample_id", "subsample_id", "r1_path", "r2_path"]
    sample_id, subsample_id, r1_filename, r2_filename = line.split()
    r1_path = os.path.join(os.path.dirname(sample_map_filename), r1_filename)
    r2_path = os.path.join(os.path.dirname(sample_map_filename), r2_filename)
    return dict(zip(keys, [sample_id, subsample_id, r1_path, r2_path]))

def get_samples_to_process(sample_file):
    with open(sample_file) as in_handle:
        return [get_sample(x, sample_file) for x in in_handle]

def get_r2_prepped_outfile(sample, alignment_dir):
    return os.path.join(alignment_dir,
                        ".".join([sample["sample_id"], sample["subsample_id"]]))

def get_aligned_outfile(fastq_file):
    return fastq_file + ".sam"

def get_cleaned_outfile(align_file):
    base, ext = os.path.splitext(align_file)
    return base + ".cleaned" + ext

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

def barcodes_to_plate_well(barcode_file):
    barcodes = OrderedDict()
    with open(barcode_file) as in_handle:
        for line in in_handle:
            tokens = line.split()
            barcodes[tokens[2]] = "_".join(tokens[0:2])
    return barcodes


if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
    parser.add_argument("--sample-map", required=True, help="Sample map file.")
    parser.add_argument("--bwa-ref", help="Path to bwa index")
    parser.add_argument("--alignment-dir", help="Output directory")
    parser.add_argument("--gtf-file", required=True, help="GTF file")
    parser.add_argument("--plate-file", required=True, help="Plate file")
    parser.add_argument("--cores", required=True, help="Number of cores to use.")
    args = parser.parse_args()

    samples = get_samples_to_process(args.sample_map)
    prepped = []
    # prep barcodes
    print "Beginning barcode preparation."
    for sample in samples:
        fq1 = sample["r1_path"]
        fq2 = sample["r2_path"]
        out_file = get_r2_prepped_outfile(sample, args.alignment_dir)
        print ("barcode-prepping %s and %s to %s." % (fq1, fq2, out_file))
        prepped.append(barcode.prep_r2_with_barcode(sample["r1_path"],
                                                    sample["r2_path"],
                                                    get_r2_prepped_outfile(sample, args.alignment_dir)))
    print "Finshed barcode preparation."

    print "Beginning alignment."
    aligned = []
    for prep in prepped:
        out_file = get_aligned_outfile(prep)
        print ("aligning %s to %s with bwa and writing to  %s." % (prep,
                                                                   args.bwa_ref,
                                                                   out_file))
        aligned.append(align.bwa_align(prep, args.bwa_ref, out_file))
    print "Finished alignment."

    print "Begin cleaning of poorly mapped reads."
    cleaned = []
    for sam_file in aligned:
        print ("Cleaning %s, removing poorly mapped reads." % align)
        cleaned.append(align.clean_align(sam_file, get_cleaned_outfile(sam_file)))
    print "Finished cleaning."

    print "Tagging reads that map to features in the GTF file."
    tagged = []
    for sam_file in cleaned:
        tagged.append(count.htseq_count(sam_file, args.gtf_file))
    print "Finished tagging reads."

    print "Reading feature names and barcodes."
    feature_names = get_feature_names(args.gtf_file)
    barcode_to_well = barcodes_to_plate_well(args.plate_file)
    print "Finished reading feature names and barcodes."

    print "Counting unique UMI mapping to features."
    counted = []
    for tag_file in tagged:
        counted.append(count.count_reads(tag_file, feature_names, barcode_to_well))
    print "Finished counting UMI."



