from argparse import ArgumentParser
from singlecell import barcode, align, count, cluster
from singlecell.barcode import barcodes_to_plate_well
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

def get_bwa_outfile(fastq_file):
    return fastq_file + ".sam"

def get_star_prefix(fastq_file):
    base, _ = os.path.splitext(fastq_file)
    return base

def get_cleaned_outfile(align_file):
    base, ext = os.path.splitext(align_file)
    return base + ".cleaned" + ext

    barcodes = OrderedDict()
    with open(barcode_file) as in_handle:
        for line in in_handle:
            tokens = line.split()
            barcodes[tokens[2]] = "_".join(tokens[0:2])
    return barcodes

if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
    parser.add_argument("--multimappers", action="store_true",
                        default=False, help="Keep multimappers")
    parser.add_argument("--sample-map", required=True, help="Sample map file.")
    parser.add_argument("--aligner", default="bwa",
                        choices=["bwa", "star"], help="Aligner to use.")
    parser.add_argument("--aligner-index", help="Path to aligner index.")
    parser.add_argument("--alignment-dir", help="Output directory")
    parser.add_argument("--gtf-file", required=True, help="GTF file")
    parser.add_argument("--plate-file", required=True, help="Plate file")
    parser.add_argument("--num-jobs", type=int,
                        default=1, help="Number of concurrent jobs to process.")
    parser.add_argument("--cores-per-job", type=int,
                        default=1, help="Number of cores to use.")
    parser.add_argument("--memory-per-job", default=2, help="Memory in GB to reserve per job.")
    parser.add_argument("--timeout", default=15, help="Time to wait before giving up starting.")
    parser.add_argument("--scheduler", default=None, help="Type of scheduler to use.",
                        choices=["lsf", "slurm", "torque", "sge"])
    parser.add_argument("--resources", default=None, help="Extra scheduler resource flags.")
    parser.add_argument("--queue", default=None, help="Queue to submit jobs to.")
    parser.add_argument("--local", action="store_true", default=False,
                        help="Run in parallel on a local machine.")

    args = parser.parse_args()

    samples = get_samples_to_process(args.sample_map)
    prepped = []

    print "Starting IPython cluster. This may take a while."
    with cluster.get_cluster_view(args) as view:
        print "IPython cluster is up."

        print "Beginning barcode preparation."
        for sample in samples:
            fq1 = sample["r1_path"]
            fq2 = sample["r2_path"]
            out_file = get_r2_prepped_outfile(sample, args.alignment_dir)
            print ("barcode-prepping %s and %s to %s." % (fq1, fq2, out_file))
            prepped.append(view.apply_async(barcode.prep_r2_with_barcode,
                                            sample["r1_path"],
                                            sample["r2_path"],
                                            get_r2_prepped_outfile(sample, args.alignment_dir)))
        prepped = cluster.wait_until_complete(prepped)
        print "Finshed barcode preparation."

        print "Beginning alignment."
        aligned = []
        for prep in prepped:
            print ("aligning %s to %s with %s." % (prep, args.aligner_index,
                                                   args.aligner))
            if args.aligner == "bwa":
                aligned.append(view.apply_async(align.bwa_align, prep, args.aligner_index,
                                                   get_bwa_outfile(prep), args.cores_per_job))
            elif args.aligner == "star":
                aligned.append(view.apply_async(align.star_align, prep, args.aligner_index,
                                                get_star_prefix(prep), args.cores_per_job))
        aligned = cluster.wait_until_complete(aligned)
        print "Finished alignment."

        print "Begin cleaning of poorly mapped reads."
        cleaned = []
        for sam_file in aligned:
            print ("Cleaning %s, removing poorly mapped reads." % sam_file)
            cleaned.append(view.apply_async(align.clean_align, sam_file, get_cleaned_outfile(sam_file)))
        cleaned = cluster.wait_until_complete(cleaned)
        print "Finished cleaning."

        print "Reading barcode to well mapping."
        barcode_to_well = barcodes_to_plate_well(args.plate_file)
        print "Finished reading feature names and barcodes."

        print "Counting unique UMI mapping to features."
        counted = []
        for sam_file in cleaned:
            counted.append(view.apply_async(count.count_umi, sam_file, args.gtf_file, barcode_to_well))
        counted = cluster.wait_until_complete(counted)
        print "Finished counting UMI."
