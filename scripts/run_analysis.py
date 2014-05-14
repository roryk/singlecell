from argparse import ArgumentParser
from singlecell import barcode

def get_sample(line):
    keys = ["sample_id", "subsample_id", "r1_path", "r2_path"]
    sample_id, subsample_id, r1_filename, r2_filename = line.split()
    r1_path = os.path.join(os.path.dirname(sample_map_filename), r1_filename)
    r2_path = os.path.join(os.path.dirname(sample_map_filename), r2_filename)
    return dict(zip(keys, [sample_id, subsample_id, r1_path, r2_path]))

def get_samples_to_process(sample_file):
    with open(sample_file) as in_handle:
        return map(get_sample, in_handle)

def get_r2_prepped_outfile(sample, alignment_dir):
    return os.path.join(alignment_dir,
                        ".".join([sample["sample_id"], sample["subsample_id"]]))



if __name__ == "__main__":
    parser = ArgumentParser(description="Run a single cell analysis.")
    parser.add_argument("--sample-map", required=True, help="Sample map file.")
    parser.add_argument("--species", help="Species to use")
    parser.add_argument("--alignment-dir", help="Output directory")
    args = parser.parse_args()

    samples = get_samples_to_process(args.sample_map)
    for sample in samples:
        barcode.prep_r2_with_barcode(sample["r1_path"], sample["r2_path"],
                                     get_r2_prepped_outfile(sample, args.alignment_dir))




