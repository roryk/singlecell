from argparse import ArgumentParser, FileType
import sys

if __name__ == "__main__":
    parser = ArgumentParser(description=("Create a dummy GTF file from a FASTA "
                                         "file of sequences."))
    parser.add_argument('infile', nargs='?', type=FileType('r'),
                       default=sys.stdin)
    args = parser.parse_args()

    gtf_format = '{feature}\t.\texon\t1\t{end}\t.\t.\t.\tgene_id "{feature}"\n'

    seq = ""
    feature = ""
    for line in args.infile:
        if line.startswith(">"):
            end = len(seq)
            if end:
                sys.stdout.write(gtf_format.format(**locals()))
            seq = ""
            end = 0
            feature = line[1:].strip()
        else:
            seq += line.strip()
