
#!/usr/bin/env python
from argparse import ArgumentParser, FileType
import sys


if __name__ == "__main__":
    parser = ArgumentParser(description=("Fix identifier of gtf_to_fasta "
                                         "converted FASTA files and produce a "
                                         "matching GTF file."))
    parser.add_argument('fasta_file', nargs='?', type=FileType('r'), default=sys.stdin)
    args = parser.parse_args()

    for line in args.infile:
        if line.startswith(">"):
            line = ">" + line.split()[1] + "\n"
        sys.stdout.write(line)
