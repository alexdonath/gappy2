#!/usr/bin/env python3
import argparse
import math
import os
import re
import sys

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from intervaltree import IntervalTree


def get_args(argv=None):
    '''Parses the command line options.\n
    Returns populated namespace.
    '''
    parser = argparse.ArgumentParser(description="gappy2 extracts splids \
        (split-inducing indels) from multiple sequence alignments.")
    parser.add_argument("-f", "--fasta",
                        help="input alignment file in FASTA format")
    parser.add_argument("-m", "--maf",
                        help="input alignment file in MAF format")
    parser.add_argument("-o", "--output-prefix", default="GAPPY2",
                        help="prefix for the output files [%(default)s]")
    parser.add_argument("-l", "--minlength", default=2, type=int,
                        help="minimum indel size [%(default)s]")
    parser.add_argument("-u", "--maxlength", default=math.inf, type=int,
                        help="maximum indel size [%(default)s]")
    parser.add_argument("-s", "--length", metavar="LENGTH", default=None,
                        type=int,
                        help="search for indels of exactly size [not set]")
    parser.add_argument("-c", "--unknownchar", metavar='CHAR', default="?",
                        help="character of unknown sites [%(default)s]")
    parser.add_argument("-z", "--fuzzy", default=False, action="store_true",
                        help="use fuzzy search [%(default)s]")
    parser.add_argument("-v", "--version", action="store_true",
                        help="output version information and exit")

    return parser.parse_args(argv)


def check_input(args):
    '''Checks provided arguments for consistency.\n
    Returns a consistent set of arguments.
    '''

    if args.version:
        sys.exit(print_version())

    if args.fasta and args.maf:
        sys.exit("[ERROR] MAF and FASTA file provided. Please select one.")

    if args.length:
        print("[INFO] Exact splid length provided. Set as min and max length.")
        args.minlength = args.length
        args.maxlength = args.length

    if args.minlength > args.maxlength:
        sys.exit("[ERROR] Minimum length larger than maximum length.")

    if len(args.unknownchar) != 1:
        sys.exit("[ERROR] Please select a single character for unknown sites.")

    # print(args)


def report_parameters(args):
    '''Prints the arguments which are used for the analysis.
    '''

    if args.fasta:
        print("[INFO] {:20s}{}".format("Input file: ", args.fasta))
        print("[INFO] {:20s}{}".format("File type: ", "FASTA"))
    elif args.maf:
        print("[INFO] {:20s}{}".format("Input file: ", args.maf))
        print("[INFO] {:20s}{}".format("File type: ", "MAF"))

    print("[INFO] {:20s}{}".format("Output prefix: ", args.output_prefix))
    print("[INFO] {:20s}{}".format("Minimum length:", args.minlength))
    print("[INFO] {:20s}{}".format("Maximum length: ", args.maxlength))
    print("[INFO] {:20s}{}".format("Unknown character:", args.unknownchar))
    print("[INFO] {:20s}{}".format("Fuzzy search:", args.fuzzy))
    print("[INFO] {:20s}{}".format("Writing PHYLIP:", args.phylip))


def delete_SNPs(tree):
    ''' Modifies an interval tree by removing single residue entries.
    '''
    items = list(sorted(tree))
    for item in items:
        if (item.end - item.begin) == 1:
            # print("SNP:",item)
            tree.remove(item)


def find_stretches(alignments, character):
    '''Finds occurrences of a character in an alignment and builds up
    an interval tree from their start and stop positions.\n
    Returns the interval tree.
    '''
    tree = IntervalTree()
    for sequence in alignments:
        # print(sequence.seq)
        find_this = re.compile(r"{}+".format(character))
        for m in re.finditer(find_this, str(sequence.seq)):
            tree.addi(m.start(), m.end(), sequence.id)
    return(tree)


def identify_candidates(tree):
    '''Identifies splid candidates in an interval tree.\n
    Candidates are intervals that are present at least twice.\n
    Returns a set of intervals.
    '''

    gaps = list(sorted(tree))
    candidates = set()
    conflicts = set()

    for gap in gaps:
        # print(gap)
        overlaps = tree.search(gap.begin, gap.end)

        for overlap in overlaps:
            # print(gap, "versus", overlap)
            if gap == overlap:
                # print("similar:", gap, overlap)
                continue
            if (gap.begin == overlap.begin) and (gap.end == overlap.end):
                candidates.add(gap)
                candidates.add(overlap)
                continue
            else:
                # print("Found conflict between",gap,"and",overlap)
                conflicts.add(overlap)

        # if conflicts :
            # print("Conflict:",conflicts)

        # if candidates:
            # print("Good", candidates)
    return candidates


def identify_splids(tree):
    '''Identifies splids in an interval tree.\n
    Returns an interval tree where all intervals are splids.
    '''

    candidates = list(sorted(tree))
    conflicts = set()

    for candidate in candidates:
        overlaps = tree.search(candidate.begin, candidate.end)

        for overlap in overlaps:
            # print(overlap, "versus", candidate)
            if candidate.begin == overlap.begin and candidate.end == overlap.end:
                continue
            else:
                # print("Found conflict between",candidate,"and",overlap)
                conflicts.add(overlap)
                conflicts.add(candidate)

    for i in conflicts:
        # tree.remove(i)
        tree.discard(i)

    return tree


def get_final_list(tree):
    '''Generates a dictionary that contains splid intervals as keys and the
    corresponding sequences as values.\n
    Returns the dictionary.
    '''
    final_list = {}

    for interval_obj in sorted(tree):
        start_stop = "{}-{}".format(interval_obj.begin, interval_obj.end)
        if start_stop in final_list:
            final_list[start_stop] += str("," + interval_obj.data)
        else:
            final_list[start_stop] = str(interval_obj.data)

    return final_list


def print_tsv(final_list, count, outtsv):
    '''Prints final intervals. Output is in BED style:\n
    Counting starts at 0.\n
    The starting position is included.\n
    The ending position is not included.\n
    '''
    outtsv.write("#aln_no\tstart\tstop\tsequences\tlength\n")
    for i, j in final_list.items():
        coordinates = [int(x) for x in i.split("-")]
        length = coordinates[1]-coordinates[0]
        outtsv.write("{:d}\t{:d}\t{:d}\t{:s}\t{:d}\n".format(
            count, coordinates[0], coordinates[1], j, length))


def create_splid_alignment(alignment, final_list, unknowntree, count):
    '''Generates presence/absence matrix.\n
    Returns a MultipleSeqAlignment object.
    '''
    record_list = []

    for record in alignment:
        newseq = SeqRecord(MutableSeq(""),
                           id=record.id,
                           description=record.description)
        record_list.append(newseq)

    align = MultipleSeqAlignment(record_list,
                                 annotations={"tool": "gappy2",
                                              "alphabet": "binary"})

    for i, j in final_list.items():

        # overlap with unknownchar
        [start, stop] = [int(x) for x in i.split("-")]
        unknown_overlap = unknowntree.search(start, stop)
        unknownseqs = []
        for region in unknown_overlap:
            unknownseqs.append(region.data)

        header = j.split(sep=',')
        for record in align:
            if record.id in header:
                record.seq += "1"
            elif record.id in unknownseqs:
                record.seq += "?"
            else:
                record.seq += "0"

    return align


def print_tree(tree):
    '''Prints an interval tree.
    '''
    print(tree)


def get_infile_basename(infile_name):
    '''Gets the basename of the input file.\n
    Returns the basename.
    '''
    infile_name = os.path.basename(infile_name)
    return infile_name


# TODO:
# * Do not overwrite existing files.
def get_outfile_name(prefix, infile_name, min, max, fuzzy):
    '''Creates the output filename.\n
    Returns the output filename.
    '''
    outfilename = prefix + "_" + infile_name + "_" + str(min) + "-" + str(max)
    if (fuzzy):
        outfilename += "_z"
    return outfilename


def filter_interval_size(tree, min=2, max=math.inf):
    '''Removes intervals from an interval tree that have a length outside
    a given min/max range.\n
    Returns an interval tree with intervals that have a length within a given
    min/max range.
    '''
    candidates = tree.items()
    for candidate in candidates:
        length = candidate.end - candidate.begin
        if length > max or length < min:
            tree.remove(candidate)
    return tree


def print_version():
    '''Prints gappy2 version.
    '''
    ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
    try:
        with open(os.path.join(ROOT_DIR, 'VERSION')) as version_file:
            version = version_file.read().strip()
    except FileNotFoundError:
        version = "0.1.0"

    print(version)


# MAIN
if __name__ == "__main__":
    ''' TODO:
    * Exception handling
    * Clean up code
    '''
    args = get_args()

    check_input(args)

    # report_parameters(args)

    alignments = {}

    input_filename = ""

    if args.fasta:
        try:
            alignments = {AlignIO.read(args.fasta, "fasta")}
            input_filename = get_infile_basename(args.fasta)
        except ValueError as err:
            print("[ERROR] Please check FASTA input.", err)
            sys.exit()
    elif args.maf:
        try:
            alignments = AlignIO.parse(args.maf, "maf")
            input_filename = get_infile_basename(args.maf)
        except ValueError as err:
            print("[ERR] Please check MAF input.", err)
            sys.exit()

    outfile_name = get_outfile_name(prefix=args.output_prefix,
                                    infile_name=input_filename,
                                    min=args.minlength,
                                    max=args.maxlength,
                                    fuzzy=args.fuzzy)
    outfile_tsv = outfile_name + ".tsv"
    outtsv = open(outfile_tsv, "w")

    counter = 0
    splid_alignments = []

    print("Number of checked alignments:", counter, sep='\n', end='')
    for alignment in alignments:

        initial_tree = find_stretches(alignment, re.escape("-"))

        unknowntree = find_stretches(alignment, re.escape(args.unknownchar))

        if (args.fuzzy):
            delete_SNPs(initial_tree)

        possible_candidates = identify_candidates(initial_tree)

        candidate_tree = IntervalTree(possible_candidates)

        splids_tree = identify_splids(candidate_tree)

        valid_splids = filter_interval_size(splids_tree,
                                            min=args.minlength,
                                            max=args.maxlength)

        # splid identification done. now we need to print the results.
        counter += 1

        final_splid_list = {}
        if not valid_splids.is_empty():
            final_splid_list = get_final_list(valid_splids)

        if final_splid_list.items():
            print_tsv(final_splid_list, counter, outtsv)
            splid_alignment = create_splid_alignment(alignment,
                                                     final_splid_list,
                                                     unknowntree,
                                                     counter)

            splid_alignments.append(splid_alignment)
        else:
            pass  # print("No splids found.")

        print("\r", counter, end='')

    print("\nDone. Bye.")

    outtsv.close()
    AlignIO.write(splid_alignments, outfile_name + ".phy", "phylip-relaxed")
