#!/usr/bin/env python3
"""Generates a database by trimming sequences in the input FASTA with
   Cutadapt using a specified primer pair and a desired number of
   mismatches. Only sequences with complete primer matches are reported.
   Ambiguous nucleotides may optionally be allowed in both the primers
   and the target region, the primers only, or in neither. A minimum and
   maximum sequence length (excluding primers) may also be specified.
   Taxonomic ranks present in the database may be specified by providing
   a highest and/or lowest rank to maintain, or by providing a list of
   desired ranks. A minimum assigned taxonomic rank for all sequences
   may also be enforced."""

import argparse
import os
import subprocess
import sys

from math import ceil, floor, log10
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def arg_fasta_file(in_fasta):
    """Checks if a specified FASTA exists and has a valid extension."""
    fasta_exts = ('.fa', '.fas', '.fasta', '.txt')
    if not os.path.isfile(in_fasta):
        raise argparse.ArgumentTypeError(
            '"{}" cannot be found.'.format(in_fasta)
        )
    if not any(in_fasta.endswith(f) for f in fasta_exts):
        raise argparse.ArgumentTypeError(
            '"{}" does not end with a valid extension of ".fa", ".fas", '
            '".fasta" or ".txt"'.format(in_fasta)
        )
    return in_fasta

def arg_primer(in_seq):
    """Checks that a primer contains only valid IUPAC DNA nucleotide
       codes and returns it as a Biopython sequence."""
    acc_nucs = {'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T',
                'V', 'W', 'Y'}
    if not all(n in acc_nucs for n in in_seq.upper()):
        raise argparse.ArgumentTypeError(
            '"{}" is not composed solely of IUPAC DNA nucleotide '
            'codes.'.format(in_seq)
        )
    return Seq(in_seq.upper())

def arg_pos_int(in_val):
    """Checks if a value is a positive integer."""
    if not in_val.isdigit() or int(in_val) == 0:
        raise argparse.ArgumentTypeError(
            '"{}" is not a positive integer.'.format(in_val)
        )
    return int(in_val)

def parse_args():
    """Handle user-supplied and default arguments for this script."""
    parser = argparse.ArgumentParser(
        description=(
            'Generates a database by trimming sequences in the input FASTA '
            'with Cutadapt using a specified primer pair and a desired number '
            'of mismatches. Only sequences with complete primer matches are '
            'reported. Ambiguous nucleotides may optionally be allowed in '
            'both the primers and the target region, the primers only, or in '
            'neither. A minimum and maximum sequence length (excluding '
            'primers) may also be specified. Taxonomic ranks present in the '
            'database may be specified by providing a highest and/or lowest '
            'rank to maintain, or by providing a list of desired ranks. A '
            'minimum assigned taxonomic rank for all sequences may also be '
            'enforced.'
        )
    )
    parser.add_argument(
        '-f', '--fasta', nargs='?', required=True,
        help='The input FASTA file to be trimmed.'
    )
    parser.add_argument(
        '-5', '--primer_5', nargs='?', type=arg_primer, required=True,
        help='The sequence of the 5\' primer.'
    )
    parser.add_argument(
        '-3', '--primer_3', nargs='?', type=arg_primer, required=True,
        help='The sequence of the 3\' primer.'
    )
    parser.add_argument(
        '-m', '--mismatches', nargs='?', type=arg_pos_int, required=True,
        help='The maximum number of mismatches to allow across both primers.'
    )
    parser.add_argument(
        '-p', '--pipeline', nargs='?', choices=('blast', 'dada2'),
        required=True,
        help='The pipeline the database FASTA should be formatted for.'
    )
    parser.add_argument(
        '-o', '--out', nargs='?', default=sys.stdout,
        help='The file where the database FASTA should be written. (Default '
             '"standard out")'
    )
    parser.add_argument(
        '--ambig', nargs='?', default='none', type=lambda a: a.lower(),
        choices=('all', 'none', 'primers'),
        help='Sequence portions where ambiguous nucleotides are allowed. '
             '(Default "none")'
    )
    parser.add_argument(
        '--mismatches_5', nargs='?', default=None, type=arg_pos_int,
        help='The maximum number of mismatches to allow for the 5\' primer. '
             '(Optional)'
    )
    parser.add_argument(
        '--mismatches_3', nargs='?', default=None, type=arg_pos_int,
        help='The maximum number of mismatches to allow for the 3\' primer. '
             '(Optional)'
    )
    parser.add_argument(
        '--min', nargs='?', default=0, type=arg_pos_int,
        help='The minimum sequence length to allow (excluding primers). '
             '(Optional)'
    )
    parser.add_argument(
        '--max', nargs='?', default=float('inf'), type=arg_pos_int,
        help='The maximum sequence length to allow (excluding primers). '
             '(Optional)'
    )
    parser.add_argument(
        '--assigned', nargs='?', default=None, type=lambda a: a.lower(),
        help='Sequences must have taxonomic assignments to at least this rank '
             'in order to be retained. (Optional, not case-sensitive)'
    )
    parser.add_argument(
        '--require', nargs='+', default=[], type=lambda a: a.lower(),
        help='One or more taxonomic levels (e.g. Viridiplantae) all sequences '
             'must contain to be retained. (Optional, not case-sensitive)'
    )
    parser.add_argument(
        '--exclude', nargs='+', default=[], type=lambda a: a.lower(),
        help='One or more taxonomic levels (e.g. Viridiplantae) that, if '
             'present, will cause a sequence to be excluded. (Optional, not '
             'case-sensitive)'
    )
    parser.add_argument(
        '--ranks', nargs='+', default=[], type=lambda a: a.lower(),
        help='One or more taxonomic ranks that, if specified, will be the '
             'only ranks retained. (Optional, not case-sensitive)'
    )
    parser.add_argument(
        '--highest', nargs='?', default=None, type=lambda a: a.lower(),
        help='The highest taxonomic rank to include. (Optional)'
    )
    parser.add_argument(
        '--lowest', nargs='?', default=None, type=lambda a: a.lower(),
        help='The lowest taxonomic rank to include. (Optional)'
    )
    parser.add_argument(
        '--primers', action='store_true',
        help='Keep primer sequences in the output FASTA. (Optional)'
    )
    args = parser.parse_args()
    args.primer_3 = args.primer_3.reverse_complement()
    if args.mismatches_5 == None:
        args.mismatches_5 = args.mismatches
    if args.mismatches_3 == None:
        args.mismatches_3 = args.mismatches
    return args

def fasta_to_pos_dict(args):
    """Returns a dictionary of every record description in the input
       FASTA."""
    fasta_pos_dict = {}
    for rec_idx, rec in enumerate(SeqIO.parse(args.fasta, 'fasta')):
        fasta_pos_dict[rec.description] = rec_idx+1
    return fasta_pos_dict

def calc_bin_size(value):
    """Determines an appropriate report bin size based on the value
       provided."""
    return max(10**(floor(log10(value))-1), 1)

def calc_error_rate(primer_mis, primer_seq):
    """Converts the number of allowed mismatches into an error rate for
       Cutadapt based on the length of the primer sequence."""
    return ceil(100*primer_mis/len(primer_seq))/100

def run_cutadapt(args):
    """Runs Cutadapt on the input FASTA using the designated primers and
       mismatches."""
    p5_err = calc_error_rate(args.mismatches_5, args.primer_5)
    p3_err = calc_error_rate(args.mismatches_3, args.primer_3)
    cmd = (
        'cutadapt -g "{};e={}" -a "{};e={}" --no-indels -n 2 '
        '--action=lowercase --discard-untrimmed --quiet {}'.format(
            args.primer_5, p5_err, args.primer_3, p3_err, args.fasta
        )
    )
    cutadapt_subproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True, shell=True
    )
    return cutadapt_subproc

def format_seq(rec, p5_len, p3_len):
    """Trims a FASTA sequence to just the region of interest and its
       accompanying primers, and reports the length."""
    target_range = [i for (i, n) in enumerate(rec.seq) if n.isupper()]
    start = 0 if target_range[0]-p5_len < 0 else target_range[0]-p5_len
    end = target_range[-1]+p3_len+1
    target_seq = rec.seq[start:end]
    return target_seq, len(target_range)

def calc_primer_mis(match_dict, primer_seq, fasta_seq):
    """Determines the number of mismatches between the primer and FASTA
       sequences."""
    mis = 0
    for primer_nuc, fasta_nuc in zip(primer_seq, fasta_seq):
        if fasta_nuc not in match_dict[primer_nuc]:
            mis += 1
    return mis

def calc_fasta_mis(args, match_dict, seq):
    """Validates that the FASTA sequence meets the mismatch requirements
       specified by the user and returns a tuple of mismatch counts if
       it does."""
    p5_fasta_seq = seq[:len(args.primer_5)]
    p3_fasta_seq = seq[-len(args.primer_3):]
    p5_mis = calc_primer_mis(match_dict, args.primer_5, p5_fasta_seq)
    p3_mis = calc_primer_mis(match_dict, args.primer_3, p3_fasta_seq)
    total_mis = p5_mis+p3_mis
    if (p5_mis <= args.mismatches_5 and p3_mis <= args.mismatches_3 and
        total_mis <= args.mismatches):
        fasta_mis = (total_mis, p5_mis, p3_mis)
    else:
        fasta_mis = ()
    return fasta_mis

def fasta_to_desc_mis_dicts(args, fasta_pos_dict):
    """Loads trimmed FASTA records into a dictionary with the record
       description as the key, and creates a dictionary containing the
       mismatch count for each accession."""
    unambig_nucs = {'A', 'C', 'G', 'T'}
    match_dict = {
      'A': ('A'),
      'C': ('C'),
      'G': ('G'),
      'T': ('T'),
      'K': ('G', 'T', 'K'),
      'M': ('A', 'C', 'M'),
      'R': ('A', 'G', 'R'),
      'S': ('C', 'G', 'S'),
      'W': ('A', 'T', 'W'),
      'Y': ('C', 'T', 'Y'),
      'B': ('C', 'G', 'T', 'K', 'S', 'Y', 'B'),
      'D': ('A', 'G', 'T', 'K', 'R', 'W', 'D'),
      'H': ('A', 'C', 'T', 'M', 'W', 'Y', 'H'),
      'V': ('A', 'C', 'G', 'M', 'R', 'S', 'V'),
      'N': ('A', 'C', 'G', 'T', 'K', 'M', 'R', 'S', 'W', 'Y', 'B', 'D', 'H',
            'V', 'N')
    }
    p5_len = len(args.primer_5)
    p3_len = len(args.primer_3)
    pos_total = 2*len(fasta_pos_dict)
    bin_report_size = calc_bin_size(pos_total)
    fasta_desc_dict = {}
    fasta_mis_dict = {}
    last_rep_bin = 0
    sys.stderr.write('Trimming {:,} sequences\n'.format(pos_total))
    for rec in SeqIO.parse(run_cutadapt(args).stdout, 'fasta'):
        rank_names = [r.split(':')[1].lower() for r in rec.id.split(';')]
        seq_nucs = set(rec.seq)
        if (
                (any(r in args.require for r in rank_names) or
                 not args.require) and
                (all(r not in args.exclude for r in rank_names) or
                 not args.exclude) and
                (not set(n.lower() for n in seq_nucs).issuperset(seq_nucs))
        ):
            rec.seq, target_len = format_seq(rec, p5_len, p3_len)
            if args.ambig == 'all':
                validation_nucs = set()
            elif args.ambig == 'none':
                validation_nucs = set(n.upper() for n in set(rec.seq))
            else:
                validation_nucs = set(n for n in set(rec.seq) if n.isupper())
            if (rec.seq and p5_len+target_len+p3_len == len(rec.seq) and
                    validation_nucs.issubset(unambig_nucs)):
                fasta_mis = calc_fasta_mis(args, match_dict, rec.seq.upper())
                if (args.min <= target_len <= args.max and fasta_mis):
                    fasta_mis_dict[rec.description.split(' ')[1]] = fasta_mis
                    if not args.primers:
                        rec.seq = Seq(''.join(n for n in rec.seq if
                                              n.isupper()))
                    if rec.description not in fasta_desc_dict:
                        fasta_desc_dict[rec.description] = rec
                    else:
                        fasta_desc_dict[rec.description] = None
                else:
                    fasta_desc_dict[rec.description] = None
        pos_bin = int(fasta_pos_dict[rec.description]/bin_report_size)
        if pos_bin > last_rep_bin:
            last_rep_bin = pos_bin
            pos_pct = '{:.2%}'.format(bin_report_size*pos_bin/pos_total)
            sys.stderr.write(
                'Trimmed {:,} sequences ({})\n'.format(bin_report_size*pos_bin,
                                                       pos_pct)
            )
    sys.stderr.write('Trimmed all sequences\n')
    return fasta_desc_dict, fasta_mis_dict

def format_subrank(parent_rank, number):
    """Formats the name for ambiguous subranks based on their parent and
       number."""
    return '{}_Subrank_{}'.format(parent_rank, number)

def parse_taxonomy(rank_dict, rank_order, lineage):
    """Parses taxonomy from a list and updates the order for how
       taxonomic ranks should be recorded."""
    ambig_ranks = ('clade', 'no_rank')
    prev_rank = ['root', 0]
    for rank, rank_name in lineage:
        if rank_name[:12] != 'unclassified':
            if rank in ambig_ranks:
                prev_rank[1] += 1
                rank = format_subrank(*prev_rank)
                if prev_rank[0] not in rank_dict:
                    rank_order.append(prev_rank[0])
                    rank_dict[prev_rank[0]] = prev_rank[1]
                elif prev_rank[1] > rank_dict[prev_rank[0]]:
                    rank_dict[prev_rank[0]] = prev_rank[1]
            else:
                if rank not in rank_dict:
                    if not rank_order:
                        rank_order.append(rank)
                    else:
                        rank_idx = max(i for i, r in enumerate(rank_order) if
                                       prev_rank[0] in r)
                        rank_order = (rank_order[:rank_idx+1] + [rank] +
                                      rank_order[rank_idx+1:])
                    rank_dict[rank] = 0
                prev_rank = [rank, 0]
    return rank_dict, rank_order

def ids_to_rank_list(args, fasta_desc_dict):
    """Returns an ordered list of ranks from the trimmed FASTA IDs."""
    rank_dict = {}
    rank_order = []
    for desc in sorted([d for d in fasta_desc_dict if fasta_desc_dict[d]],
                       key=lambda r: r.split(' ')[1]):             
        lineage = [r.split(':') for r in fasta_desc_dict[desc].id.split(';')]
        rank_dict, rank_order = parse_taxonomy(rank_dict, rank_order, lineage)
    if args.ranks:
        rank_list = [r for r in rank_order if r.lower() in args.ranks]
    else:
        rank_list = []
        for rank in rank_order:
            if rank != 'root':
                rank_list.append(rank)
            for rank_num in range(rank_dict[rank]):
                rank_list.append(format_subrank(rank, rank_num+1))
    if args.highest and args.highest in [r.lower() for r in rank_list]:
        high_idx = [r.lower() for r in rank_list].index(args.highest)
        rank_list = rank_list[high_idx:]
    if args.lowest and args.lowest in [r.lower() for r in rank_list]:
        low_idx = [r.lower() for r in rank_list].index(args.lowest)
        rank_list = rank_list[:low_idx+1]
    return rank_list

def comp_tax(rank_list, lineage):
    """Fills out a complete taxonomy for a sequence's lineage."""
    tax = []
    lineage_dict = {}
    prev_rank = ['root', 0]
    for rank, rank_name in lineage:
        if rank == 'no_rank':
            prev_rank[1] += 1
            rank = format_subrank(*prev_rank)
            lineage_dict[rank] = rank_name
        else:
            prev_rank = [rank, 0]
            lineage_dict[rank] = rank_name
    rank_indices = [rank_list.index(r) for r in lineage_dict if r in rank_list]
    lowest_rank = max(rank_indices) if rank_indices else -1
    for rank_idx, rank in enumerate(rank_list):
        if rank in lineage_dict:
            name = lineage_dict[rank]
        else:
            if rank_idx > lowest_rank:
                name = 'Unknown'
            else:
                name = 'Unclassified'
        tax.append((rank, name))
    return tax

def desc_to_seq_dict(args, fasta_desc_dict, rank_list):
    """Loads Cutadapt-trimmed FASTA records from the description
       dictionary and stores them in a new dictionary by their sequence.
       This enables accessions with the same taxonomy and sequence to
       be merged to avoid overrepresentation in the database FASTA."""
    fasta_seq_dict = {}
    for desc in sorted([d for d in fasta_desc_dict if fasta_desc_dict[d]],
                       key=lambda r: r.split(' ')[1]):
        lineage = [r.split(':') for r in fasta_desc_dict[desc].id.split(';')]
        tax = comp_tax(rank_list, lineage)
        if (
            (args.assigned in [r[0] for r in tax] and
             (args.assigned, 'Unknown') not in tax) or
            not args.assigned
        ):
            tax = ';'.join(':'.join((t[0].capitalize(), t[1])) for t in tax)
            acc = desc.split(' ')[1]
            seq = str(fasta_desc_dict[desc].seq)
            if seq not in fasta_seq_dict:
                fasta_seq_dict[seq] = {}
            if tax not in fasta_seq_dict[seq]:
                fasta_seq_dict[seq][tax] = set()
            fasta_seq_dict[seq][tax].add(acc)
    return fasta_seq_dict

def mis_counts_to_str(fasta_mis_set):
    """Generates a condensed mismatch count string for a FASTA
       description."""
    fasta_mis_strs = []
    for fasta_mis in sorted(fasta_mis_set):
        fasta_mis_strs.append(','.join(str(f) for f in fasta_mis))
    fasta_mis_str = ';'.join(fasta_mis_strs)
    return fasta_mis_str

def seq_mis_dicts_to_db_fasta(args, fasta_seq_dict, fasta_mis_dict):
    """Writes records in the sequence dictionary to a database FASTA.
       Entries with the same sequence and taxonomy are merged, while
       entries with different sequences but the same taxonomy are
       flagged with a '|' and a number (provided the DADA2 pipeline
       option is not specified)."""
    seq_recs = []
    seq_count_dict = {}
    for seq in sorted(fasta_seq_dict):
        for tax in sorted(fasta_seq_dict[seq]):
            accs = ';'.join(sorted(fasta_seq_dict[seq][tax]))
            mmcs = mis_counts_to_str(set(fasta_mis_dict[a] for a in
                                         fasta_seq_dict[seq][tax]))
            desc = '|'.join([accs, mmcs])
            seq_recs.append([tax.split(';'), desc, seq])
            if tax not in seq_count_dict:
                seq_count_dict[tax] = 1
            else:
                seq_count_dict[tax] += 1
    seq_mult_dict = {}
    for tax in seq_count_dict:
        if seq_count_dict[tax] > 1:
            seq_mult_dict[tax] = 0
    out_handle = open(args.out, 'w') if args.out != sys.stdout else sys.stdout
    for seq_rec in sorted(seq_recs):
        tax, desc, seq = seq_rec
        if args.pipeline == 'blast':
            tax = ';'.join(tax)
            if tax in seq_mult_dict:
                seq_mult_dict[tax] += 1
                tax = '{}|{}'.format(tax, seq_mult_dict[tax])
            rec = SeqRecord(Seq(seq), id=tax, description=desc)
            SeqIO.write(rec, out_handle, 'fasta')
        else:
            tax = '{};'.format(';'.join([t.split(':')[1] for t in tax if
                                         t.split(':')[1]!= 'Unknown']))
            rec = SeqRecord(Seq(seq), id=tax, description='')
            SeqIO.write(rec, out_handle, 'fasta-2line')
    if out_handle != sys.stdout:
        out_handle.close()

def main():
    """The main sequence of the script."""
    args = parse_args()
    fasta_pos_dict = fasta_to_pos_dict(args)
    fasta_desc_dict, fasta_mis_dict = fasta_to_desc_mis_dicts(args,
                                                              fasta_pos_dict)
    rank_list = ids_to_rank_list(args, fasta_desc_dict)
    fasta_seq_dict = desc_to_seq_dict(args, fasta_desc_dict, rank_list)
    seq_mis_dicts_to_db_fasta(args, fasta_seq_dict, fasta_mis_dict)

if __name__ == '__main__':
    main()
