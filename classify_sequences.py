#!/usr/bin/env python3
"""Assigns taxonomy to sequences in a DADA2 CSV by aligning them to a
   database FASTA using BLAST."""

import argparse
import os
import subprocess
import sys

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

def arg_perc(in_val):
    """Checks if a value is valid percent number, and returns it as a
       value between 0 and 1."""
    if not in_val.replace('.', '').isdigit() or float(in_val) > 100:
        raise argparse.ArgumentTypeError(
            '"{}" is not a valid percent number.'.format(in_val)
        )
    return float(in_val)/100

def arg_pos_float(in_val):
    """Checks if a value is a positive number."""
    if not in_val.replace('.', '').isdigit():
        raise argparse.ArgumentTypeError(
            '"{}" is not a positive number.'.format(in_val)
        )
    return float(in_val)

def arg_pos_int(in_val):
    """Checks if a value is a positive integer."""
    if not in_val.isdigit():
        raise argparse.ArgumentTypeError(
            '"{}" is not a positive integer.'.format(in_val)
        )
    return int(in_val)

def parse_args():
    """Handle user-supplied and default arguments for this script."""
    parser = argparse.ArgumentParser(
        description=(
            'Assigns taxonomy to sequences in a DADA2 CSV by aligning them to '
            'a database FASTA using BLAST.'
        )
    )
    parser.add_argument(
        '-c', '--csv', nargs='?', required=True,
        help='The input CSV to use.'
    )
    parser.add_argument(
        '-d', '--database', nargs='?', required=True,
        help='The database FASTA to use.'
    )
    parser.add_argument(
        '-o', '--out', nargs='?', default=sys.stdout,
        help='The CSV file to write results to. (Default "standard out")'
    )
    parser.add_argument(
        '--taxa', nargs='?',
        help='A file of taxa to preferentially retain when multiple best hits '
             'are present. (Optional)'
    )
    parser.add_argument(
        '--consensus', nargs='?', default=0.9, type=arg_perc,
        help='The minimum percent consensus for taxonomy to be assigned. '
             '(Optional, default 90)'
    )
    parser.add_argument(
        '--length', nargs='?', default=0.9, type=arg_perc,
        help='When considering the length of a database sequence and the '
             'alignment length reported for a BLAST hit, the smaller of these '
             'two values must be at least this percentage of the larger '
             'value. (Optional, default 90)'
    )
    parser.add_argument(
        '--identity', nargs='?', default=0.9, type=arg_perc,
        help='The minimum percent identity for sequence hits to be '
             'considered. (Optional, default 90)'
    )
    parser.add_argument(
        '--e_value', nargs='?', default=1.0, type=arg_pos_float,
        help='The minimum E-value for sequence hits to be considered. '
             '(Optional, default 1.0)'
    )
    parser.add_argument(
        '--threads', nargs='?', default=1, type=arg_pos_int,
        help='The number of threads to use for BLAST searches. (Optional, '
             'default 1)'
    )
    args = parser.parse_args()
    return args

def build_database(args):
    """Builds an NCBI database if it does not already exist."""
    db_exts = ('nhr', 'nin', 'nsq')
    temp_files = []
    if not all(os.path.isfile('{}.{}'.format(args.database, d)) for d in
               db_exts):
        cmd = 'makeblastdb -in {} -dbtype nucl'.format(args.database)
        subprocess.call(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            universal_newlines=True, shell=True
        )
        for db_ext in db_exts:
            temp_files.append('{}.{}'.format(args.database, db_ext))
    return temp_files

def dada2_to_fasta(args):
    """Generates a temporary FASTA file from the input DADA2 CSV."""
    with open(args.csv, 'r') as in_handle:
        for line_idx, line in enumerate(in_handle):
            if line_idx == 0:
                seqs = [s for s in line.rstrip().replace('"', '').split(',') if
                        s]
    seq_id_width = len(str(len(seqs)))
    temp_fasta = '{}.temp'.format(args.csv)
    with open(temp_fasta, 'w') as temp_handle:
        for seq_idx, seq in enumerate(seqs):
            seq_id = '{:0>{width}}'.format(seq_idx+1, width=seq_id_width)
            seq_rec = SeqRecord(Seq(seq), id=seq_id, description='')
            SeqIO.write(seq_rec, temp_handle, 'fasta')
    return temp_fasta

def taxa_to_set(args):
    """Reads in a colon-delimited file of rank and taxa to preferentially
       retain."""
    taxa_set = set()
    if args.taxa:
        with open(args.taxa, 'r') as in_handle:
            for line in in_handle:
                taxa_set.add(line.rstrip())
    return taxa_set

def run_blastn(args):
    """Aligns sequences against the BLAST database using blastn."""
    cmd = (
        'blastn -task "blastn-short" -db {} -query {}.temp -evalue {} '
        '-word_size 4 -perc_identity {} -outfmt "6 qseqid sseqid pident '
        'length evalue" -num_threads {}'.format(
            args.database, args.csv, args.e_value, 100*args.identity,
            args.threads
        )
    )
    blastn_subproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr = subprocess.DEVNULL,
        universal_newlines=True, shell=True
    )
    return blastn_subproc

def blast_to_dict(args, db_fasta_dict, taxa_set):
    """Loads the best BLAST hits into a dictionary."""
    blast_dict = {}
    for line in run_blastn(args).stdout:
        line = line.rstrip().split('\t')
        seq_id = line[0]
        seq_hit = line[1]
        # BLAST drops the last character if it is a '.', so add it
        # back in.
        if seq_hit not in db_fasta_dict:
            seq_hit = '{}.'.format(seq_hit)
        if any(t in taxa_set for t in seq_hit.split(';')):
            taxa_score = 1
        else:
            taxa_score = 0
        ident = float(line[2])
        q_len = int(line[3])
        # Use the negative of the E-value so greater than
        # comparisons work as expected for all values.
        e_val = -1*float(line[4])
        r_len = len(db_fasta_dict[seq_hit])
        cov = min(q_len/r_len, r_len/q_len)
        vals = [e_val, ident, cov]
        seq_hit = seq_hit.split('|')[0]
        if args.length <= cov <= 1:
            if seq_id not in blast_dict:
                blast_dict[seq_id] = vals + [taxa_score, set([seq_hit])]
            elif all(v >= p for v, p in zip(vals, blast_dict[seq_id])):
                if any(v > p for v, p in zip(vals, blast_dict[seq_id])):
                    blast_dict[seq_id] = vals + [taxa_score, set([seq_hit])]
                elif taxa_score > blast_dict[seq_id][3]:
                    blast_dict[seq_id] = vals + [taxa_score, set([seq_hit])]
                elif taxa_score == blast_dict[seq_id][3]:
                    blast_dict[seq_id][4].add(seq_hit)
    return blast_dict

def get_cons_tax(args, taxa, cutoff_val):
    """Returns the shared consensus of taxa information, as well as the
       proportion of taxa in agreement at the last rank for which a
       consensus can be determined."""
    cons_tax = []
    cons_pro = 0
    taxa = [t[::-1] for t in taxa]
    for idx in range(len(taxa[0])):
        tax_count_dict = {}
        for tax in taxa:
            tax_count_dict[tax[idx]] = tax_count_dict.get(tax[idx], 0) + 1
        max_count_taxa = []
        max_count = 0
        for tax in tax_count_dict:
            if tax_count_dict[tax] > max_count:
                max_count_taxa = [tax]
                max_count = tax_count_dict[tax]
            elif tax_count_dict[tax] == max_count:
                max_count_taxa.append(tax)
        if (len(max_count_taxa) == 1 and
            tax_count_dict[max_count_taxa[0]] >= cutoff_val):
            if (max_count_taxa[0] in ('Unclassified', 'Unknown') and
                not cons_pro):
                cons_tax.append('NA')
            else:
                cons_tax.append('"{}"'.format(max_count_taxa[0]))
                if not cons_pro:
                    cons_pro = tax_count_dict[max_count_taxa[0]]/len(taxa)
        else:
            cons_tax.append('NA')
    cons_tax = cons_tax[::-1]
    return cons_tax, cons_pro

def parse_ranks_taxa(ranks_taxa):
    """Parses ranks and taxa names from a semicolon- and colon-delimited
       rank and taxa string."""
    ranks_taxa = [r.split(';') for r in ranks_taxa]
    ranks = [r.split(':')[0] for r in ranks_taxa[0]]
    taxa = []
    for rank_tax in ranks_taxa:
        taxa.append([r.split(':')[1] for r in rank_tax])
    return ranks, taxa

def write_dada2_table(args, dada2_fasta_dict, blast_dict):
    """Writes a DADA2 CSV table."""
    rank_len = None
    lines = []
    out_handle = open(args.out, 'w') if args.out != sys.stdout else sys.stdout
    for seq_id in sorted(dada2_fasta_dict):
        if seq_id in blast_dict:
            e_val, ident, cov, taxa_score, ranks_taxa = blast_dict[seq_id]
            ranks, taxa = parse_ranks_taxa(ranks_taxa)
            if not rank_len:
                header = (['Sequence'] + ranks +
                          ['Consensus', 'Coverage', 'Identity', 'E-value'])
                rank_len = len(ranks) + 4
                if args.taxa:
                    header.append('Filtered')
                    rank_len += 1
            cons_tax, cons_pro = get_cons_tax(args, taxa,
                                              args.consensus*len(taxa))
            line = '"{}",{},{},{},{},{}'.format(
                    str(dada2_fasta_dict[seq_id].seq), ','.join(cons_tax),
                    100*cons_pro, 100*cov, ident, -1*e_val
            )
            if args.taxa:
                filtered = 'Yes' if taxa_score == 1 else 'No'
                line = '{},{}'.format(line, filtered)
            lines.append(line)
        else:
            lines.append(((str(dada2_fasta_dict[seq_id].seq)), None))
    out_handle.write('{}\n'.format(','.join(header)))
    for line in lines:
        if isinstance(line, tuple):
            out_handle.write('"{}",{}\n'.format(line[0],
                                                ','.join(rank_len*['NA'])))
        else:
            out_handle.write('{}\n'.format(line))
    if out_handle != sys.stdout:
        out_handle.close()

def remove_temp_files(temp_files):
    """Cleans up temporary files produced by this script."""
    for temp_file in temp_files:
        os.remove(temp_file)

def main():
    """The main sequence of the script."""
    args = parse_args()
    temp_files = build_database(args)
    temp_fasta = dada2_to_fasta(args)
    db_fasta_dict = SeqIO.to_dict(SeqIO.parse(args.database, 'fasta'))
    dada2_fasta_dict = SeqIO.to_dict(SeqIO.parse(temp_fasta, 'fasta'))
    taxa_set = taxa_to_set(args)
    blast_dict = blast_to_dict(args, db_fasta_dict, taxa_set)
    write_dada2_table(args, dada2_fasta_dict, blast_dict)
    remove_temp_files(temp_files + [temp_fasta])

if __name__ == '__main__':
    main()
