#!/usr/bin/env python3
"""Downloads sequence and taxonomy information from the NCBI database
   and generates a FASTA file containing those sequences, their reverse
   complement, the full taxonomy of those sequences (set as the FASTA
   record ID) and the accession of each sequence (set as the FASTA
   description)."""

import argparse
import os
import subprocess
import sys
import xml.etree.ElementTree as ElementTree

from math import floor, log10
from time import sleep
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
            'Downloads sequence and taxonomy information from the NCBI '
            'database and generates a FASTA file containing those sequences, '
            'their reverse complement, the full taxonomy of those sequences '
            '(set as the FASTA record ID) and the accession of each sequence '
            '(set as the FASTA description).'
        )
    )
    parser.add_argument(
        '-q', '--query', nargs='?', required=True,
        help='The query to be submitted to NCBI.'
    )
    parser.add_argument(
        '-o', '--out', nargs='?', default=sys.stdout,
        help='The FASTA file to write sequence and taxonomy data to. (Default '
             '"standard out")'
    )
    parser.add_argument(
        '--records', nargs='?', default=1000, type=arg_pos_int,
        help='The maximum number of records to download at one time. This '
             'value may be lowered if the script runs out of memory. Values '
             'well under 100 or over 1000 are generally not recommended. '
             '(Optional, default 1000)'
    )
    args = parser.parse_args()
    return args

def run_eutils_uids(args):
    """Runs the NCBI E-utilities to download the unique identifier
       numbers for all sequences returned by the query."""
    esearch_cmd = 'esearch -db nucleotide -query "{}"'.format(args.query)
    efetch_cmd = 'efetch -format uid'
    eutil_subproc = subprocess.Popen(
        ' | '.join([esearch_cmd, efetch_cmd]), stdout=subprocess.PIPE,
        universal_newlines=True, shell=True
    )
    return eutil_subproc

def run_efetch_ids(id_type, id_slice):
    """Runs the NCBI EFetch utility to download sequences or taxonomy
       information in the XML format."""
    if id_type == 'sequences':
        cmd = ('efetch -db nucleotide -id "{}" -format fasta '
               '-mode xml'.format(','.join(str(i) for i in id_slice)))
    else:
        cmd = ('efetch -db taxonomy -id "{}" '
               '-format xml'.format(','.join(str(i) for i in id_slice)))
    efetch_subproc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        universal_newlines=True, shell=True
    )
    return efetch_subproc

def xml_rec_to_fasta(xml_rec):
    """Converts individual XML records into FASTA sequences, using the
       taxid as the sequence ID and the accession as the description."""
    acc = xml_rec.find('TSeq_accver').text
    taxid = xml_rec.find('TSeq_taxid').text
    seq = xml_rec.find('TSeq_sequence').text
    return SeqRecord(Seq(seq), id=taxid, description=acc)

def format_name(rank_name):
    """Reformats a name to avoid reserved characters."""
    rep_tups = ((' ', '_'), (':', '-'), (';', ','), ('|', '_'))
    for rep_tup in rep_tups:
        orig_char, rep_char = rep_tup
        rank_name = rank_name.replace(orig_char, rep_char)
    return rank_name

def xml_rec_to_taxid_and_lineage(xml_rec):
    """Returns the taxid and lineage for an XML record. The lineage is
       reported as a semicolon-delimited string, with each rank and its
       assignment stored as a colon-delimited string."""
    taxids = [xml_rec.find('TaxId').text]
    if xml_rec.find('AkaTaxIds'):
        for taxid in xml_rec.find('AkaTaxIds'):
            taxids.append(taxid.text)
    org_name = format_name(xml_rec.find('ScientificName').text)
    lineage_ex = xml_rec.find('LineageEx')
    lineage_list = []
    for taxon in lineage_ex:
        rank_name = format_name(taxon.find('ScientificName').text)
        rank = format_name(taxon.find('Rank').text)
        lineage_list.append((rank, rank_name))
    if (not any(r[0] == 'species' for r in lineage_list) and
            any(r[0] == 'genus' for r in lineage_list)):
        lineage_list.append(('species', org_name))
    lineage = ';'.join(':'.join(r) for r in lineage_list)
    return taxids, lineage

def parse_xml_rec(id_type, return_object, temp_handle, xml_rec):
    """Parses XML entries downloaded by the NCBI based on whether they
       are sequence or taxonomic records."""
    if id_type == 'sequences':
        fasta_rec = xml_rec_to_fasta(xml_rec)
        return_object.add(int(fasta_rec.id))
        SeqIO.write(fasta_rec, temp_handle, 'fasta')
    else:
        taxids, lineage = xml_rec_to_taxid_and_lineage(xml_rec)
        for taxid in taxids:
            return_object[taxid] = lineage
    return return_object

def download_and_format_xml_data(id_type, return_object, temp_handle,
                                 id_slice):
    """Downloads data in the XML format from NCBI, and retries the
       download after 60 seconds if errors are encountered. XML records
       are processed according to their ID type if the download is
       successful."""
    try:
        efetch_subproc = run_efetch_ids(id_type, id_slice)
        for xml_rec in ElementTree.fromstringlist(efetch_subproc.stdout):
            return_object = parse_xml_rec(id_type, return_object, temp_handle,
                                          xml_rec)
        return return_object
    except ElementTree.ParseError:
        sys.stderr.write('Warning: Network issue encountered, waiting 60 '
                         'seconds to resume download\n')
        sleep(60)
        download_and_format_xml_data(id_type, return_object, temp_handle,
                                     id_slice)

def calc_bin_size(value):
    """Determines an appropriate report bin size based on the value
       provided."""
    return max(10**(floor(log10(value))-1), 1)

def download_ncbi_data(args, id_set, id_type):
    """Downloads information for ID numbers in blocks of a specified
       size. Sequence records are converted to the FASTA format and
       written to a temporary file, and the taxids are returned in a
       set; taxonomy information is returned in a dictionary."""
    return_object = set() if id_type == 'sequences' else {}
    if id_type == 'sequences':
        if args.out == sys.stdout:
            temp_handle = open('temp.fasta', 'w')
        else:
            temp_handle = open('{}.temp'.format(args.out), 'w')
    else:
        temp_handle = None
    id_list = sorted(id_set)
    id_total = len(id_list)
    id_report_size = calc_bin_size(id_total)
    last_id_bin = 0
    sys.stderr.write('Downloading {:,} {}\n'.format(id_total, id_type))
    for id_idx_start in range(0, id_total, min(id_total, args.records)):
        id_idx_end = id_idx_start+min(id_total, args.records)
        download_and_format_xml_data(id_type, return_object, temp_handle,
                                     id_list[id_idx_start:id_idx_end])
        id_bin = int(id_idx_end/id_report_size)
        if id_bin > last_id_bin:
            last_id_bin = id_bin
            id_pct = '{:.2%}'.format(id_idx_end/id_total)
            sys.stderr.write(
                'Downloaded {:,} {} ({})\n'.format(id_idx_end, id_type, id_pct)
            )
    sys.stderr.write('Downloaded all {}\n'.format(id_type))
    if temp_handle:
        temp_handle.close()
    return return_object

def ncbi_data_to_fasta(args):
    """Downloads sequence and taxonomy data from NCBI and combines them
       into a single FASTA. A temporary FASTA file, which is later
       removed, is used to store sequences while downloading. Each
       sequence, as well as its reverse complement, is included in the
       final file."""
    eutil_uids = set(int(l.rstrip()) for l in run_eutils_uids(args).stdout)
    eutil_taxids = download_ncbi_data(args, eutil_uids, 'sequences')
    tax_dict = download_ncbi_data(args, eutil_taxids, 'taxonomies')
    out_handle = open(args.out, 'w') if args.out != sys.stdout else sys.stdout
    if args.out == sys.stdout:
        temp_fasta = 'temp.fasta'
    else:
        temp_fasta = '{}.temp'.format(args.out)
    rec_total = len(eutil_uids)
    rec_report_size = calc_bin_size(rec_total)
    sys.stderr.write('Finalizing {:,} records\n'.format(rec_total))
    for rec_idx, rec in enumerate(SeqIO.parse(temp_fasta, 'fasta')):
        rec_idx += 1
        rec.description = rec.description.split(' ')[1]
        if rec.id in tax_dict:
            rec.id = tax_dict[rec.id]
            SeqIO.write(rec, out_handle, 'fasta')
            rec.seq = rec.seq.reverse_complement()
            SeqIO.write(rec, out_handle, 'fasta')
        else:
            sys.stderr.write('Warning: Skipping accession {} due to invalid '
                             'taxonomy data\n'.format(rec.description))
        if not rec_idx % rec_report_size:
            rec_pct = '{:.2%}'.format(rec_idx/rec_total)
            sys.stderr.write('Finalized {:,} records ({})\n'.format(rec_idx,
                                                                    rec_pct))
    sys.stderr.write('Finalized all records\n')
    if out_handle != sys.stdout:
        out_handle.close()
    os.remove(temp_fasta)

def main():
    """The main sequence of the script."""
    args = parse_args()
    ncbi_data_to_fasta(args)

if __name__ == '__main__':
    main()
