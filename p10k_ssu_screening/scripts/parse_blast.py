#!/usr/bin/env python3
"""
Parse BLAST results and extract sequences.
"""

import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO

def parse_blast_results(blast_results_file, input_fasta_file, output_txt_file, output_fasta_file, output_fasta_file_for_pplacer, cutoff_length_for_pplacer):
    """
    Parse BLAST results and extract sequences.
    
    Args:
        blast_results_file: Path to BLAST results file
        input_fasta_file: Path to input FASTA file (for reference)
        output_txt_file: Path to output parsed results text file
        output_fasta_file: Path to output parsed sequences FASTA file
    """
    # Dictionary to store the longest sequence for each unique header
    sequences = {}
    blast_lines = {}
    sstart_end_dict = {}
    # Record strand/orientation for subject sequences (True if reverse)
    sstrand_dict = {}
    
    with open(blast_results_file, 'r') as infile:
        for line in infile:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sseq, qcovs, qcovhsp, sstrand = line.strip().split('\t')
            
            qcovs = float(qcovs)
            if qcovs > 30:
                # Use the subject sequence directly from BLAST results
                qstart = int(qstart)
                qend = int(qend)

                # Determine strand/orientation. BLAST's sstrand can be 'plus'/'minus' or similar.
                strand_flag = False
                if sstrand:
                    ss = sstrand.strip().lower()
                    if ss in ('minus', '-', 'rev', 'reverse', 'r'):
                        strand_flag = True

                # Also consider qstart>qend as an indicator of reverse orientation
                if int(qstart) > int(qend):
                    strand_flag = True

                # Record orientation for this subject sequence (if any record indicates reverse,
                # keep it True).
                sstrand_dict[sseqid] = sstrand_dict.get(sseqid, False) or strand_flag

                # Use the subject sequence directly from BLAST results; apply reverse complement
                # if the alignment indicates reverse orientation.
                if sstrand_dict.get(sseqid, False):
                    sequence = str(Seq(sseq).reverse_complement())
                    print(f'RC: {sseqid} length {len(sequence)}')
                else:
                    sequence = sseq
                    print(f'Not RC: {sseqid} length {len(sequence)}')

                # Keep track of the longest sequence for each unique sseqid
                if sseqid not in sequences or len(sequence) > len(sequences[sseqid]):
                    sequences[sseqid] = sequence
                    blast_lines[sseqid] = line
            
            if sseqid in sstart_end_dict:
                if int(sstart) < sstart_end_dict[sseqid][0]:
                    sstart_end_dict[sseqid] = (int(sstart), sstart_end_dict[sseqid][1])
                if int(send) > sstart_end_dict[sseqid][1]:
                    sstart_end_dict[sseqid] = Seq((sstart_end_dict[sseqid][0], int(send))).reverse_complement
            else:
                sstart_end_dict[sseqid] = (int(sstart), int(send))
    
    with open(input_fasta_file, 'r') as fasta_file:
        records = SeqIO.parse(fasta_file, 'fasta')
        for record in records:
            if record.id in sstart_end_dict and record.id in blast_lines:
                start, end = sstart_end_dict[record.id]
                # Check if coordinates indicate reverse orientation (start > end)
                needs_rc = start > end
                if needs_rc:
                    start, end = end, start
                
                extracted_seq = record.seq[start-1:end]  # Adjust for 0-based indexing
                # Apply reverse complement if needed based on coordinates or stored strand info
                if needs_rc or sstrand_dict.get(record.id, False):
                    try:
                        extracted_seq = extracted_seq.reverse_complement()
                    except Exception:
                        # If extracted_seq is a string fallback to Seq wrapper
                        extracted_seq = Seq(str(extracted_seq)).reverse_complement()

                # Use empty string length if no BLAST-derived sequence exists yet
                current_len = len(sequences.get(record.id, ''))
                if len(extracted_seq) > current_len:
                    sequences[record.id] = str(extracted_seq)

    # Ensure output directories exist
    os.makedirs(os.path.dirname(output_txt_file), exist_ok=True)
    os.makedirs(os.path.dirname(output_fasta_file), exist_ok=True)
    
    # Check if any sequences passed the coverage threshold
    if not sequences:
        print('Warning: No sequences found with coverage > 30%. Output files will be empty.')
        with open(output_txt_file, 'w') as out_txt:
            out_txt.write('')
        with open(output_fasta_file, 'w') as out_fasta:
            out_fasta.write('')
        with open(output_fasta_file_for_pplacer, 'w') as out_fasta:
            out_fasta.write('')
        return
    
    # Write output files
    with open(output_txt_file, 'w') as out_txt, open(output_fasta_file, 'w') as out_fasta, open(output_fasta_file_for_pplacer, 'w') as out_fasta_pplacer:
        for sseqid, sequence in sequences.items():
            out_txt.write(blast_lines[sseqid])
            out_fasta.write(f'>{sseqid}\n{sequence}\n')
            if len(sequence) >= cutoff_length_for_pplacer:
                out_fasta_pplacer.write(f'>{sseqid}\n{sequence}\n')
    
    print(f"Successfully parsed {len(sequences)} sequences from BLAST results.")

if __name__ == '__main__':
    blast_results_file = sys.argv[1]
    input_fasta_file = sys.argv[2]
    output_txt_file = sys.argv[3]
    output_fasta_file = sys.argv[4]
    output_fasta_file_for_pplacer = sys.argv[5]
    cutoff_length_for_pplacer = int(sys.argv[6])
    
    parse_blast_results(blast_results_file, input_fasta_file, output_txt_file, output_fasta_file, output_fasta_file_for_pplacer, cutoff_length_for_pplacer)