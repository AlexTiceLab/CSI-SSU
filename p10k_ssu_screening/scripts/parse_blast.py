#!/usr/bin/env python3
"""
Parse BLAST results and extract full subject sequences,
with correct orientation based on BLAST subject coordinates.
"""

import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO


def parse_blast_results(blast_results_file, input_fasta_file,
                        output_txt_file, output_fasta_file,
                        output_fasta_file_for_pplacer, cutoff_length_for_pplacer):

    # Store merged subject start/end per sseqid
    sregions = {}

    # Store the BLAST line for each sseqid (longest HSP kept)
    blast_lines = {}

    # Store extracted sequences (longest kept)
    sequences = {}

    # Store strand orientation for each sseqid
    # True = reverse complement needed
    sstrand_dict = {}

    # --------------------------
    # 1. FIRST PASS: Parse BLAST
    # --------------------------
    with open(blast_results_file) as infile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) < 16:
                continue

            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, \
            sstart, send, evalue, bitscore, sseq, qcovs, qcovhsp, sstrand = parts

            qcovs = float(qcovs)
            if qcovs < 30:
                continue

            sstart = int(sstart)
            send = int(send)

            # Determine strand: BLAST rule → reverse if sstart > send
            is_reverse = sstart > send
            sstrand_dict[sseqid] = is_reverse

            # Merge HSP regions on subject coordinates
            if sseqid not in sregions:
                sregions[sseqid] = [min(sstart, send), max(sstart, send)]
            else:
                sregions[sseqid][0] = min(sregions[sseqid][0], sstart, send)
                sregions[sseqid][1] = max(sregions[sseqid][1], sstart, send)

            # Keep longest BLAST HSP line for output_txt_file
            if sseqid not in blast_lines:
                blast_lines[sseqid] = line
            else:
                # replace if this HSP is longer
                old_len = int(blast_lines[sseqid].split("\t")[3])
                if int(length) > old_len:
                    blast_lines[sseqid] = line

    # -------------------------------------
    # 2. SECOND PASS: Extract from FASTA
    # -------------------------------------
    fasta_records = SeqIO.to_dict(SeqIO.parse(input_fasta_file, "fasta"))

    for sseqid, (start, end) in sregions.items():

        if sseqid not in fasta_records:
            continue  # Subject not found in FASTA

        full_seq = fasta_records[sseqid].seq
        extracted = full_seq[start-1:end]

        # Reverse complement if BLAST indicates minus strand
        if sstrand_dict.get(sseqid, False):
            extracted = extracted.reverse_complement()

        # Keep the longest extraction per id
        if sseqid not in sequences or len(extracted) > len(sequences[sseqid]):
            sequences[sseqid] = str(extracted)

    # --------------------------
    # 3. WRITE OUTPUT FILES
    # --------------------------

    os.makedirs(os.path.dirname(output_txt_file), exist_ok=True)
    os.makedirs(os.path.dirname(output_fasta_file), exist_ok=True)

    # If no sequences survive, produce empty files
    if not sequences:
        open(output_txt_file, "w").close()
        open(output_fasta_file, "w").close()
        open(output_fasta_file_for_pplacer, "w").close()
        print("Warning: no sequences with qcovs > 30 found.")
        return

    with open(output_txt_file, "w") as out_txt, \
         open(output_fasta_file, "w") as out_fasta, \
         open(output_fasta_file_for_pplacer, "w") as out_pplacer:

        for sseqid, seq in sequences.items():
            out_txt.write(blast_lines.get(sseqid, ""))

            out_fasta.write(f">{sseqid}\n{seq}\n")

            if len(seq) >= cutoff_length_for_pplacer:
                out_pplacer.write(f">{sseqid}\n{seq}\n")

    print(f"Successfully extracted {len(sequences)} oriented sequences.")


if __name__ == "__main__":
    blast_results_file = sys.argv[1]
    input_fasta_file = sys.argv[2]
    output_txt_file = sys.argv[3]
    output_fasta_file = sys.argv[4]
    output_fasta_file_for_pplacer = sys.argv[5]
    cutoff_length_for_pplacer = int(sys.argv[6])

    parse_blast_results(blast_results_file, input_fasta_file,
                        output_txt_file, output_fasta_file,
                        output_fasta_file_for_pplacer,
                        cutoff_length_for_pplacer)
