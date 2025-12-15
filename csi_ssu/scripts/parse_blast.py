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

    # Store strand orientation per HSP and HSP ranges per sseqid
    # Each entry is a dict: {'start': int, 'end': int, 'is_reverse': bool}
    s_hsps = {}

    # Maximum allowed extracted sequence length (bp)
    MAX_EXTRACT_LEN = 2500

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

            # Normalize HSP coordinates to forward orientation
            hsp_start = min(sstart, send)
            hsp_end = max(sstart, send)
            if sseqid not in s_hsps:
                s_hsps[sseqid] = []
            s_hsps[sseqid].append({
                'start': hsp_start,
                'end': hsp_end,
                'is_reverse': is_reverse,
            })

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
        # Ensure coordinates are within the sequence bounds
        seq_len = len(full_seq)
        if start < 1:
            start = 1
        if end > seq_len:
            end = seq_len

        # Determine orientation consistency between the HSP that defines
        # the merged `start` and the HSP that defines the merged `end`.
        # If they come from different HSPs with differing orientations,
        # skip this subject.
        region_is_reverse = False
        hsp_list = s_hsps.get(sseqid, [])
        if not hsp_list:
            continue

        # Find HSPs that contain the merged start/end (or whose boundary equals them)
        start_hsp = next((h for h in hsp_list if h['start'] <= start <= h['end']), None)
        end_hsp = next((h for h in hsp_list if h['start'] <= end <= h['end']), None)

        if start_hsp is None or end_hsp is None:
            # As a fallback, try matching exact boundary values
            start_hsp = next((h for h in hsp_list if h['start'] == start), start_hsp)
            end_hsp = next((h for h in hsp_list if h['end'] == end), end_hsp)

        if start_hsp is None or end_hsp is None:
            # Can't determine contributing HSPs reliably; skip
            continue

        if start_hsp is end_hsp or start_hsp['is_reverse'] == end_hsp['is_reverse']:
            region_is_reverse = start_hsp['is_reverse']
        else:
            # Start and end come from HSPs with different orientations: skip
            continue

        # If extraction is too long, try stepping `end` down to the next
        # BLAST HSP `send` coordinates (descending). If none of the HSP
        # ends produce a short enough region, exclude the region.
        if (end - start + 1) > MAX_EXTRACT_LEN:
            candidate_ends = sorted(set(h['end'] for h in hsp_list))
            # consider only ends smaller than current merged end
            candidates = [e for e in candidate_ends if e < end]
            trimmed = False
            for c in sorted(candidates, reverse=True):
                if (c - start + 1) <= MAX_EXTRACT_LEN:
                    end = c
                    trimmed = True
                    break

            if not trimmed:
                # no HSP end shrank it enough; exclude this region
                continue

        # If trimming removed the entire region, skip
        if end <= start:
            continue

        extracted = full_seq[start-1:end]

        # Reverse complement according to the region orientation
        if region_is_reverse:
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
