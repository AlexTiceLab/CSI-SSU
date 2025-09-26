in_fasta = config['in_fasta']
in_taxonomy = config['in_taxonomy']
out_dir = config['out_dir']
ref_pckg = config['ref_pckg']
query_fasta = config['query_fasta']
ref_aln = config['ref_aln']

# Extract basename of input fasta file for database naming
import os
input_basename = os.path.splitext(os.path.basename(in_fasta))[0]


rule all:
    input:
        f'{out_dir}/summary/summary_report.txt'

rule make_blast_db:
    input:
        in_fasta
    output:
        db = f'{out_dir}/blast_db/{input_basename}.nin',
        other_files = multiext(f'{out_dir}/blast_db/{input_basename}', '.nhr', '.nsq')
    params:
        db_name = f'{out_dir}/blast_db/{input_basename}'
    log:
        f'{out_dir}/logs/make_blast_db.log'
    shell:
        '''
        mkdir -p {out_dir}/blast_db
        makeblastdb -in {input} -dbtype nucl -out {params.db_name} &> {log}
        '''

rule blast:
    input:
        db = f'{out_dir}/blast_db/{input_basename}.nin',
        query = query_fasta
    output:
        f'{out_dir}/blast/results.txt'
    params:
        outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq',
        db_name = f'{out_dir}/blast_db/{input_basename}'
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/blast.log'
    shell:
        '''
        mkdir -p {out_dir}/blast
        blastn -db {params.db_name} -query {input.query} -out {output} -outfmt "{params.outfmt}" -num_threads {threads} 2> {log}
        '''

rule parse_blast:
    input:
        f'{out_dir}/blast/results.txt'
    output:
        txt = f'{out_dir}/parsed_blast/parsed_results.txt',
        fasta = f'{out_dir}/parsed_blast/parsed_sequences.fasta'
    log:
        f'{out_dir}/logs/parse_blast/parse_blast.log'
    run:
        # Dictionary to store the longest sequence for each unique header
        sequences = {}
        blast_lines = {}
        
        with open(input[0], 'r') as infile:
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) < 13:
                    continue
                qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sseq = parts
                
                # Keep track of the longest sequence for each unique sseqid
                if sseqid not in sequences or len(sseq) > len(sequences[sseqid]):
                    sequences[sseqid] = sseq
                    blast_lines[sseqid] = line
        
        # Write output files
        with open(output.txt, 'w') as out_txt, open(output.fasta, 'w') as out_fasta:
            for sseqid, sequence in sequences.items():
                out_txt.write(blast_lines[sseqid])
                out_fasta.write(f'>{sseqid}\n{sequence}\n')

rule cdhit:
    input:
        f'{out_dir}/parsed_blast/parsed_sequences.fasta'
    output:
        f'{out_dir}/cd-hit/clustered_sequences.fasta'
    params:
        threshold = 0.99
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/cd-hit.log'
    shell:
        '''
        cd-hit -i {input} -o {output} -c {params.threshold} -T {threads} &> {log}
        '''

rule mafft:
    input:
        f'{out_dir}/cd-hit/clustered_sequences.fasta'
    output:
        f'{out_dir}/mafft/aligned.fasta'
    params:
        ref = ref_aln
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/mafft.log'
    shell:
        '''
        mafft --add {input} --keeplength --thread {threads} {params.ref} > {output} 2> {log}
        '''

rule pplacer:
    input:
        f'{out_dir}/mafft/aligned.fasta'
    output:
        f'{out_dir}/pplacer/placement.jplace'
    params:
        refpkg = ref_pckg
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/pplacer.log'
    shell:
        '''
        pplacer -c {params.refpkg} -o {output} -j {threads} {input} &> {log}
        '''

rule guppy:
    input:
        f'{out_dir}/pplacer/placement.jplace'
    output:
        f'{out_dir}/guppy/guppy_results.txt'
    log:
        f'{out_dir}/logs/guppy.log'
    shell:
        '''
        guppy classify {input} > {output} 2> {log}
        ''' 

rule summarize:
    input:
        f'{out_dir}/guppy/guppy_results.txt'
    output:
        f'{out_dir}/summary/summary_report.txt'
    log:
        f'{out_dir}/logs/summarize.log'
    shell:
        '''
        cat {input} > {output} 2> {log}
        '''
