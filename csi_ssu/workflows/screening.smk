in_fasta = config['in_fasta']
out_dir = config['out_dir']
ref_pckg = config['ref_pckg']
query_fasta = config['query_fasta']
ref_aln = config['ref_aln']
supergroup_of_interest = config['supergroup_of_interest']
pplacer_cutoff_length = config['pplacer_cutoff_length']

import os
import time

input_basename = os.path.splitext(os.path.basename(in_fasta))[0]

def choose_targets(wildcards=None):
    ck = checkpoints.parse_blast.get()
    parsed = ck.output.pplacer_fasta

    # Wait briefly for the checkpoint output to be available
    wait_seconds = 20
    for _ in range(wait_seconds * 2):
        if os.path.exists(parsed):
            break
        time.sleep(0.5)

    # If file does not exist or is empty / has no FASTA headers, treat as no-hits
    if not os.path.exists(parsed):
        no_hits = True
    else:
        try:
            size = os.path.getsize(parsed)
            if size == 0:
                no_hits = True
            else:
                # count fasta headers
                with open(parsed, 'r') as fh:
                    headers = sum(1 for line in fh if line.startswith('>'))
                no_hits = (headers == 0)
        except Exception:
            no_hits = True

    if no_hits:
        return [f'{out_dir}/summary/__NO_HITS__']
    else:
        summary_files = [
            f'{out_dir}/summary/taxonomy_summary.csv',
            f'{out_dir}/summary/sequence_classifications.csv',
            f'{out_dir}/summary/rank_counts.csv',
            f'{out_dir}/summary/placement_tree.pdf'
        ]
        return summary_files


rule all:
    input:
        choose_targets

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
        outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq qcovs qcovhsp sstrand',
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

checkpoint parse_blast:
    input:
        f'{out_dir}/blast/results.txt'
    output:
        txt = f'{out_dir}/parsed_blast/parsed_results.txt',
        fasta = f'{out_dir}/parsed_blast/parsed_sequences.fasta',
        pplacer_fasta = f'{out_dir}/parsed_blast/parsed_sequences_for_pplacer.fasta'
    params:
        pplacer_cutoff_length = pplacer_cutoff_length
    log:
        f'{out_dir}/logs/parse_blast/parse_blast.log'
    shell:
        '''
        mkdir -p {out_dir}/parsed_blast
        python3 {workflow.basedir}/../scripts/parse_blast.py {input} {in_fasta} {output.txt} {output.fasta} {output.pplacer_fasta} {params.pplacer_cutoff_length} &> {log}
        '''

rule mafft:
    input:
        pplacer_fasta = lambda wc: checkpoints.parse_blast.get().output.pplacer_fasta,
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
        mkdir -p {out_dir}/mafft
        mafft --add {input.pplacer_fasta} --keeplength --thread {threads} {params.ref} > {output} 2> {log}
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

rule rppr:
    input:
        f'{out_dir}/pplacer/placement.jplace'
    output:
        f'{out_dir}/rppr/results.db'
    params:
        refpkg = ref_pckg
    log:
        f'{out_dir}/logs/rppr.log'
    shell:
        '''
        rppr prep_db -c {params.refpkg} --sqlite {output} {input} &> {log}
        '''

rule guppy:
    input:
        f'{out_dir}/pplacer/placement.jplace',
        f'{out_dir}/rppr/results.db'
    output:
        f'{out_dir}/guppy/placement.tree'
    params:
        refpkg = ref_pckg
    log:
        f'{out_dir}/logs/guppy.log'
    shell:
        '''
        guppy classify -c {params.refpkg} --sqlite {input[1]} {input[0]} &> {log}
        guppy tog -o {output} {input[0]} &>> {log}
        '''

rule summarize:
    input:
        rppr = f'{out_dir}/rppr/results.db',
        guppy = f'{out_dir}/guppy/placement.tree',
        parsed = lambda wc: checkpoints.parse_blast.get().output.fasta
    output:
        taxonomy_summary = f'{out_dir}/summary/taxonomy_summary.csv',
        sequence_classifications = f'{out_dir}/summary/sequence_classifications.csv',
        rank_counts = f'{out_dir}/summary/rank_counts.csv',
        tree_pdf = f'{out_dir}/summary/placement_tree.pdf'
    log:
        f'{out_dir}/logs/summarize.log'
    params:
        supergroup = supergroup_of_interest
    shell:
        '''
        mkdir -p {out_dir}/summary
        python3 {workflow.basedir}/../scripts/summarize_results.py {input.rppr} {output.taxonomy_summary} {output.sequence_classifications} {output.rank_counts} {input_basename} &> {log}
        python3 {workflow.basedir}/../scripts/plot_tree.py {input.guppy} {output.tree_pdf} {input.rppr} {params.supergroup} "" {input.parsed} &>> {log}
        '''

# Rule to create summary files when no hits were found
rule summarize_no_hits:
    input:
        parsed = lambda wc: checkpoints.parse_blast.get().output.fasta
    output:
        marker = f'{out_dir}/summary/__NO_HITS__'
    log:
        f'{out_dir}/logs/summarize.log'
    shell:
        '''
        mkdir -p {out_dir}/summary
        echo "No hits found meeting coverage threshold" > {log}
        touch {output.marker}
        '''
