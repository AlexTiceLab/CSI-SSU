"""
Combined workflow for SSU retrieval and phylogenetic placement.
This workflow can run both parts or skip retrieval if SSU sequences are provided.
"""

in_fasta = config['in_fasta']
out_dir = config['out_dir']
ref_pckg = config['ref_pckg']
query_fasta = config['query_fasta']
ref_aln = config['ref_aln']
supergroup_of_interest = config['supergroup_of_interest']
pplacer_cutoff_length = config['pplacer_cutoff_length']
busco_mode = config['busco_mode']
busco_downloads = config['busco_downloads']
skip_retrieval = config.get('skip_retrieval', False)
ssu_input = config.get('ssu_input', None)
workflow_mode = config.get('workflow_mode', 'full')

import os
import time

input_basename = os.path.splitext(os.path.basename(in_fasta))[0]

# Determine the SSU fasta file to use and whether to run each part
if workflow_mode == 'placement' or skip_retrieval:
    ssu_fasta_for_placement = ssu_input or in_fasta
    run_retrieval = False
    run_placement = True
else:
    ssu_fasta_for_placement = f'{out_dir}/parsed_blast/parsed_sequences_for_pplacer.fasta'
    run_retrieval = True
    run_placement = workflow_mode == 'full'

def choose_targets(wildcards=None):
    """Determine target files based on workflow mode."""
    # For retrieval-only mode, return retrieval outputs
    if workflow_mode == 'retrieval':
        return [
            f'{out_dir}/parsed_blast/parsed_sequences_for_pplacer.fasta',
            f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.txt'
        ]
    
    # Start building target list - always include BUSCO if running retrieval
    targets = []
    if run_retrieval:
        targets.append(f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.txt')
    
    # For placement or full mode, check SSU sequences
    if not run_retrieval:
        # When skipping retrieval, check user-provided SSU file
        ssu_file = ssu_input or in_fasta
    else:
        # When running retrieval, check the checkpoint output
        ck = checkpoints.parse_blast.get()
        ssu_file = ck.output.pplacer_fasta

        # Wait briefly for the checkpoint output to be available
        wait_seconds = 20
        for _ in range(wait_seconds * 2):
            if os.path.exists(ssu_file):
                break
            time.sleep(0.5)

    # Check if file exists and has content
    if not os.path.exists(ssu_file):
        no_hits = True
    else:
        try:
            size = os.path.getsize(ssu_file)
            if size == 0:
                no_hits = True
            else:
                # count fasta headers
                with open(ssu_file, 'r') as fh:
                    headers = sum(1 for line in fh if line.startswith('>'))
                no_hits = (headers == 0)
        except Exception:
            no_hits = True

    if no_hits:
        # When no hits, just create marker file
        # BUSCO results are already available in busco/ directory if run_retrieval was True
        targets.append(f'{out_dir}/summary/__NO_HITS__')
    else:
        targets.extend([
            f'{out_dir}/summary/taxonomy_summary.csv',
            f'{out_dir}/summary/sequence_classifications.csv',
            f'{out_dir}/summary/rank_counts.csv',
            f'{out_dir}/summary/placement_tree.pdf'
        ])
        # Only include BUSCO summary if we ran retrieval
        if run_retrieval:
            targets.append(f'{out_dir}/summary/busco_summary.txt')
    
    return targets


# Create list of static targets that don't depend on checkpoint
static_targets = []
if run_retrieval:
    static_targets.append(f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.txt')

rule all:
    input:
        static_targets,
        choose_targets

if run_retrieval:
    rule bacterial_busco:
        input:
            in_fasta
        output:
            f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.txt'
        params:
            lineage=f'{busco_downloads}/lineages/bacteria_odb12',
            mode=busco_mode,
            out_dir=f'{out_dir}/busco'
        threads:
            workflow.cores
        log:
            f'{out_dir}/logs/busco.log'
        shell:
            '''
            busco -i {input} -o {params.out_dir} -l {params.lineage} -m {params.mode} --cpu {threads} -f --offline &> {log}
            '''

if run_retrieval:
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

if run_retrieval:
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

if run_retrieval:
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
            f'{out_dir}/logs/parse_blast.log'
        shell:
            '''
            mkdir -p {out_dir}/parsed_blast
            python3 {workflow.basedir}/../scripts/parse_blast.py {input} {in_fasta} {output.txt} {output.fasta} {output.pplacer_fasta} {params.pplacer_cutoff_length} &> {log}
            '''

rule mafft:
    input:
        pplacer_fasta = ssu_fasta_for_placement if not run_retrieval else lambda wc: checkpoints.parse_blast.get().output.pplacer_fasta,
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
        jplace = f'{out_dir}/pplacer/placement.jplace',
        db = f'{out_dir}/rppr/results.db'
    output:
        f'{out_dir}/guppy/placement.tree'
    params:
        refpkg = ref_pckg
    log:
        f'{out_dir}/logs/guppy.log'
    shell:
        '''
        guppy classify -c {params.refpkg} --sqlite {input.db} {input.jplace} &> {log}
        guppy tog -o {output} {input.jplace} >> {log} 2>> {log}
        '''

rule summarize:
    input:
        rppr = f'{out_dir}/rppr/results.db',
        guppy = f'{out_dir}/guppy/placement.tree',
        parsed = ssu_input if not run_retrieval else lambda wc: checkpoints.parse_blast.get().output.fasta,
        busco = f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.txt' if run_retrieval else []
    output:
        taxonomy_summary = f'{out_dir}/summary/taxonomy_summary.csv',
        sequence_classifications = f'{out_dir}/summary/sequence_classifications.csv',
        rank_counts = f'{out_dir}/summary/rank_counts.csv',
        tree_pdf = f'{out_dir}/summary/placement_tree.pdf',
        busco_summary = f'{out_dir}/summary/busco_summary.txt'
    log:
        f'{out_dir}/logs/summarize.log'
    params:
        supergroup = supergroup_of_interest,
        busco_file = f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.txt'
    run:
        shell('''
        mkdir -p {out_dir}/summary
        python3 {workflow.basedir}/../scripts/summarize_results.py {input.rppr} {output.taxonomy_summary} {output.sequence_classifications} {output.rank_counts} {input_basename} &> {log}
        python3 {workflow.basedir}/../scripts/plot_tree.py {input.guppy} {output.tree_pdf} {input.rppr} {params.supergroup} "" {input.parsed} >> {log} 2>> {log}
        ''')
        if run_retrieval:
            shell('cp {params.busco_file} {output.busco_summary} >> {log} 2>> {log}')
        else:
            shell('echo "BUSCO not run - skipped retrieval" > {output.busco_summary}')

# Rule to create summary files when no hits were found
rule summarize_no_hits:
    input:
        busco = f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.txt' if run_retrieval else []
    output:
        marker = f'{out_dir}/summary/__NO_HITS__'
    log:
        f'{out_dir}/logs/summarize.log'
    shell:
        '''
        mkdir -p {out_dir}/summary
        echo "No SSU sequences found or provided" > {log}
        touch {output.marker}
        '''
