"""
Combined workflow for SSU retrieval and phylogenetic placement.
This workflow can run both parts or skip retrieval if SSU sequences are provided.
"""

in_fasta = config['in_fasta']
out_dir = config['out_dir']
ref_pckg = config['ref_pckg']
query_fasta = config['query_fasta']
ref_aln = config['ref_aln']
trim_ref_aln = config['trim_ref_aln']
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
            f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.json'
        ]
    
    # Start building target list - always include BUSCO if running retrieval
    targets = []
    if run_retrieval:
        targets.append(f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.json')
    
    # For placement or full mode, check SSU sequences
    if not run_retrieval:
        # When skipping retrieval, check user-provided SSU file
        ssu_file = ssu_input or in_fasta
    else:
        # When running retrieval, check the checkpoint output
        ck = checkpoints.parse_mafft_trim.get()
        ssu_file = ck.output[0]

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
        targets.append(f'{out_dir}/summary/__NO_HITS__')
    else:
        # Check if placements were successful
        placement_status_file = f'{out_dir}/pplacer/placement_status.txt'
        
        # Trigger the checkpoint
        ck_placement = checkpoints.check_placements.get()
        placement_status_file = ck_placement.output[0]
        
        # Wait for file to be available
        wait_seconds = 10
        for _ in range(wait_seconds * 2):
            if os.path.exists(placement_status_file):
                break
            time.sleep(0.5)
        
        # Check placement status
        has_placements = False
        if os.path.exists(placement_status_file):
            with open(placement_status_file, 'r') as f:
                status = f.read().strip()
                has_placements = (status == 'HAS_PLACEMENTS')
        
        if not has_placements:
            # No valid placements - create marker
            targets.append(f'{out_dir}/summary/__NO_VALID_PLACEMENTS__')
        else:
            # Valid placements - run full pipeline
            targets.extend([
                f'{out_dir}/summary/taxonomy_summary.csv',
                f'{out_dir}/summary/sequence_classifications.csv',
                f'{out_dir}/summary/rank_counts.csv',
                f'{out_dir}/summary/placement_tree.pdf'
            ])
            # Only include BUSCO summary if we ran retrieval
            if run_retrieval:
                targets.append(f'{out_dir}/summary/busco_summary.json')
    
    return targets


# Create list of static targets that don't depend on checkpoint
static_targets = []
if run_retrieval:
    static_targets.append(f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.json')

rule all:
    input:
        static_targets,
        choose_targets

if run_retrieval:
    rule bacterial_busco:
        input:
            in_fasta
        output:
            f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.json'
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
    rule parse_blast:
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

rule mafft_trim:
    input:
        f'{out_dir}/parsed_blast/parsed_sequences_for_pplacer.fasta'
    output:
        f'{out_dir}/mafft_trim/parsed_sequences_for_pplacer.aln'
    params:
        ref = trim_ref_aln
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/mafft_trim.log'
    shell:
        '''
        mkdir -p {out_dir}/mafft_trim
        # Check if input file is empty or has no sequences
        if [ ! -s {input} ] || ! grep -q "^>" {input}; then
            echo "No sequences to align - input file is empty" > {log}
            touch {output}
        else
            mafft --add {input} --globalpair --keeplength --thread {threads} {params.ref} > {output} 2> {log}
        fi
        '''

checkpoint parse_mafft_trim:
    input:
        f'{out_dir}/mafft_trim/parsed_sequences_for_pplacer.aln',
        trim_ref_aln
    output:
        f'{out_dir}/mafft_trim/parsed_sequences_for_pplacer.fasta'
    log:
        f'{out_dir}/logs/parse_mafft_trim.log'
    run:
        from Bio import SeqIO
        ref_records = set()
        with open(input[1], 'r') as ref_file:
            for record in SeqIO.parse(ref_file, 'fasta'):
                ref_records.add(record.id)
        
        aln_rerecords = set()
        with open(input[0], 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                aln_rerecords.add(record.id)
        
        to_keep = aln_rerecords - ref_records
        with open(output[0], 'w') as out_f:
            for record in SeqIO.parse(input[0], 'fasta'):
                if record.id in to_keep and len(record.seq.replace('-', '')) >= 500:
                    SeqIO.write(record, out_f, 'fasta')


rule mafft:
    input:
        f'{out_dir}/mafft_trim/parsed_sequences_for_pplacer.fasta'
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
        mafft --add {input[0]} --keeplength --thread {threads} {params.ref} > {output} 2> {log}
        '''

rule pplacer:
    input:
        f'{out_dir}/mafft/aligned.fasta'
    output:
        f'{out_dir}/pplacer/placement.jplace'
    params:
        refpkg = ref_pckg,
        ref_aln = ref_aln
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/pplacer.log'
    shell:
        '''
        python3 {workflow.basedir}/../scripts/run_pplacer.py \
          -i {input} \
          -c {params.refpkg} \
          -o {output} \
          -r {params.ref_aln} \
          -j {threads} \
          -l {log} &> {log}
        '''

checkpoint check_placements:
    input:
        f'{out_dir}/pplacer/placement.jplace'
    output:
        f'{out_dir}/pplacer/placement_status.txt'
    log:
        f'{out_dir}/logs/check_placements.log'
    run:
        import json
        
        has_placements = False
        try:
            with open(input[0], 'r') as f:
                jplace_data = json.load(f)
                placements = jplace_data.get('placements', [])
                has_placements = len(placements) > 0
        except Exception as e:
            with open(log[0], 'w') as log_f:
                log_f.write(f"Error reading jplace file: {e}\n")
            has_placements = False
        
        with open(output[0], 'w') as f:
            if has_placements:
                f.write('HAS_PLACEMENTS\n')
            else:
                f.write('NO_PLACEMENTS\n')
        
        with open(log[0], 'a') as log_f:
            log_f.write(f"Status: {'HAS_PLACEMENTS' if has_placements else 'NO_PLACEMENTS'}\n")

rule rppr:
    input:
        jplace = f'{out_dir}/pplacer/placement.jplace',
        status = f'{out_dir}/pplacer/placement_status.txt'
    output:
        f'{out_dir}/rppr/results.db'
    params:
        refpkg = ref_pckg
    log:
        f'{out_dir}/logs/rppr.log'
    run:
        # Check if we have valid placements
        with open(input.status, 'r') as f:
            status = f.read().strip()
        
        if status != 'HAS_PLACEMENTS':
            # Skip rppr for empty placements
            with open(log[0], 'w') as log_f:
                log_f.write("No valid placements - skipping rppr\n")
            shell("touch {output}")
        else:
            shell("rppr prep_db -c {params.refpkg} --sqlite {output} {input.jplace} &> {log}")

rule guppy:
    input:
        jplace = f'{out_dir}/pplacer/placement.jplace',
        db = f'{out_dir}/rppr/results.db',
        status = f'{out_dir}/pplacer/placement_status.txt'
    output:
        f'{out_dir}/guppy/placement.tree'
    params:
        refpkg = ref_pckg
    log:
        f'{out_dir}/logs/guppy.log'
    run:
        # Check if we have valid placements
        with open(input.status, 'r') as f:
            status = f.read().strip()
        
        if status != 'HAS_PLACEMENTS':
            # Skip guppy for empty placements
            with open(log[0], 'w') as log_f:
                log_f.write("No valid placements - skipping guppy\n")
            shell("touch {output}")
        else:
            shell("""
            guppy classify -c {params.refpkg} --sqlite {input.db} {input.jplace} &> {log}
            guppy tog -o {output} {input.jplace} >> {log} 2>> {log}
            """)

rule summarize:
    input:
        rppr = f'{out_dir}/rppr/results.db',
        guppy = f'{out_dir}/guppy/placement.tree',
        parsed = ssu_input if not run_retrieval else f'{out_dir}/parsed_blast/parsed_sequences.fasta',
        busco = f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.json' if run_retrieval else [],
        status = f'{out_dir}/pplacer/placement_status.txt'
    output:
        taxonomy_summary = f'{out_dir}/summary/taxonomy_summary.csv',
        sequence_classifications = f'{out_dir}/summary/sequence_classifications.csv',
        rank_counts = f'{out_dir}/summary/rank_counts.csv',
        tree_pdf = f'{out_dir}/summary/placement_tree.pdf',
        busco_summary = f'{out_dir}/summary/busco_summary.json'
    log:
        f'{out_dir}/logs/summarize.log'
    params:
        supergroup = supergroup_of_interest,
        busco_file = f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.json'
    run:
        # Check if we have valid placements
        with open(input.status, 'r') as f:
            status = f.read().strip()
        
        if status != 'HAS_PLACEMENTS':
            # Should not reach here due to choose_targets, but handle gracefully
            with open(log[0], 'w') as log_f:
                log_f.write("No valid placements - cannot generate full summary\n")
            for out in [output.taxonomy_summary, output.sequence_classifications, 
                       output.rank_counts, output.tree_pdf, output.busco_summary]:
                shell(f"touch {out}")
        else:
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
        busco = f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.json' if run_retrieval else []
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

# Rule to create summary when no valid placements were generated
rule summarize_no_valid_placements:
    input:
        placement_status = f'{out_dir}/pplacer/placement_status.txt',
        busco = f'{out_dir}/busco/short_summary.specific.bacteria_odb12.busco.json' if run_retrieval else []
    output:
        marker = f'{out_dir}/summary/__NO_VALID_PLACEMENTS__'
    log:
        f'{out_dir}/logs/summarize.log'
    shell:
        '''
        mkdir -p {out_dir}/summary
        echo "No valid placements generated - all query sequences failed pplacer" > {log}
        echo "Check {out_dir}/mafft/aligned_filtered.fasta and {out_dir}/mafft/removed_sequences.txt for details" >> {log}
        touch {output.marker}
        '''
