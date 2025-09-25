in_fasta = config['in_fasta']
in_taxonomy = config['in_taxonomy']
out_dir = config['out_dir']
ref_pckg = config['ref_pckg']

rule all:
    input:
        f'{out_dir}/summary/summary_report.txt'

rule make_blast_db:
    input:
        in_fasta
    output:
        db = in_fasta + '.nin'
    log:
        f'{out_dir}/logs/make_blast_db.log'
    shell:
        '''
        makeblastdb -in {input} -dbtype nucl &> {log}
        '''

rule blast:
    input:
        db = in_fasta + '.nin',
        query = in_fasta
    output:
        f'{out_dir}/blast/blast_results.txt'
    params:
        outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/blast/blast.log'
    shell:
        '''
        blastn -db {input.db} -query {input.query} -out {output} -outfmt "{params.outfmt}" -num_threads {threads} 2> {log}
        '''

rule parse_blast:
    input:
        f'{out_dir}/blast/blast_results.txt'
    output:
        f'{out_dir}/parsed/parsed_results.txt'
    log:
        f'{out_dir}/logs/parse_blast/parse_blast.log'
    run:
        # check/correct orientation
        # reconcile the same hit from different queries
        with open(output[0], 'w') as f:
            f.write("# Placeholder for parsed blast results\n")

rule mafft_add:
    input:
        ref = f'{out_dir}/reference/reference{ref_pckg}.fasta',
        query = f'{out_dir}/parsed/parsed_results.txt'
    output:
        f'{out_dir}/alignment/aligned.fasta'
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/mafft/mafft.log'
    shell:
        '''
        mafft --add {input.query} --keeplength --thread {threads} {input.ref} > {output} 2> {log}
        '''

rule pplacer:
    input:
        f'{out_dir}/alignment/aligned.fasta'
    output:
        f'{out_dir}/placement/placement.jplace'
    params:
        refpkg = lambda wildcards: config.get('pplacer_refpkg', f'reference_packages/SSU_Amoebozoa_Metazoa{ref_pckg}.refpkg')
    threads:
        workflow.cores
    log:
        f'{out_dir}/logs/pplacer/pplacer.log'
    shell:
        '''
        pplacer -c {params.refpkg} -o {output} -j {threads} {input} &> {log}
        '''

rule guppy:
    input:
        f'{out_dir}/placement/placement.jplace'
    output:
        f'{out_dir}/guppy/guppy_results.txt'
    log:
        f'{out_dir}/logs/guppy/guppy.log'
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