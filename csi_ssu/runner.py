"""Snakemake runner for CSI-SSU screening tool."""

import os
import sys
import subprocess
from pathlib import Path


class SnakemakeRunner:
    """Handles execution of Snakemake workflows for SSU screening."""
    
    def __init__(self, fasta, output_dir, supergroup, threads, 
                 dry_run, pplacer_cutoff_length, busco_mode, 
                 workflow_mode='full', skip_retrieval=False):
        """Initialize the runner with configuration parameters."""
        self.fasta = fasta
        self.output_dir = output_dir
        self.supergroup = supergroup
        self.threads = threads
        self.dry_run = dry_run
        self.pplacer_cutoff_length = pplacer_cutoff_length
        self.busco_mode = busco_mode
        self.workflow_mode = workflow_mode
        self.skip_retrieval = skip_retrieval
        
        # Auto-detect if SSU sequences already exist and skip retrieval if so
        self._auto_detect_existing_sequences()
        
        # Get the Snakefile path using a simple file system approach
        self.snakefile_path = self._get_snakefile_path()
    
    def _auto_detect_existing_sequences(self):
        """
        Auto-detect if SSU sequences have already been extracted.
        If found, automatically use them and skip retrieval to save resources.
        """
        # Only check if we're in full mode and not explicitly told to skip
        if self.workflow_mode == 'full' and not self.skip_retrieval:
            # Check for existing parsed sequences in the output directory
            parsed_fasta = Path(self.output_dir) / 'parsed_blast' / 'parsed_sequences_for_pplacer.fasta'
            
            if parsed_fasta.exists():
                # Check if file has content (at least one sequence)
                try:
                    with open(parsed_fasta, 'r') as f:
                        has_sequences = any(line.startswith('>') for line in f)
                    
                    if has_sequences:
                        print(f"Note: Found existing SSU sequences at {parsed_fasta}")
                        print("      Skipping retrieval step to save computational resources.")
                        print("      Using existing sequences for phylogenetic placement.")
                        self.skip_retrieval = True
                        self.existing_ssu_file = str(parsed_fasta.resolve())
                except Exception:
                    # If we can't read the file, don't skip retrieval
                    pass
    
    def _get_snakefile_path(self):
        """Get the path to the Snakemake workflow file."""
        # Use screening.smk for all modes - it handles all workflows
        current_dir = Path(__file__).parent
        workflow_file = 'screening.smk'
            
        snakefile_path = current_dir / 'workflows' / workflow_file

        if snakefile_path.exists():
            return str(snakefile_path.resolve())
        else:
            raise FileNotFoundError(f"Could not find Snakemake workflow file at: {snakefile_path}")
    
    def make_config(self):
        """
        Make config list to be passed to Snakemake.
            
        Returns:
            str: Space-separated config string for Snakemake
        """
        current_dir = Path(__file__).parent

        # Use the required supergroup parameter
        supergroup_of_interest = self.supergroup

        ref_pckg = f'{current_dir}/data/reference_packages/{supergroup_of_interest}.refpkg'
        query_file = f'{current_dir}/data/queries/pr2_version_5.1.0_18S_divisions_query.fasta'
        ref_aln = f'{ref_pckg}/{supergroup_of_interest}.aln'
        trim_ref_aln = f'{current_dir}/data/reference_alignments/ref_aln.fasta'
        busco_downloads = f'{current_dir}/data/busco_downloads'

        # Base config items
        config_items = [
            f'out_dir={self.output_dir}',
            f'ref_pckg={ref_pckg}',
            f'ref_aln={ref_aln}',
            f'trim_ref_aln={trim_ref_aln}',
            f'supergroup_of_interest={supergroup_of_interest}',
            f'pplacer_cutoff_length={self.pplacer_cutoff_length}',
            f'busco_mode={self.busco_mode}',
            f'busco_downloads={busco_downloads}',
            f'in_fasta={self.fasta}',
            f'query_fasta={query_file}',
            f'workflow_mode={self.workflow_mode}',
        ]
        
        # Add mode-specific config
        if self.workflow_mode == 'placement' or self.skip_retrieval:
            # For placement mode or when skipping retrieval
            ssu_file = getattr(self, 'existing_ssu_file', self.fasta)
            config_items.extend([
                f'skip_retrieval=True',
                f'ssu_input={ssu_file}',
            ])
        else:
            config_items.append(f'skip_retrieval=False')

        return ' '.join(config_items)
    
    def run_snakemake(self):
        """
        Combine snakemake command fragments and run snakemake.

        """
        # Build Snakemake command
        smk_frags = [
            'snakemake',
            f'-s {self.snakefile_path}',
            f'--config {self.make_config()}',
            f'--cores {self.threads}',
            '--rerun-incomplete',
            '--keep-going',
            '--nolock',
        ]
        
        # Add optional flags
        if self.dry_run:
            smk_frags.append('--dry-run')
        
        smk_cmd = ' '.join(smk_frags)
        
        try:
            subprocess.run(smk_cmd, shell=True, executable='/bin/bash', check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error: Command '{e.cmd}' failed with return code {e.returncode}.", file=sys.stderr)
            raise
    
    def run(self):
        """Run the complete screening workflow."""
        # Ensure output directory exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        
        # Run 18S analysis
        self.run_snakemake()
