"""Snakemake runner for P10K SSU screening tool."""

import os
import sys
import subprocess
from pathlib import Path


class SnakemakeRunner:
    """Handles execution of Snakemake workflows for SSU screening."""
    
    def __init__(self, fasta, output_dir, supergroup, threads=1, 
                 dry_run=False, pplacer_cutoff_length=200):
        """Initialize the runner with configuration parameters."""
        self.fasta = os.path.abspath(fasta)
        self.output_dir = os.path.abspath(output_dir)
        self.supergroup = supergroup
        self.threads = threads
        self.dry_run = dry_run
        self.pplacer_cutoff_length = pplacer_cutoff_length
        
        # Get the Snakefile path using a simple file system approach
        self.snakefile_path = self._get_snakefile_path()
    
    def _get_snakefile_path(self):
        """Get the path to the Snakemake workflow file."""
        # Look for the Snakefile relative to this module
        current_dir = Path(__file__).parent
        snakefile_path = current_dir / 'workflows' / 'screening.smk'

        if snakefile_path.exists():
            return str(snakefile_path.resolve())
        else:
            raise FileNotFoundError(f"Could not find Snakemake workflow file at: {snakefile_path}")
    
    def make_config(self):
        """
        Make config list to be passed to Snakemake.
        
        Args:
            eighteen_s (bool): Whether to use 18S configuration
            
        Returns:
            str: Space-separated config string for Snakemake
        """
        current_dir = Path(__file__).parent

        # Use the required supergroup parameter
        supergroup_of_interest = self.supergroup

        ref_pckg = current_dir / 'data' / 'reference_packages' / f'{supergroup_of_interest}.refpkg'
        query_file = current_dir / 'data' / 'queries' / 'pr2_version_5.1.0_18S_divisions_query.fasta'
        ref_aln = ref_pckg / f'{supergroup_of_interest}.aln'

        config_items = [
            f'in_fasta={self.fasta}',
            f'out_dir={self.output_dir}',
            f'ref_pckg={ref_pckg}',
            f'query_fasta={query_file}',
            f'ref_aln={ref_aln}',
            f'supergroup_of_interest={supergroup_of_interest}',
            f'pplacer_cutoff_length={self.pplacer_cutoff_length}'
        ]

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
