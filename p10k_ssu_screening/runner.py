"""Snakemake runner for P10K SSU screening tool."""

import os
import sys
import subprocess
from pathlib import Path


class SnakemakeRunner:
    """Handles execution of Snakemake workflows for SSU screening."""
    
    def __init__(self, fasta, taxonomy, output_dir, threads=1, 
                 only_18s=False, only_16s=False, dry_run=False, 
                 keep_temp=False, verbose=False):
        """Initialize the runner with configuration parameters."""
        self.fasta = os.path.abspath(fasta)
        self.taxonomy = os.path.abspath(taxonomy)
        self.output_dir = os.path.abspath(output_dir)
        self.threads = threads
        self.only_18s = only_18s
        self.only_16s = only_16s
        self.dry_run = dry_run
        self.keep_temp = keep_temp
        self.verbose = verbose
        
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
    
    def make_config(self, eighteen_s=True):
        """
        Make config list to be passed to Snakemake.
        
        Args:
            eighteen_s (bool): Whether to use 18S configuration
            
        Returns:
            str: Space-separated config string for Snakemake
        """
        if eighteen_s:
            ref_pckg = ''
        else:
            ref_pckg = '_16S'

        config_items = [
            f'in_fasta={self.fasta}',
            f'in_taxonomy={self.taxonomy}',
            f'out_dir={self.output_dir}',
            f'ref_pckg={ref_pckg}',
        ]

        return ' '.join(config_items)
    
    def run_snakemake(self, eighteen_s=True):
        """
        Combine snakemake command fragments and run snakemake.
        
        Args:
            eighteen_s (bool): Whether to run 18S analysis
        """
        # Build Snakemake command
        smk_frags = [
            'snakemake',
            f'-s {self.snakefile_path}',
            f'--config {self.make_config(eighteen_s)}',
            f'--cores {self.threads}',
            '--rerun-incomplete',
            '--keep-going',
            '--nolock',
        ]
        
        # Add optional flags
        if self.dry_run:
            smk_frags.append('--dry-run')
        
        if self.verbose:
            smk_frags.append('--verbose')
        
        if not self.keep_temp:
            smk_frags.append('--delete-temp-output')
        
        smk_cmd = ' '.join(smk_frags)
        
        if self.verbose:
            print(f"Running: {smk_cmd}")
        
        try:
            subprocess.run(smk_cmd, shell=True, executable='/bin/bash', check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error: Command '{e.cmd}' failed with return code {e.returncode}.", file=sys.stderr)
            raise
    
    def run(self):
        """Run the complete screening workflow."""
        # Ensure output directory exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        
        if self.verbose:
            print(f"Starting P10K SSU screening...")
            print(f"Input FASTA: {self.fasta}")
            print(f"Taxonomy: {self.taxonomy}")
            print(f"Output directory: {self.output_dir}")
            print(f"Threads: {self.threads}")
        
        # Determine which analyses to run
        run_18s = not self.only_16s
        run_16s = not self.only_18s
        
        try:
            if run_18s:
                if self.verbose:
                    print("Running 18S SSU screening...")
                self.run_snakemake(eighteen_s=True)
            
            if run_16s:
                if self.verbose:
                    print("Running 16S SSU screening...")
                self.run_snakemake(eighteen_s=False)
                
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Snakemake execution failed: {e}")
        
        if self.verbose:
            print("Screening workflow completed successfully!")