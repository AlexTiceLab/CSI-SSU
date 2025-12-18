#!/usr/bin/env python3
"""Command-line interface for CSI-SSU screening tool."""

import argparse
import sys
import os
from pathlib import Path
from .runner import SnakemakeRunner
from . import __version__

def create_parser():
    """Create and configure argument parser."""
    parser = argparse.ArgumentParser(
        prog='csi-ssu',
        description='CSI-SSU screening tool for genomic and transcriptomic data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full pipeline (SSU retrieval + phylogenetic placement)
  csi-ssu genome.fasta --supergroup Amoebozoa --data-type genome -o results/ -t 8
  
  # SSU retrieval only
  csi-ssu genome.fasta --supergroup Amoebozoa --data-type genome --mode retrieval -o results/ -t 8
  
  # Phylogenetic placement only (using pre-extracted SSU sequences)
  csi-ssu ssu_sequences.fasta --supergroup Amoebozoa --mode placement -o results/ -t 8
        """
    )
    
    # Version
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    
    # Required arguments
    parser.add_argument('fasta', help='Input FASTA file path (genome/transcriptome or pre-extracted SSU sequences)')

    parser.add_argument('--supergroup', 
                       type=str,
                       required=True,
                       help='Supergroup of interest (e.g., Amoebozoa, Excavata, TSAR)')
    
    parser.add_argument('--data-type', 
                       type=str,
                       required=False,
                       help='Data type (genome or transcriptome). Required for retrieval mode, optional for placement-only mode.')
    
    # Workflow mode arguments
    parser.add_argument('--mode',
                       type=str,
                       choices=['full', 'retrieval', 'placement'],
                       default='full',
                       help='Workflow mode: full (both parts), retrieval (SSU extraction + BUSCO), placement (phylogenetic placement only)')
    
    # Optional arguments
    parser.add_argument('-o', '--output-dir', 
                       default='screening_tool_output',
                       help='Output directory (default: screening_tool_output)')
    
    parser.add_argument('--pplacer-cutoff-length',
                       type=int,
                       default=500,
                       help='Cutoff length for pplacer (default: 500)')
    
    parser.add_argument('-t', '--threads', 
                       type=int, default=1,
                       help='Number of threads to use (default: 1)')
    
    parser.add_argument('--dry-run', 
                       action='store_true',
                       help='Show what would be run without executing')
    
    return parser

def validate_args(args):
    """Validate command-line arguments."""
    # Check input files exist
    if not os.path.isfile(args.fasta):
        print(f"Error: Input FASTA file '{args.fasta}' not found", file=sys.stderr)
        return False
    
    # Validate required supergroup
    valid_supergroups = ['Amoebozoa', 'Obazoa', 'Excavata', 'TSAR', 'Archaeplastida', 'Cryptista', 'Haptista', 'CRuMs', 'Provora']
    if args.supergroup not in valid_supergroups:
        print(f"Error: Invalid supergroup '{args.supergroup}'. Valid options: {', '.join(valid_supergroups)}", file=sys.stderr)
        return False
    
    # Validate data type
    valid_data_types = ['genome', 'transcriptome']
    
    # Data type is required for retrieval mode, but optional for placement mode
    if args.mode != 'placement':
        if not args.data_type:
            print(f"Error: --data-type is required for retrieval mode. Valid options: {', '.join(valid_data_types)}", file=sys.stderr)
            return False
        if args.data_type not in valid_data_types:
            print(f"Error: Invalid data type '{args.data_type}'. Valid options: {', '.join(valid_data_types)}", file=sys.stderr)
            return False
    elif args.data_type and args.data_type not in valid_data_types:
        # If data type is provided in placement mode, validate it
        print(f"Error: Invalid data type '{args.data_type}'. Valid options: {', '.join(valid_data_types)}", file=sys.stderr)
        return False
    
    # Validate workflow mode
    if args.mode == 'placement':
        print("Note: Assuming input file contains pre-extracted SSU sequences", file=sys.stderr)
    
    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    return True

def main():
    """Main CLI entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    # Validate arguments
    if not validate_args(args):
        sys.exit(1)
    
    # Configure runner
    runner = SnakemakeRunner(
        fasta=args.fasta,
        output_dir=args.output_dir,
        supergroup=args.supergroup,
        threads=args.threads,
        dry_run=args.dry_run,
        pplacer_cutoff_length=args.pplacer_cutoff_length,
        busco_mode=args.data_type,
        workflow_mode=args.mode,
        skip_retrieval=(args.mode == 'placement')
    )
    
    try:
        runner.run()
    except Exception as e:
        sys.exit(1)

if __name__ == '__main__':
    main()