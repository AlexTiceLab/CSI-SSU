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
  csi-ssu input.fasta -o results/ -t 8 --supergroup Amoebozoa
        """
    )
    
    # Version
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    
    # Required arguments
    parser.add_argument('fasta', help='Input FASTA file path')

    parser.add_argument('--supergroup', 
                       type=str,
                       required=True,
                       help='Supergroup of interest (e.g., Amoebozoa, Excavata, TSAR)')
    
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
    valid_supergroups = ['Amoebozoa', 'Obazoa', 'Excavata', 'TSAR', 'Archaeplastida', 'Cryptista', 'Haptista', 'Eukaryota_X', 'CRuMs', 'Provora']
    if args.supergroup not in valid_supergroups:
        print(f"Error: Invalid supergroup '{args.supergroup}'. Valid options: {', '.join(valid_supergroups)}", file=sys.stderr)
        return False
    
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
        pplacer_cutoff_length=args.pplacer_cutoff_length
    )
    
    try:
        runner.run()
    except Exception as e:
        sys.exit(1)

if __name__ == '__main__':
    main()