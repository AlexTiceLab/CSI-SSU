#!/usr/bin/env python3
"""Command-line interface for P10K SSU screening tool."""

import argparse
import sys
import os
from pathlib import Path
from .runner import SnakemakeRunner
from . import __version__

def create_parser():
    """Create and configure argument parser."""
    parser = argparse.ArgumentParser(
        prog='p10k-ssu-screen',
        description='SSU screening tool for P10K project genomic and transcriptomic data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  p10k-ssu-screen input.fasta taxonomy.txt
  p10k-ssu-screen input.fasta taxonomy.txt -o results/ -t 8
  p10k-ssu-screen --18s-only input.fasta taxonomy.txt
        """
    )
    
    # Version
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    
    # Required arguments
    parser.add_argument('fasta', help='Input FASTA file path')
    parser.add_argument('taxonomy', help='Input taxonomic classification file')
    
    # Optional arguments
    parser.add_argument('-o', '--output-dir', 
                       default='screening_tool_output',
                       help='Output directory (default: screening_tool_output)')
    
    parser.add_argument('-t', '--threads', 
                       type=int, default=1,
                       help='Number of threads to use (default: 1)')
    
    parser.add_argument('--18s-only', 
                       action='store_true',
                       help='Run only 18S SSU screening')
    
    parser.add_argument('--16s-only', 
                       action='store_true',
                       help='Run only 16S SSU screening')
    
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
        
    if not os.path.isfile(args.taxonomy):
        print(f"Error: Taxonomy file '{args.taxonomy}' not found", file=sys.stderr)
        return False
    
    # Check for conflicting options
    if getattr(args, '18s_only', False) and getattr(args, '16s_only', False):
        print("Error: Cannot specify both --18s-only and --16s-only", file=sys.stderr)
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
        taxonomy=args.taxonomy,
        output_dir=args.output_dir,
        threads=args.threads,
        only_18s=getattr(args, '18s_only', False),
        only_16s=getattr(args, '16s_only', False),
        dry_run=args.dry_run,
    )
    
    try:
        runner.run()
    except Exception as e:
        sys.exit(1)

if __name__ == '__main__':
    main()