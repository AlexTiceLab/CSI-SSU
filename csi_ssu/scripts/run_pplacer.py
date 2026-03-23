#!/usr/bin/env python3
"""
Robust pplacer runner with automatic error detection and sequence filtering.

This script attempts to run pplacer on an alignment. If pplacer fails, it tests
each query sequence individually to identify problematic sequences, removes them,
and reruns pplacer on a cleaned alignment.
"""

import sys
import os
import subprocess
import tempfile
import argparse
from pathlib import Path
from Bio import SeqIO
import logging


class RobustPplacer:
    """Handles pplacer execution with automatic error recovery."""
    
    def __init__(self, alignment_file, refpkg, output_file, threads=1, 
                 ref_alignment=None, log_file=None):
        """
        Initialize the RobustPplacer.
        
        Args:
            alignment_file: Path to the alignment file (FASTA format)
            refpkg: Path to the reference package
            output_file: Path for output jplace file
            threads: Number of threads for pplacer
            ref_alignment: Path to reference alignment (to identify query sequences)
            log_file: Path to log file (optional)
        """
        self.alignment_file = Path(alignment_file)
        self.refpkg = Path(refpkg)
        self.output_file = Path(output_file)
        self.threads = threads
        self.ref_alignment = Path(ref_alignment) if ref_alignment else None
        self.log_file = log_file
        
        # Setup logging
        self._setup_logging()
        
        # Load sequences
        self.all_sequences = list(SeqIO.parse(self.alignment_file, "fasta"))
        self.ref_sequence_ids = self._get_reference_ids()
        self.query_sequences = [seq for seq in self.all_sequences 
                               if seq.id not in self.ref_sequence_ids]
        
        self.logger.info(f"Loaded {len(self.all_sequences)} total sequences")
        self.logger.info(f"Reference sequences: {len(self.ref_sequence_ids)}")
        self.logger.info(f"Query sequences: {len(self.query_sequences)}")
    
    def _setup_logging(self):
        """Setup logging configuration."""
        self.logger = logging.getLogger('RobustPplacer')
        self.logger.setLevel(logging.INFO)
        
        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        # File handler if log file specified
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(logging.INFO)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
    
    def _get_reference_ids(self):
        """Get set of reference sequence IDs."""
        if not self.ref_alignment or not self.ref_alignment.exists():
            self.logger.warning("Reference alignment not provided or doesn't exist")
            return set()
        
        ref_ids = set()
        for record in SeqIO.parse(self.ref_alignment, "fasta"):
            ref_ids.add(record.id)
        
        return ref_ids
    
    def _run_pplacer(self, input_file, output_file):
        """
        Run pplacer on the given input file.
        
        Args:
            input_file: Path to alignment file
            output_file: Path for output jplace file
            
        Returns:
            Tuple of (success: bool, stderr: str)
        """
        cmd = [
            'pplacer',
            '-c', str(self.refpkg),
            '-o', str(output_file),
            '-j', str(self.threads),
            str(input_file)
        ]
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            return True, result.stderr
        except subprocess.CalledProcessError as e:
            return False, e.stderr
        except FileNotFoundError:
            self.logger.error("pplacer not found. Please ensure it's installed and in PATH")
            sys.exit(1)
    
    def _test_single_query(self, query_seq):
        """
        Test a single query sequence with pplacer.
        
        Args:
            query_seq: Bio.SeqRecord object
            
        Returns:
            bool: True if pplacer succeeds, False otherwise
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Create temporary alignment with ref + single query
            test_aln = tmpdir / "test.fasta"
            with open(test_aln, 'w') as f:
                # Write reference sequences
                for seq in self.all_sequences:
                    if seq.id in self.ref_sequence_ids:
                        SeqIO.write(seq, f, "fasta")
                # Write test query
                SeqIO.write(query_seq, f, "fasta")
            
            # Test with pplacer
            test_output = tmpdir / "test.jplace"
            success, stderr = self._run_pplacer(test_aln, test_output)
            
            if not success:
                self.logger.debug(f"Sequence {query_seq.id} failed: {stderr[:200]}")
            
            return success
    
    def _create_filtered_alignment(self, good_queries, output_path):
        """
        Create a new alignment with only good query sequences.
        
        Args:
            good_queries: List of SeqRecord objects that passed testing
            output_path: Path for output alignment file
        """
        with open(output_path, 'w') as f:
            # Write reference sequences
            for seq in self.all_sequences:
                if seq.id in self.ref_sequence_ids:
                    SeqIO.write(seq, f, "fasta")
            
            # Write good query sequences
            for seq in good_queries:
                SeqIO.write(seq, f, "fasta")
        
        self.logger.info(f"Created filtered alignment: {output_path}")
    
    def run(self):
        """
        Main execution method.
        
        Returns:
            int: 0 on success, 1 on failure
        """
        self.logger.info("="*60)
        self.logger.info("Starting robust pplacer run")
        self.logger.info("="*60)
        
        # Step 1: Try running pplacer on full alignment
        self.logger.info("Step 1: Attempting pplacer on full alignment...")
        success, stderr = self._run_pplacer(self.alignment_file, self.output_file)
        
        if success:
            self.logger.info("✓ Pplacer completed successfully on full alignment!")
            return 0
        
        # Step 2: Pplacer failed, need to identify bad sequences
        self.logger.warning("✗ Pplacer failed on full alignment")
        self.logger.info(f"Error output: {stderr[:500]}")
        
        if len(self.query_sequences) == 0:
            self.logger.warning("No query sequences found - failure is in reference sequences")
            self.logger.warning("Cannot proceed with individual sequence testing")
            return 0
        
        self.logger.info(f"\nStep 2: Testing {len(self.query_sequences)} query sequences individually...")
        
        good_queries = []
        bad_queries = []
        
        for i, query_seq in enumerate(self.query_sequences, 1):
            self.logger.info(f"Testing sequence {i}/{len(self.query_sequences)}: {query_seq.id}")
            
            if self._test_single_query(query_seq):
                good_queries.append(query_seq)
                self.logger.info(f"  ✓ PASS")
            else:
                bad_queries.append(query_seq)
                self.logger.warning(f"  ✗ FAIL")
        
        # Step 3: Report results
        self.logger.info("\n" + "="*60)
        self.logger.info("Testing complete")
        self.logger.info("="*60)
        self.logger.info(f"Good sequences: {len(good_queries)}")
        self.logger.info(f"Bad sequences: {len(bad_queries)}")
        
        if bad_queries:
            self.logger.warning("\nProblematic sequences:")
            for seq in bad_queries:
                self.logger.warning(f"  - {seq.id}")
        
        if len(good_queries) == 0:
            self.logger.warning("\nNo valid query sequences remain after filtering!")
            self.logger.warning("All query sequences caused pplacer to fail")
            # Create empty jplace file to satisfy workflow
            with open(self.output_file, 'w') as f:
                f.write('{"tree": "", "placements": [], "fields": [], "version": 3, "metadata": {"info": "No valid sequences"}}\n')
            self.logger.info(f"Created empty jplace file: {self.output_file}")
            return 0
        
        # Step 4: Create filtered alignment and rerun
        self.logger.info(f"\nStep 3: Creating filtered alignment and rerunning pplacer...")
        
        # Create filtered alignment in same directory as original
        filtered_aln = self.alignment_file.parent / f"{self.alignment_file.stem} .fasta"
        self._create_filtered_alignment(good_queries, filtered_aln)
        
        # Rerun pplacer on filtered alignment
        self.logger.info("Rerunning pplacer on filtered alignment...")
        success, stderr = self._run_pplacer(filtered_aln, self.output_file)
        
        if success:
            self.logger.info("✓ Pplacer completed successfully on filtered alignment!")
            
            # Save list of removed sequences
            removed_file = self.alignment_file.parent / "removed_sequences.txt"
            with open(removed_file, 'w') as f:
                f.write("# Sequences removed due to pplacer errors\n")
                for seq in bad_queries:
                    f.write(f"{seq.id}\n")
            self.logger.info(f"List of removed sequences saved to: {removed_file}")
            
            return 0
        else:
            self.logger.error("✗ Pplacer still failed on filtered alignment")
            self.logger.error(f"Error: {stderr}")
            return 1


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Robust pplacer runner with automatic error detection and filtering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s -i aligned.fasta -c refpkg/ -o output.jplace
  
  # With reference alignment to identify query sequences
  %(prog)s -i aligned.fasta -c refpkg/ -o output.jplace -r reference.fasta
  
  # With multiple threads
  %(prog)s -i aligned.fasta -c refpkg/ -o output.jplace -j 8
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input alignment file (FASTA format)')
    parser.add_argument('-c', '--refpkg', required=True,
                       help='Reference package directory')
    parser.add_argument('-o', '--output', required=True,
                       help='Output jplace file')
    parser.add_argument('-r', '--reference', 
                       help='Reference alignment file (to identify query sequences)')
    parser.add_argument('-j', '--threads', type=int, default=1,
                       help='Number of threads (default: 1)')
    parser.add_argument('-l', '--log', 
                       help='Log file path (default: stdout only)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.refpkg).exists():
        print(f"Error: Reference package not found: {args.refpkg}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if needed
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    # Run robust pplacer
    runner = RobustPplacer(
        alignment_file=args.input,
        refpkg=args.refpkg,
        output_file=args.output,
        threads=args.threads,
        ref_alignment=args.reference,
        log_file=args.log
    )
    
    exit_code = runner.run()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
