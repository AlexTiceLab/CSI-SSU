# CSI-SSU

A command-line tool for screening SSU (Small Subunit ribosomal RNA) sequences in genomic and transcriptomic data. This tool helps identify and classify SSU sequences using phylogenetic placement methods.

## Features

- **Automated SSU screening** using BLAST, MAFFT, and pplacer
- **Phylogenetic placement** for accurate taxonomic classification
- **18S and 16S analysis** with dedicated reference packages
- **Flexible output formats** including summary reports
- **Command-line interface** for easy integration into pipelines
- **Configurable parameters** for different analysis needs

## Installation

### From PyPI (when published)
```bash
pip install csi-ssu
```

### From Source
```bash
git clone https://github.com/AlexTiceLab/CSI-SSU.git
cd CSI-SSU
pip install -e .
```

### Dependencies

The tool requires the following external programs to be installed and available in your PATH:

- **BLAST+** (makeblastdb, blastn)
- **MAFFT** (for sequence alignment)
- **pplacer** (for phylogenetic placement)
- **cd-hit** (for clustering sequences)
- **Snakemake** (workflow management)

You can install these using conda:
```bash
conda install -c bioconda blast mafft pplacer snakemake-minimal
```

## Quick Start

### Basic Usage
```bash
# Run screening on FASTA file with taxonomy information
csi-ssu input_sequences.fasta taxonomy_info.txt

# Specify output directory and number of threads
csi-ssu input_sequences.fasta taxonomy_info.txt -o results/ -t 8

# Run only 18S analysis
csi-ssu --18s-only input_sequences.fasta taxonomy_info.txt
```

## Command-Line Options

```
positional arguments:
  fasta                 Input FASTA file path
  taxonomy              Input taxonomic classification file

optional arguments:
  -h, --help           show this help message and exit
  --version            show program's version number and exit
  -o, --output-dir     Output directory (default: screening_tool_output)
  -t, --threads        Number of threads to use (default: 1)
  --dry-run            Show what would be run without executing
```

## Input Files

### FASTA File
Standard FASTA format containing sequences to be screened:
```
>sequence1
ATCGATCGATCGATCG...
>sequence2
GCTAGCTAGCTAGCTA...
```

### Taxonomy File (FIX THIS ONCE WE HAVE FINALIZED THE FORMAT)
Tab-separated file with taxonomic information:
```
sequence_id    taxonomy_string
sequence1      Eukaryota;Amoebozoa;...
sequence2      Bacteria;Proteobacteria;...
```

## Output

The tool creates a structured output directory:

```
screening_tool_output/
├── logs/              # Log files for each step
├── blast_db/          # BLAST database files
├── blast/             # BLAST results
├── parsed_blast/      # Parsed BLAST hits
├── cd-hit/            # CD-HIT clustered sequences
├── mafft/             # MAFFT alignments
├── pplacer/           # pplacer phylogenetic placements
├── guppy/             # Classification results
└── summary/           # Final summary reports
```

### Main Output Files

- **summary_report.txt**: Main results with taxonomic classifications
- **placement.jplace**: Phylogenetic placement results (JSON format)
- **guppy_results.txt**: Detailed classification information

## Examples

### Example 1: Basic Screening
```bash
csi-ssu genome_assembly.fasta taxonomy.txt
```

### Example 2: Multi-threaded Analysis
```bash
csi-ssu transcriptome.fasta taxonomy.txt -o ssu_results -t 16 -v
```

### Example 3: 18S-only Analysis
```bash
csi-ssu eukaryotic_samples.fasta taxonomy.txt --18s-only -o euk_ssu
```

### Example 4: Dry Run (Check What Would Execute)
```bash
csi-ssu input.fasta taxonomy.txt --dry-run -v
```


### Getting Help

- **GitHub Issues**: [Report bugs or request features](https://github.com/AlexTiceLab/CSI-SSU/issues)
- **Dry run**: Use `--dry-run` to see what commands would be executed

## License

This project is licensed under the MIT License - see the LICENSE file for details.
