# P10K SSU Screening Tool

A command-line tool for screening SSU (Small Subunit ribosomal RNA) sequences in genomic and transcriptomic data. This tool is part of the P10K project and helps identify and classify SSU sequences using phylogenetic placement methods.

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
pip install p10k-ssu-screening-tool
```

### From Source
```bash
git clone https://github.com/AlexTiceLab/P10K-SSU-screening-tool.git
cd P10K-SSU-screening-tool
pip install -e .
```

### Dependencies

The tool requires the following external programs to be installed and available in your PATH:

- **BLAST+** (makeblastdb, blastn)
- **MAFFT** (for sequence alignment)
- **pplacer** (for phylogenetic placement)
- **guppy** (for classification post-processing)
- **Snakemake** (workflow management)

You can install these using conda:
```bash
conda install -c bioconda blast mafft pplacer snakemake-minimal
```

## Quick Start

### Basic Usage
```bash
# Run screening on FASTA file with taxonomy information
p10k-ssu-screen input_sequences.fasta taxonomy_info.txt

# Specify output directory and number of threads
p10k-ssu-screen input_sequences.fasta taxonomy_info.txt -o results/ -t 8

# Run only 18S analysis
p10k-ssu-screen --18s-only input_sequences.fasta taxonomy_info.txt
```

### Alternative Execution Methods
```bash
# As a Python module
python -m p10k_ssu_screening input_sequences.fasta taxonomy_info.txt

# From source directory
python -m p10k_ssu_screening.cli input_sequences.fasta taxonomy_info.txt
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
  --18s-only           Run only 18S SSU screening
  --16s-only           Run only 16S SSU screening
  --dry-run            Show what would be run without executing
  --keep-temp          Keep temporary files
  -v, --verbose        Verbose output
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

### Taxonomy File
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
├── blast/             # BLAST results
├── alignment/         # MAFFT alignments
├── placement/         # pplacer phylogenetic placements
├── guppy/            # Classification results
└── summary/          # Final summary reports
```

### Main Output Files

- **summary_report.txt**: Main results with taxonomic classifications
- **placement.jplace**: Phylogenetic placement results (JSON format)
- **guppy_results.txt**: Detailed classification information

## Examples

### Example 1: Basic Screening
```bash
p10k-ssu-screen genome_assembly.fasta taxonomy.txt
```

### Example 2: Multi-threaded Analysis
```bash
p10k-ssu-screen transcriptome.fasta taxonomy.txt -o ssu_results -t 16 -v
```

### Example 3: 18S-only Analysis
```bash
p10k-ssu-screen eukaryotic_samples.fasta taxonomy.txt --18s-only -o euk_ssu
```

### Example 4: Dry Run (Check What Would Execute)
```bash
p10k-ssu-screen input.fasta taxonomy.txt --dry-run -v
```

## Configuration

The tool supports optional configuration files for advanced settings:

### Config File Locations (searched in order)
1. `~/.p10k_ssu_screening.yaml`
2. `./p10k_config.yaml`
3. `./.p10k_config.yaml`

### Example Configuration
```yaml
threads: 8
blast:
  evalue: 1e-5
  max_target_seqs: 10
mafft:
  algorithm: auto
  maxiterate: 1000
pplacer:
  keep_at_most: 7
```

## Development

### Running Tests
```bash
pip install -e ".[dev]"
pytest tests/
```

### Code Style
```bash
black p10k_ssu_screening/
flake8 p10k_ssu_screening/
```

## Troubleshooting

### Common Issues

1. **Command not found**: Ensure external dependencies (BLAST, MAFFT, pplacer) are installed and in PATH
2. **Permission errors**: Check write permissions for output directory
3. **Memory issues**: Reduce thread count or use smaller input files
4. **Reference package errors**: Ensure reference packages are properly installed

### Getting Help

- **GitHub Issues**: [Report bugs or request features](https://github.com/AlexTiceLab/P10K-SSU-screening-tool/issues)
- **Verbose output**: Use `-v` flag for detailed execution information
- **Dry run**: Use `--dry-run` to see what commands would be executed

## Citation

If you use this tool in your research, please cite:

```
P10K SSU Screening Tool
Jones, R.E. (2024)
GitHub: https://github.com/AlexTiceLab/P10K-SSU-screening-tool
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Authors

- **Robert E. Jones** - *Initial work* - [robert.ervin.jones@gmail.com](mailto:robert.ervin.jones@gmail.com)

## Acknowledgments

- P10K project team
- Alex Tice Lab
- Contributors to BLAST, MAFFT, pplacer, and Snakemake tools