from setuptools import setup, find_packages
import os

# Read long description from README
here = os.path.abspath(os.path.dirname(__file__))
try:
    with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "SSU screening tool for P10K project genomic and transcriptomic data"

setup(
    name='p10k-ssu-screening-tool',
    version='1.0.0',
    description='SSU screening tool for P10K project genomic and transcriptomic data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/AlexTiceLab/P10K-SSU-screening-tool',
    author='Robert E. Jones',
    author_email='robert.ervin.jones@gmail.com',
    license='MIT',
    
    # Package discovery
    packages=find_packages(),
    
    # Python version requirement
    python_requires='>=3.13',   
    
    # Command-line entry point
    entry_points={
        'console_scripts': [
            'p10k-ssu-screen=p10k_ssu_screening.cli:main',
        ],
    },
    
    # Include data files
    package_data={
        'p10k_ssu_screening': [
            'workflows/*.smk',
            'data/reference_packages/**/*',
            'data/queries/*',
        ],
    },
    include_package_data=True,
    
    # Keywords for discovery
    keywords='bioinformatics genomics transcriptomics SSU screening phylogeny P10K',
    
    # Project URLs
    project_urls={
        'Bug Reports': 'https://github.com/AlexTiceLab/P10K-SSU-screening-tool/issues',
        'Source': 'https://github.com/AlexTiceLab/P10K-SSU-screening-tool',
    },
)
