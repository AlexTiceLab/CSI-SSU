from setuptools import setup, find_packages
import os

# Read long description from README
here = os.path.abspath(os.path.dirname(__file__))
try:
    with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = f.read()
except FileNotFoundError:
    long_description = "SSU screening tool for genomic and transcriptomic data"

setup(
    name='csi-ssu',
    version='1.0.0',
    description='CSI CSI SSU screening tool for genomic and transcriptomic data',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/AlexTiceLab/CSI-SSU',
    author='Robert E. Jones',
    author_email='robert.ervin.jones@gmail.com',
    license='MIT',
    
    # Package discovery
    packages=find_packages(),
    
    # Command-line entry point
    entry_points={
        'console_scripts': [
            'csi-ssu=csi_ssu_screening.cli:main',
        ],
    },
    
    # Include data files
    package_data={
        'csi_ssu_screening': [
            'workflows/*.smk',
            'data/reference_packages/**/*',
            'data/queries/*',
        ],
    },
    include_package_data=True,
    
    # Keywords for discovery
    keywords='bioinformatics genomics transcriptomics SSU-screening',
    
    # Project URLs
    project_urls={
        'Bug Reports': 'https://github.com/AlexTiceLab/CSI-SSU/issues',
        'Source': 'https://github.com/AlexTiceLab/CSI-SSU',
    },
)
