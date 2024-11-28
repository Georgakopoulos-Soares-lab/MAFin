from setuptools import setup, find_packages

setup(
    name='MAFin',  # Package name
    version='0.6',  # Version number
    description='Motif Detection in Multiple Alignment Files: CLI tool to process and analyze MAF files ( Multiple alignment format ) ',  # Short description
    long_description=open('README.md').read(),  # Full description from README.md
    long_description_content_type='text/markdown',  # Content type for long description
    author='Patsakis Michail , Provatas Kimon, Ilias Georgakopoulos Soares, Ioannis Mouratidis',  # Your name
    author_email='kap6605@psu.edu , mpp5977@psu.edu',  # Your email
    url='https://github.com/Georgakopoulos-Soares-lab/MAFin',  # URL to your project (optional)
    packages=find_packages(),  # Automatically find all packages in the project
    entry_points={
        'console_scripts': [
            'MAFin=mafin.mafin:main',  # Register CLI command
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # Specify Python version requirement
    install_requires=[
        'biopython' ,
        'pyahocorasick'
    ],
)