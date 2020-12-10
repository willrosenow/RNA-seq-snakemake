# RNA-seq-snakemake
Automated RNA-seq analysis pipeline using Snakemake

# Installation instructions:
1. Clone this repository to your computer using the command `git clone https://github.com/willrosenow/RNA-seq-snakemake.git`
2. This pipeline uses a conda environment to ensure package dependencies work correctly. Create a conda environment by running: `conda env create -f environment.yaml`
3. Once this is complete, activate the conda environment using: `conda activate rnaseq`
4. Finally, we must install Snakemake using mamba. To install snakemake, we must activate the conda environment as shown above. Then type `mamba install snakemake`
5. Once this is complete your conda envioronment is set up and ready to run the pipeline.

# Pipeline Instructions:
1. Hisat2 requires index files. I download these from the hisat2 website at http://daehwankimlab.github.io/hisat2/download/. Find the genome_snp index for your organism/build and copy the link address. To do this from the command line, type `get--content-disposition https://cloud.biohpc.swmed.edu/index.php/s/grcm38_snp/download`. This will download a .tar.gz file.
2. Unpack this file by running: `tar -xzvf grcm38_snp.tar.gzgrcm38_snp.tar.gz`. This will create a directory called grcm38_snp with all the index files in it.
3. 

