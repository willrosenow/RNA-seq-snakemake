# RNA-seq-snakemake
This automated RNA-seq analysis pipeline uses Hisat2 for alignment and Stringie for transcript assembly and quantification. It outputs gtf and abundance files for each sample, along with overall gene count and trascript count matricies that can be read into DESeq2 for differential expression analysis. This workflow takes trimmed fastq files as inputs, therefore trimming of adaptor sequences must be done prior. 

# Installation Instructions:
1. Open your terminal and clone this repository to your computer using the command: `git clone https://github.com/willrosenow/RNA-seq-snakemake.git`
2. This pipeline uses a conda environment to ensure package dependencies work correctly. Create a conda environment named rnaseq by running: `conda env create -f environment.yaml`
3. Once this is complete, activate the conda environment using: `conda activate rnaseq`. To deactivate the environment simply type: `conda deactivate`
4. Finally, we must install Snakemake using mamba. To install snakemake, we must activate the conda environment as shown above. Then type: `mamba install -c conda-forge -c bioconda snakemake`
5. Once this is complete your conda envioronment is set up and ready to run the pipeline.

# Pipeline Instructions:
1. Hisat2 requires index files for alignment. I download these from the [Hisat2 website](http://daehwankimlab.github.io/hisat2/download/). Find the genome_snp index for your organism/build and copy the link address. To do this from the command line, type: `get--content-disposition https://cloud.biohpc.swmed.edu/index.php/s/grcm38_snp/download`. This will download a .tar.gz file.
2. Unpack this file by running: `tar -xzvf grcm38_snp.tar.gzgrcm38_snp.tar.gz`. This will create a directory called grcm38_snp with all the index files in it.
3. Next, we need to download the GTF file for Stringtie. An easy way to do this is using the [UCSC website](https://hgdownload.soe.ucsc.edu/downloads.html#mouse). To download the mm39 GTF file, type: `wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz`. 
4. Finally, create a data directory to store your trimmed fastq files as inputs. You can do this by typing: `mkdir data`. This pipeline does not run trimmomatic, so make sure you trim adaptors from your files before running this pipeline. Each sample should have a R1 and R2 file associated with it. 

## Edit the config.yaml file:
5. Edit the config.yaml to insert your sample names and paths, GTF path, Hisat index path, threads, etc. There is no need to make changes to the Snakefile. 
## Dry-run:
6. Now you are ready to run the pipeline. It's a good idea to do a dry-run, which will just print the commands to screen. This is helpful to make sure all the paths are correct. To do a dry-run type: `snakemake -np`
## Run snakemake:
7. Once everything is correct, make sure your rnaseq conda environment is activated and type: `snakemake -p -j 1`. -j specifies the number of cores, so set this to match what you put in the config.yaml file. The -p option simply prints the commands to the screen. 
## Outputs:
8. The pipeline will create directories for each sample `stringtie/quant/exp1`. Each directory will contain mulitple output files including the gtf and abundance files. Overall gene count and transcript count matricies will be ouput in `differential_expression/`. These can be directly read into DESeq2 for differential expression analysis.

### Feel free to send any questions to Will Rosenow at wr8yp@virginia.edu

