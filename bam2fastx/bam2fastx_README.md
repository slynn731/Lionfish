# bam2fastx
  Converting PacBio bam output files to fasta or fastq files.

  You may first need to make sure [__Miniconda__](https://docs.conda.io/en/latest/miniconda.html) is installed on the machine.

  If so, perform the following in Linux:

>curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

>sh Miniconda3-latest-Linux-x86_64.sh

  Follow the installation instructions that are provided from this command.

  Once installation is complete, exit and log back into the system.

  You will then need to install the [__Bioconda__](https://bioconda.github.io/) package.

To do that run these commands:
>conda config --add channels defaults

>conda config --add channels bioconda

>conda config --add channels conda-forge

Then perform:

> conda install bam2fastx

A list of available packages can be found [__here__](https://github.com/PacificBiosciences/pbbioconda)
