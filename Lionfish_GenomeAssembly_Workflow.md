# Tackling the Lionfish Genome

_Goal:_ For the Vertebrate Genome Project, to have a high quality, reference level genome for one individual representing every major group of vertebrates. For MGL, to potentially use gene editing technologies for population control of a highly invasive species throughout most of the Atlantic.

_Data:_ 60X PacBio, 100X Illumina, Hi-C, RNA-Seq

1. Quality Control

   - Purpose: Screen all data (of all data types) for outlier libraries, sequencing runs and contamination

   - Tools:

 [__Mash__](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x), reduces large sequences/sequence sets to compressed representative sketches
[Github](https://github.com/marbl/mash)

 Start: to run start with 21mers to generate sketches with sketch size of 10,000, compare between sequencing runs and then sequencing sets

2. Genome size, repeat content, and heterozygosity estimates

   - Purpose: Will affect how we choose to do most downstream analysis

   - Tools:

	 __scaff10x,__ trim barcodes during preprocessing
[Github](https://github.com/wtsi-hpag/Scaff10X)

	 __Meryl,__ kmer counter, sorts kmers
[Github](https://github.com/marbl/meryl)

	__[GenomeScope,](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5870704/pdf/btx153.pdf)__ quickly estimates genome size, heterozygosity and repeat content *(Rhie et al. mention that they got ~double the genome size with this tool so be sure to check)*
[Github](https://github.com/tbenavi1/genomescope2.0)

3. Contigging & Scaffolding

  - Tools:

	PacBio LR

[__FALCON & FALCON-Unzip__](https://www.pacb.com/products-and-services/analytical-software/smrt-analysis/), denovo assembler for PBLR, generate primary contigs (use __Arrow__ to fill gaps, details below), decide min read length (*Rhie et al. use 2kbp*), __FALCON-Unzip__ uses the primary contigs from __FALCON__ to generate haplotigs (divergent haplotypes from primary contigs). Run with default params first until calculating contig overlap.

[__Background__](https://pb-falcon.readthedocs.io/en/latest/about.html)

[Github](https://github.com/PacificBiosciences/FALCON) and [FALCON-integrate](https://github.com/PacificBiosciences/FALCON-integrate)

Requirements: [SMRTanalysis suite](https://www.pacb.com/products-and-services/analytical-software/smrt-analysis/)

‘nucmer’

‘samtools’

‘pbalign’

‘variantCaller’

‘pysam’ (python)

[__Canu__](https://canu.readthedocs.io/en/latest/quick-start.html), specific to assemblies for PacBio and Oxford Nanopore
[Github](https://github.com/marbl/canu)

[__Flye__](https://www.nature.com/articles/s41587-019-0072-8), also specific to PacBio and ONT data, this worked best with our beetle ONT data. We could run all three and compare initial assembly stats (wtdbg2 could be another option, but I think is outperformed by Flye in nearly every comparison plus Flye runs faster).
	[Github](https://github.com/fenderglass/Flye)

[*Kolmogorov et al.](https://www.nature.com/articles/s41587-019-0072-8.pdf) show Flye and Canu producing the most contiguous assemblies, but with Canu having the highest number of missasemblies when compared with Flye and FALCON. Flye also seems to be the fastest of the three, particularly as genome size increases, so there may be a better indicator of the best option for assembly once we know genome size for the lionfish.

Illumina SR assembly

[__SoapDeNovo2,__](https://www.animalgenome.org/bioinfo/resources/manuals/SOAP.html) SR assembler (newer)
	 [Github](https://github.com/aquaskyline/SOAPdenovo2)

[__Supernova,__](https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance) SR assembler (worked quite well in my experience, although do we have 10X data for the lionfish?)
[Github](https://github.com/10XGenomics/supernova)

Hi-C reads
Map to PacBio assembly to scaffold and order contigs, creating chromosome length scaffolds.

ArimaGenomics has a mapping pipeline for processing alignments
[Github](https://github.com/ArimaGenomics/mapping_pipeline)

[__Salsa2.0,__](https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1007273&type=printable) used to scaffold LR with Hi-C data (HiRise and Proximo are also options, but produce similar scaffolding stats to Salsa and are not open access, according to Rhie et al.)
[Github](https://github.com/marbl/SALSA)

4. Gap Filling

 - Tools:

 	[__BAMTOOLS,__](https://bioinformatics.readthedocs.io/en/latest/bamtools/) extract FASTQ sequences from PacBio BAM files

	[__PBJelly,__](https://github.com/esrice/PBJelly) for LR (Rhie et al. found that the # of gaps closed with PBJelly was similar to that with Arrow, so PBJelly may not be a necessary step in the pipeline, but PBJelly is more reported in the lit so might be best to use it anyways)

	[__Abyss Sealer,__](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0663-4) SR assembler and gap sealer (what a great name), works particularly well for large genomes
[Github](https://github.com/bcgsc/abyss/tree/sealer-release)

5. Polishing

 - Tools:

	[__FreeBayes__](https://arxiv.org/pdf/1207.3907.pdf), SR polishing ([__Pilon__](https://github.com/broadinstitute/pilon) is a good option too but not great for large [genomes](http://michaelalonge.com/comparing-short-read-polishers.html), so another thing we can troubleshoot when we have a genome size estimate)

	__Arrow,__ for PacBio LR polishing

6. Annotation with RNA-Seq Data

	[__MAKER,__](https://genome.cshlp.org/content/18/1/188.full.pdf+html) bioinformatic tool for genome annotation
[Github](https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2)

__Other Resources:__
-	[Cool camel paper that had a similar dataset to the Lionfish](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6618069/pdf/MEN-19-1015.pdf)


-	[Yellowbelly pufferfish genome assembly with similar dataset](https://www.nature.com/articles/s41597-019-0279-z.pdf)


-	[A Roadmap to denovo genome assembly](https://cehg.stanford.edu/sites/g/files/sbiybj9906/f/assembly_workshop_stanford_2017_v_1.pdf)


-	[There are really just so many options](https://warrenlr.github.io/papers/DeNovoAssemblyBTL.pdf)


__Example workflows:__

[Rhie et al. 2020, Vertebrate Genome Project (VGP) workflow](https://www.biorxiv.org/content/10.1101/2020.05.22.110833v1.full.pdf)



[Zhou et al. 2019, Yellowbelly Pufferfish workflow](https://www.nature.com/articles/s41597-019-0279-z.pdf)
