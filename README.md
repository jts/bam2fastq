bam2fastq
=========
    
bam2fastq is a program to extract sequences and qualities from a BAM file.

The original version can be found [here](http://gsl.hudsonalpha.org/information/software/bam2fastq).
It has subsequently been modified to handle BAM files with mixtures of paired and unpaired reads and write to stdout.

------------
INSTALLATION

Requires make, gcc, and the zlib compression libraries - all of which should
be present on most Unix-like systems (including Macs).

First clone the repository from github, including its dependencies:

```
git clone --recursive https://github.com/jts/bam2fastq
```

Then change to the directory and run make:

```
cd bam2fastq
make
```

------
AUTHOR

This program was originally developed at HudsonAlpha by Phillip Dexheimer. 
It was packaged on github and slightly modified by Jared Simpson.
