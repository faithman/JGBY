#!/usr/bin/env bash

nextflow run joint_sv.nf -profile quest --reference "/projects/b1059/projects/Ye/snpEff/data/Chicken.GRCg6a/genomes/Chicken.GRCg6a.fa" --bamdir '/projects/b1059/projects/Ye/repo/ch-calling-nf/CH-20190730/BAM/sv_bam/*.{bam,bam.bai}' --email "yewangfaith@gmail.com" -resume
