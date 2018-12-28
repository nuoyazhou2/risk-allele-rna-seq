#!/usr/bin/env python3

import os
import glob


for sampdir in glob.glob("/staging/wl/caim/RNA_seq_KO/1.rawdata/*/"):
	sampname = os.path.basename(sampdir.strip("/"))
	print(sampdir)
	slurmf = "align_{}.slurm".format(sampname)
	with open(slurmf, "w") as f:
		f.write("#!/bin/bash\n\n")

		f.write("#SBATCH --job-name=align_{}\n".format(sampname))
		f.write("#SBATCH --ntasks=10\n")
		f.write("#SBATCH --time=23:59:59\n\n")

		f.write("FQ1={}_1.fq.gz\n".format(sampname))
		f.write("FQ2={}_2.fq.gz\n".format(sampname))
		f.write("OUT={}\n".format(sampname))
		f.write("TOPHAT=/home/rcf-proj/wl/caim/usr/tophat-2.1.1.Linux_x86_64/tophat2\n\n")

		f.write("cp /staging/wl/caim/seqlib/Homo_sapiens/GENCODE/gencode.v19.annotation.gtf $SCRATCHDIR/\n")
		f.write("cp /staging/wl/caim/seqlib/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome* $SCRATCHDIR/\n")
		f.write("cp /staging/wl/caim/RNA_seq_KO/1.rawdata/{}/$FQ1 $SCRATCHDIR/\n".format(sampname))
		f.write("cp /staging/wl/caim/RNA_seq_KO/1.rawdata/{}/$FQ2 $SCRATCHDIR/\n\n".format(sampname))

		f.write("cd /staging/wl/caim/analysis/22Rv1_KO_RNA_Seq/\n")
		f.write("$TOPHAT --library-type fr-firststrand -p 10 -G $SCRATCHDIR/gencode.v19.annotation.gtf -o $OUT $SCRATCHDIR/genome $SCRATCHDIR/$FQ1 $SCRATCHDIR/$FQ2\n\n")

		f.write("fastqc -o FastQC -f fastq $SCRATCHDIR/$FQ1 $SCRATCHDIR/$FQ2\n")
