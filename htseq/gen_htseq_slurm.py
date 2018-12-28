#!/usr/bin/env python3

import os
import glob


for sampdir in glob.glob("/staging/wl/caim/RNA_seq_KO/1.rawdata/*/"):
	sampname = os.path.basename(sampdir.strip("/"))
	print(sampdir)
	slurmf = "htseq_{}.slurm".format(sampname)
	with open(slurmf, "w") as f:
		f.write("#!/bin/bash\n\n")

		f.write("#SBATCH --job-name=htseq_{}\n".format(sampname))
		f.write("#SBATCH --ntasks=2\n")
		f.write("#SBATCH --mem-per-cpu=4GB\n")
		f.write("#SBATCH --time=23:59:59\n\n")

		f.write("SAMPLE={}\n".format(sampname))
		f.write("SAMTOOLS=/home/rcf-proj/wl/caim/usr/samtools-1.3.1/samtools\n")
		f.write("HTSEQ_COUNT=/home/rcf-proj/wl/caim/usr/anaconda2/bin/htseq-count\n")
		f.write("cd /staging/wl/caim/analysis/22Rv1_KO_RNA_Seq/\n")
		f.write("cp /staging/wl/caim/seqlib/Homo_sapiens/GENCODE/gencode.v19.annotation.gtf $SCRATCHDIR/\n")
		f.write("$SAMTOOLS sort -n -o $SAMPLE/accepted_hits_sort.bam $SAMPLE/accepted_hits.bam\n")
		f.write("$HTSEQ_COUNT -f bam -r name -s reverse -t exon -i gene_name $SAMPLE/accepted_hits_sort.bam $SCRATCHDIR/gencode.v19.annotation.gtf > $SAMPLE/htseq_$SAMPLE.txt\n")
