# Convert genome coordinates/assembly 
CrossMap - genome coordinates conversion between different assemblies ((https://crossmap.sourceforge.net/)

Supports BED/BedGraph, GFF/GTF, BAM/SAM/CRAM, BigWig/Wig, VCF, and MAF format files

```bash
#Load the software
module load gcc/12.2.0 python/3.10.10
source /mnt/apps/custom/python-envs/crossmap/bin/activate

# Software usage
CrossMap -h

# Convert genome coordinates using the appropriate chain file
CrossMap bed path/to/oviAri4ToGCF_016772045.1.over.chain.gz path/to/test.bed
