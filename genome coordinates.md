# Convert genome coordinates/assembly 
CrossMap (https://crossmap.sourceforge.net/#convert-bed-format-files)
Supports BED/BedGraph, GFF/GTF, BAM/SAM/CRAM, BigWig/Wig, VCF, and MAF format files

```bash
#Load the software
module load gcc/12.2.0 python/3.10.10
source /mnt/apps/custom/python-envs/crossmap/bin/activate

# Software usage
CrossMap -h

cd path/to/your/input/files
CrossMap bed oviAri4ToGCF_016772045.1.over.chain.gz test.bed
