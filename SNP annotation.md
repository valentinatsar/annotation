# SNP annotation - Variant Effect Predictor 
(running VEP https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)

```bash
# Load miniconda and create a new environment
module load gcc/13.2.0-nbog6z2 miniconda3/24.3.0-lxaq4g2
source $CONDA_PROFILE/conda.sh

conda create --name annotation
conda activate annotation

#Check if the software is available
conda list | grep ensembl-vep

# Install VEP
conda install bioconda::ensembl-vep=113.3

#Install a species specific cache file
vep_install -a cf -s ovis_aries -y /path/to/ARS-UI_Ramb_v2.0 --CONVERT #For sheep (genome assembly: ARS-UI_Ramb_v2.0)
vep_install -a cf -s capra_hircus -y /path/to/ARS1 --CONVERT #For goat (genome assembly: ARS1)```

3. Prepare Your SNP List - Input files

For example, in default VEP format (https://www.ensembl.org/info/docs/tools/vep/vep_formats.html#default)
