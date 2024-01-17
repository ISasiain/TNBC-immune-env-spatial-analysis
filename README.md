# Spatial analysis of the tumour immune microenvironment in Triple-Negative Breast Cancer

* Author: Iñaki Sasiain Casado
* Supervisors: Johan Staaf and Suze Roostee

### ANALYSING TMArQ PROCESSED SINGLE-PLEX SAMPLES

1. Generating spatial_experiment files from coordinates files.

* Preprocessing files' names to homogenize their format

```bash
for file in $(ls ./*_coordinates.txt | grep "CD8_"); 
    do new_filename=$(echo ${file} | sed 's/CD8_/CD8/');
       mv ./${file} ./${new_filename};
    done;

for file in $(ls ./*_coordinates.txt | grep "CD4_"); 
    do new_filename=$(echo ${file} | sed 's/CD4_/CD4/');
       mv ./${file} ./${new_filename};
    done;

for file in $(ls ./*_coordinates.txt | grep "CD3_"); 
    do new_filename=$(echo ${file} | sed 's/CD3_/CD3/');
       mv ./${file} ./${new_filename};
    done;


for file in $(ls ./*_coordinates.txt | grep "CD3TNBC"); 
    do new_filename=$(echo ${file} | sed 's/CD3TNBC/CD3_TNBC/');
       mv ./${file} ./${new_filename};
    done;


for file in *\ *; do
    # Check if the file has a space in its name
    if [[ -f "$file" ]]; then
        # Remove spaces by replacing them with an underscore
        newname="${file// /}"
        
        # Rename the file
        mv "$file" "$newname"
        
        # Display the changes made
        echo "Renamed '$file' to '$newname'"
    fi
done;
```

* Generating spatial experiment objects (Run in corsaire)

``` bash
conda install -c conda-forge r-optparse;
conda install -c conda-forge r-magick;
conda install -c conda-forge r-biocmanager;
conda install -c bioconda bioconductor-spatialexperiment;
conda install -c conda-forge r-terra;
conda install -c r r-raster;
conda install -c conda-forge r-ggplot2;
conda install -c conda-forge r-r.utils;
conda install r-alphahull;
conda install -c conda-forge r-tidyr;
conda install -c bioconda bioconductor-complexheatmap;



```R
install.packages("plotly")
```

```bash
cd /home/Illumina/Iñaki_Sasiain/immune_spatial/spe_objects; 

Rscript ../scripts/create_spe_objects.r -m p53,CD3,CD4,CD68,CD8,FOXP3,CD20 -a ../annotation/supplData_withimages.csv -p ../coordinates/;
```

2. Preliminary data analysis

* TILs and survival

* C



3. Analysing metrics in Basal or NonBasal TNBC


4. Anaysing clustering

- Generating DIN matrices for every sample

```bash
cd /home/Illumina/Iñaki_Sasiain/immune_spatial/analyse_clusters/DIN; 

#Create a variable with the comma separated paths to spe objects
spe_paths=$(find ../../spe_objects | tr "\n" ",");

#Generating DIN matrices for all the samples
Rscript ../../DIN_calculator.r -c 40 -o ${spe_paths::-1};

```