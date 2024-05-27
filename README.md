# Spatial analysis of the tumour immune microenvironment in Triple-Negative Breast Cancer

* Author: Iñaki Sasiain Casado (inaki.sasiain_casado@med.lu.se)
* Supervisors: Johan Staaf and Suze Roostee
* Master's Programme in Bioinformatics
* BINP 51 (45 ECTS)

### SUMMARY

This GitHub repository contains the scripts used and workflow followed for the completion of the "Spatial analysis of the Tumour immune microenvironment in Triple-Negative Breast Cancer" master's thesis (Master's degree in Bioinformatics, Lund University). This thesis consists of the analysis of both single-plex IHC images and multiplex IHC images using already existing computational tools such as TMArQ, or spatial data analysis tools such as SPIAT, and other generated new methods. 

### SCRIPTS

The main scripts included and used for this project are the following: 

* **create_spe_objects.r**. This R script transform the output produced by TMArQ into Spatial Experiment objects for them to be compatible with spatianl analysis tools such as SPIAT.
* **spe_objects_from_vectra.r**. This R script generates SPE objects from the output of the mIHC cell segmentation and phenoyping script designed for PhenoImager (Vectra Polaris) images.
* **segmentation_and_phenotyping.py**. This Python script performs annotation-free cell segmentation and phenotyoing for PhenoImager mIHC images.
* **clustering.r**. This R script perform cell environment detection and phenotyping based on the DIN values determined from every mIHC or sIHC cores.
* **DIN_calculator.r**. This R script calculates Density In the Neighbourhood (DIN) matrices for every core (SPE object) to be used to identify cell environments using the clustering.r script.
* **compare_counts_with_annotations.r**. This R script compares the output (cell counts) generated from the segmented and phenotyed PhenoIMager mIHC images with the sample annotations, analysing survival, counts per subtype, TILs, pathologist's B cell scores, etc. and generating plots.
* **comparing_counts.r**. This R script compares the counts obtained using the designed cell segmentation and phenotyping approach with the values obtained by TMArQ.
* **cell_counter.r**. This R script generates a summary dataframe with the cell counts detected per core in each of the TMAs analysed.
* **analyse_sIHC_clusters.r**. This script analyses the composition of the detected sIHC envoironmnets.

### SETTING UP CONDA ENVIRONMENT

```bash
# Create a conda environment that replicates the one used for the original analysis
conda env create -f spatial_analysis.yml
```

### EXPERIMENTAL PROCEDURE

#### Analysing TMArQ processed sIHC data

1. Generating Spatial Experiment (SPE) objects from coordinates files

* Preprocessing file names to homogenize their format

```bash
# CD8+ cell cordinate files
for file in $(ls ./*_coordinates.txt | grep "CD8_"); 
    do new_filename=$(echo ${file} | sed 's/CD8_/CD8/');
       mv ./${file} ./${new_filename};
    done;


# CD4+ cell cordinate files
for file in $(ls ./*_coordinates.txt | grep "CD4_"); 
    do new_filename=$(echo ${file} | sed 's/CD4_/CD4/');
       mv ./${file} ./${new_filename};
    done;


# CD3+ cell cordinate files
for file in $(ls ./*_coordinates.txt | grep "CD3_"); 
    do new_filename=$(echo ${file} | sed 's/CD3_/CD3/');
       mv ./${file} ./${new_filename};
    done;

for file in $(ls ./*_coordinates.txt | grep "CD3TNBC"); 
    do new_filename=$(echo ${file} | sed 's/CD3TNBC/CD3_TNBC/');
       mv ./${file} ./${new_filename};
    done;

# Removing spaces
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
* Generating SPE objects (It was run on the server)

```bash
cd /home/Illumina/Iñaki_Sasiain/immune_spatial/spe_objects; 

Rscript ../scripts/create_spe_objects.r -m p53,CD3,CD4,CD68,CD8,FOXP3,CD20,H2AXp,CKPAN -a ../annotation/supplData_withimages.csv -p ../coordinates/;
```

2. Determining and annotating cell neighbourhoods from the spe objects

* Determining density in the neighbourhood (DIN) matrices 

```bash
cd /home/Illumina/Iñaki_Sasiain/immune_spatial/analyse_clusters/DIN; 

#Create a variable with the comma separated paths to spe objects
spe_paths=$(find ../../spe_objects/* | tr "\n" ";");

#Generating DIN matrices for all the samples. Using 75 pixels as the radius
nohup Rscript ../../scripts/DIN_calculator.r -c 33 -o ${spe_paths::-1} -r 100 -n r100_DIN_sIHC;
```

* Running cluster detection and annotation

```bash
cd /home/Illumina/Iñaki_Sasiain/immune_spatial/analyse_clusters/detected_clusters; 
# Determining and annotating clusters
nohup Rscript ../../scripts/clustering.r -d ../DIN/r100_DIN_sIHC.rds -m p53,CD3,CD20,CD8,CD4 -a H2AXp,CKPAN,CD68,FOXP3 -c 33 -t p53 -n r100_DIN_sIHC_clusters;
```

3. Analysing results. Generating plots

> This section was performed using the analyse_sIHC_clusters.r R script


#### Analysing PhenoImager mIHC images

1. Running cell segmentation and phenotyping

```bash
# Getting paths to tif files
files=$(find /media/isc/Data1/PhenoImager_cores/SCAN-B_TNBC_TMA_1A/*]_component_data.tif | tr "\n" ":");

# Running cell segmentation and phenotyping
python segmentation_and_phenotyping.py -p ${files::-1} -o ./media/isc/Data1/Processed_cores/SCAN-B_TNBC_TMA_1A/;
```

2. Obtaining cell count matrices for each core.

> This section was done using the cell_counter.r R script.


3. Comparing results with TMArQ processed data.

> This section was performed using the comparing_counts.r R script.

4. Transforming results into SPE objects

```bash
# Transforming files from every slide into spe objects. 
blocks=("SCAN-B_TNBC_TMA_1A" "SCAN-B_TNBC_TMA_2A" "SCAN-B_TNBC_TMA_3A" "SCAN-B_TNBC_TMA_4A" "SCAN-B_TNBC_TMA_5A")

for block in ${blocks[@]};
    do Rscript spe_objects_from_vectra.r -m PAN.CK,CD4,CD8,CD20,CD68,FOXP3 -i /media/isc/Data1/Processed_cores/${block}/ -o /media/isc/Data1/spe_objects/;
    done;
```

5. Clustering analysis.

```bash
# Calculating DIN matrices per core
spes=$(find ../spe_objects/* | tr "\n" ";");

Rscript ../../../scripts/DIN_calculator.r -o ${spes} -g FALSE -r 100 -n r75_all_markers -c 33;

# Clustering cores
Rscript ../../../scripts/clustering.r -t PAN-CK -d ../DIN_matrices/r100_all_markers.rds -m PAN-CK,CD4,CD8,CD20 -a CD4_FOXP3,CD8_FOXP3,Other -c 33 -n r100_all_marker_clusters;
```

6. Analysing counts vs annotations and identified clusters

> This section was performed using the compare_counts_with_annotations.r script.

7. Analysing cell-to-cell distances vs subgroups

> This section was performed using the spatial_metrics_SPIAT.r script.