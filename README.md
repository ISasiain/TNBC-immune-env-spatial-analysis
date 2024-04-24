# Spatial analysis of the tumour immune microenvironment in Triple-Negative Breast Cancer

* Author: Iñaki Sasiain Casado (inaki.sasiain_casado@med.lu.se)
* Supervisors: Johan Staaf and Suze Roostee
* Master's Programme in Bioinformatics
* BINP 51 (45 ECTS)


### WORKFLOW

### SCRIPTS

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
Rscript ../../scripts/DIN_calculator.r -c 33 -o ${spe_paths::-1} -r 75 -n r75_DIN;
```

* Running cluster detection and annotation

```bash
# Determining and annotating clusters
nohup Rscript ../../scripts/clustering.r -d ../DIN/r75_DIN.rds -m p53,CD3,CD20,CD8,CD4 -a H2AXp,CKPAN,CD68,FOXP3 -c 33 -t p53 -n r75_TMArQ_clusters;
```

3. Analysing results

* 

#### Analysing PhenoImager mIHC images

1. Running cell segmentation and phenotyping

```bash
# Getting paths to tif files
files=$(find /media/isc/Data1/PhenoImager_cores/SCAN-B_TNBC_TMA_1A/*]_component_data.tif | tr "\n" ":");

# Running cell segmentation and phenotyping
python segmentation_and_phenotyping.py -p ${files::-1} -o ./media/isc/Data1/Processed_cores/SCAN-B_TNBC_TMA_1A/;
```

2. Analysing obtained results

3. Determining and annotating cell neighborhoods

* Detecting stromal/tumour areas

* Unsupervised neighbourhood detection

* Analysing determined spatial thresholds


















2. Preliminary data analysis

* TILs and survival


* Counts and survival




3. Analysing metrics in Basal or NonBasal TNBC

* Counts, distances and densities

4. Anaysing clustering

- Generating DIN matrices for every sample

```bash
cd /home/Illumina/Iñaki_Sasiain/immune_spatial/analyse_clusters/DIN; 

#Create a variable with the comma separated paths to spe objects
spe_paths=$(find ../../spe_objects/* | tr "\n" ";");

#Generating DIN matrices for all the samples. Using 75 pixels as the radius
Rscript ../../scripts/DIN_calculator.r -c 33 -o ${spe_paths::-1} -r 75 -n r75_DIN;

# Determining and annotating clusters
rscript ../../scripts/clustering.r -d ../DIN/r75_DIN.rds -m p53,CD3,CD20,CD8,CD4 -a H2AXp,CKPAN,CD68,FOXP3 -c 30 -t p53 -n r75_clusters;






#Identifying cell clusters
cd /home/Illumina/Iñaki_Sasiain/immune_spatial/analyse_clusters/detected_clusters;
nohup Rscript ../../scripts/cluster_detector.r -d ../DIN/all_samples_DIN.rds -c 40 -m p53,CD3,CD20,CD8;

#Preliminary analysis of cell cluster density
cd ../analyse_clustering/;
Rscript ../../scripts/analyse_denisty_clusters.r -d ../DIN/all_samples_DIN.rds -l ../detected_clusters/optimal_clusters.rds -a p53,CD3,CD4,CD68,CD8,FOXP3,CD20,H2AXp,CKPAN -n ../../annotation/supplData_withimages.csv;

#Cluster annotation
Rscript ../../scripts/cluster_classificator.r -d ../DIN/all_samples_DIN.rds -l ../detected_clusters/optimal_clusters.rds -u p53,CD3,CD20,CD8 -n ../../annotation/supplData_withimages.csv -a H2AXp,CKPAN,CD4,CD68,FOXP3;




# Merging the scripts in a single one.
nohup Rscript ../../scripts/clustering.r -d ../DIN/all_samples_DIN.rds -m p53,CD3,CD20,CD8,CD4 -a H2AXp,CKPAN,CD4,CD68,FOXP3 -c 30;

# Comparing the effect of the homogenity cutoff
thresholds=(0.9 1 1.1 1.2 1.3 1.4 1.5)

for num in ${thresholds[@]}; 
    do Rscript ../../scripts/clustering.r -d ../DIN/all_samples_DIN.rds -m p53,CD3,CD20,CD8 -a H2AXp,CKPAN,CD4,CD68,FOXP3 -c 40 -M ${num} -n ${num}_cutoff.clustering;
    done;


# Analysing results
thresholds=(0.9 1 1.1 1.2 1.3 1.4 1.5)

for num in ${thresholds[@]}; 
    do Rscript ../../scripts/clustering_analysis.r -d ../DIN/all_samples_DIN.rds -c optimal_clusters.rds -s ${num}_cutoff.clustering.rds -n ../../annotation/supplData_withimages.csv -m p53,CD3,CD20,CD8 -a H2AXp,CKPAN,CD4,CD68,FOXP3 -o ${num};
    done;

```

### ANALYSING OF VECTRA VERIS MULTIPLEX IHC DATA
