# ClassiCOL version 2

<img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/isoBLAST_ClassiCOL_logos_.png" width="1000" height="350" />

## ClassiCOL version 2

What is **new in ClassiCOL2**:
- Up to 10X faster analysis time
- The removal of sequences with isoBLASTs to keratin has been scrapped. (Trypsin isoBLASTs are still excluded.)
- Indicative mixture deconvolution for multiple species samples and non-database taxa.
- New plots: Taxonomic slope plot, multiple alignment plot, order uniqueness overview plot, measured data-based tree plot 
- Updated **bird collagen database**. This database includes 209 bird species (8399 sequences). Old bird database has been removed.
- Update: ClassiCOL_2_0_1 can now import **fragpipe** results files. (Only psm.tsv are used)

## Previous updates of ClassiCOL version 1.0.0
UPDATE version 1_0_2, bugs in output file fixed
UPDATE 12th of June 2025: added additional species to the original database
    
## Code and User's Guide
Welcome to the user guide to ClassiCOL. Here we explain how to use the algorithm and how to interpret the results. If you have any additional questions please contact ian.engels@ugent.be

### Citation
When using ClassiCOL please cite:

I. Engels, A. Burnett, P. Robert, C. Pironneau, G. Abrams, R. Bouwmeester, P. Van der Plaetsen, K. Di Modica, M. Otte, L. G. Straus, V. Fischer, F. Bray, B. Mesuere, I. De Groote, D. Deforce, S. Daled, M. Dhaenens, Classification of Collagens via Peptide Ambiguation, in a Paleoproteomic LC-MS/MS-Based Taxonomic Pipeline. J. Proteome Res. 24, 1907–1925 (2025).

Ian Engels, Tristan Dedrie, Synnøve M. Saugen, Simon Van de Vyver, Thijs Vandenbroucke, Kévin Di Modíca, Jan Decher, Alice Toso, Dieter Deforce, Simon Daled, Alexandra Burnett, Grégory Abrams, Maarten Dhaenens, Palaeoproteomic deconvolution of physical and genetic collagen mixtures. bioRxiv 2026.06.17.732552; doi: https://doi.org/10.64898/2026.06.17.732552

### Installation

1. Download the code in this repository. This includes:
- The ClassiCOL python script (ClassiCOL_v2).
- The Demo folder (if you want to run the demo).
- The MISC folder (contains distance csv and the unimod database).
- The BoneDB folder, which contains the curated ClassiCOL collagen fasta files.
- Download the requirements.txt file to install all additional packages.
   **Put all these folders in the ClassiCOL_version_x_x_x folder downloaded from GitHub**

2. Download the **UniProt taxonomy database**:
- Go to https://www.uniprot.org.
- Navigate to Taxonomy or use this link adress: https://www.uniprot.org/taxonomy?query=*.
- Download the taxonomy by clicking the download button and choosing TSV as file format. You should include the **Common name, Synonyms, Other names, Scientific name, Lineage, and Rank** columns in the download.
- Not recommended but possible: you can restrict the taxonomy download file if you do not want to include all taxa.
- Put the tsv file in the **MISC** folder.

<img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/ClassiCOL-files-layout.PNG"/>

3. Open Anaconda command Prompt and navigate to the location of the folder to where you downloaded the ClassiCOL folders. This is used as the base directory
4. Install the required packages using `pip install -r requirements.txt`.

### Usage
Use the following command to start the algorithm with the demo data:
```sh
$ python ClassiCOL.py -d path_to_the_script -l path_to_folder_containing_your_search_results -s MASCOT -t Mammalia -c number_of_CPUs -b S
```
**Warning: Windows has an upperlimit for the path/directory name length, keep this in mind when running ClassiCOL!**

You can use the arguments as follows:
  - `-l` folder location containing your personal Mascot \*.csv, MaxQuant \*.txt, or Manual \*.csv output files. In case you want to test the algorithm, a MASCOT output file is provided in the Demo folder. Accessible by using `-l Demo`
```sh
$ python ClassiCOL.py -d path_to_the_script -l Demo -s MASCOT
```
  - `-s` `MASCOT`, `MaxQuant`, `PEAKS` or `Manual` (search engine parameter)
  - `-t` (optional) you can restrict the taxonomy by specifying it, e.g., Pecora or for species: Bos_taurus or both: Homo_sapiens/Canis
  - `-m` specify the fixed modification used during protein extraction, e.g., C,45.987721 or multiple with C,45.98/M,...
  - `-f` (optional) location of the folder containing a custom database in fasta format
  - `-d` the directory to where the ClassiCOL algorithm is located on your computer
  - `-c` The number of CPUs you want to use (default = 3 less than available on your computer)
  - `-b` Either use S for single bone analysis or M for Mixtures. Default = mixtures
    
### A dummy example
1. **Input files:**
    - MASCOT.csv:
     Download your results directly from MASCOT in csv format
    - MaxQuant.txt
     Use the output datafile containing peptides and PTM data from MaxQuant in txt format (e.g. evidence.txt)
    - PEAKS.csv:
     Use the peptide.csv output file
    - Fragpipe.tsv:
     Use the psm.tsv output file
    - Manual.csv:
     A manual csv can be made and used as input. This file should include a sequence and (if present) the modification/s with localisation. N-term location =0, first amino acid has location 1, and C-term uses -1 as location number e.g.:
```csv
seq,modifications
GAAGLPGPK,6|Oxidation
GFSGLDGAK,
AGPPGPPGPAGK,3|Oxidation|9|Oxidation
```

2. **Batch searches:**
   - MASCOT: Place all MASCOT csv files in the same folder. The algorithm will automatically analyse all files in this folder
   - MaxQuant: As with MASCOT, you can place all files in the same folder. Additionally if 1 output file contains multiple experiments, the algorithm will automatically recognise this and analyse each experiment individually
   - Manual: Same as MASCOT

3. **The ClassiCOL output:**
ClassiCOL will put all results in the folder 'ClassiCOL_outputs', here each experiment will get its own folder for easy access. This will contain the heatmap, sunburst plot, sunburst plot with species missingness, rescored_barplot, rescored_lineplot, temporary csv output files and the final csv output file. For batch searches there will be a summary output file in the ClassiCOL_output folder.
In version 2 there will be an additional folder called 'mixtures' which contains the new plots listed below.

5. **Interpretation of the results:**
ClassiCOL will provide an estimation of taxonomy based on the available sequences in the ClassiCOL database and peptides from your search engine. It is always up to the user to interpret what these results mean.

- **The Heatmap**:
  The heatmap shows the path the algorithm will take given the NCBI taxonomy (y axis) and how the proteins relate to each other (x axis). The colours show the number of peptides assigned to each protein after isoBLAST.
  
  <img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/heatmap_example_html.png" width="1500" height="1500" />
  
- **The sunburst**:
  This figure shows an interactive overview of the output of your ClassiCOL search. A colour scheme is used to highlight the most likely classification/s (the more yellow = the more likely). By hovering over the sunburst plot you can see the taxonomic score, the number of attributed peptides and the number of isoBLASTed peptides. You can zoom in by clicking on the sunburst plot, and zoom out by clicking on the center node (or by refreshing).

<img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/sunburst_example_html.png"/>
  
- **The sunburst with missingness**:
  This plot shows exactly the same results as the sunburst plot, however now it includes all known species in NCBI that were not present during the ClassiCOL analysis. Only branches neighboring the main branch are shown up to the Order level. e.g. attached to the Family node, every missing genus (no represenative in the database used) will be shown. All species relevant to the output which are present in NCBI taxonomy but which are not included in the database used are shown in grey. This can be used to guide interpretation in cases of extinct/hybrid/non-database species.

  <img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/sunburst_with_missingness_example_html.png" />
  
- **The temporary output csv**:
  This csv is generated after the initial classification. The species/taxa are ranked to likelihood and peptides-proteins that were used during the classification are shown.
  
- **Rescored barplot**:
  For each classification, the top results are rescored. This rescoring is based on uniqueness within the top scoring group of species, meaning that all peptides shared amongst these species will be excluded. The overlap that has uniqueness is shown in this barplot.

  <img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/barplot_example_html.png"/>
  
- **Rescored lineplot**:
  This lineplot shows how the scoring changes amongst top-scoring candidates. When a dropoff is noticed after rescoring, lower-scoring candidates can be considered as discardable. When no drop-off is noticable, the sample may be comprised of a physical and/or genetic mixture.

  <img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/lineplot_example_html.png" width="400" height="500" />
  
- **The final output csv**:
  This is an easy-to-navigate output after rescoring. This includes peptide-protein information and classification information.
  
- **The batch summary csv**:
  This is a minimal information file that gives an overview of the top results alongside some metadata from the batch search.

- **The mixture multiple sequence alignment**:
  These plots show per taxonomic 'order' the multiple sequence alignment of the selected proteins within the samples.

- **The before mixture analysis plot**:
  This plot shows the taxa that have been discarded from the mixture analysis. There needs to be enough 'order' level uniqueness in comparison to other candidates

- **Tree plot**:
  A tree plot is constructed based on distances amoung species based on measured data only. Theoretical sequences resemble potential missing taxa from the database reconstructed by the mixture algorithm.

- **Taxonomic slope plots**:
  This plot shows the remaining candidates after mixture deconvolution. Per theoretical species, the taxonomic distances are calculated to give a better understanding of the location within the taxonomic tree for a species not present in the database.

**WARNING:** The algorithm can use a substantial amount of the available CPU and memory. When not enough is free, there is a chance the algorithm will go into error. 
