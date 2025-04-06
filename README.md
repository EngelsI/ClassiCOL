# ClassiCOL version 1.0.2

<img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/ClassiCOL-logo5_whitebackground.png" width="1000" height="350" />

## Updates since ClassiCOL version 1.0.0
UPDATE version 1_0_2, bugs in output file fixed
UPDATE 27March2025: Distance matrix saves faster
1. Addition of general input file (for non-MASCOT/MaxQuant users)
2. User defined usage of CPUs now possible
3. Batch search remembers isoBLASTED peptides (decrease in computing time throughout batch searches)
4. Protein distance calculations are now faster, precalcutated distance file for pecora can be found in the MISC folder
5. Addition of sunburst plot, containing species not present in the ClassiCOL database
6. Addition of easy to navigate csv output file, including rescored values
7. Addition of summary output file for batch searches
8. '_Bos javanicus_','_Bubalus kerabau_','_Capricornis sumarensis_', '_Daubentonia madagascariensis_', '_Eulemur rufifrons_', '_Macaca thibetana thibetana_', '_Mustela lutreola_', '_Mustela nigripes_', '_Ovis canadensis_', '_Ovibos moschatus_', '_Petaurus breviceps papuanus_', and '_Tachyglossus aculeatus_' were added to the ClassiCOL collagen database. _Homo sapiens_ COL1A1 and COL1A2 Uniprot reference sequences were exchanged.
9. Improved speed of alignment and missingness-plot building
    
## Code and User's Guide
Welcome to the user guide to ClassiCOL. Here will be explained how to use the algorithm and how to interprete the results. If you have any additional questions please contact maarten.dhaenens@ugent.be

### Citation
When using ClassiCOL please cite:
Engels, I. et al. ClassiCOL: LC-MS/MS analysis for ancient species Classification via Collagen peptide ambiguation. bioRxiv 2024.10.01.616034 (2024) doi:10.1101/2024.10.01.616034.

### Installation

1. Download the code in this repository. This includes:
     a) The ClassiCOL python script
     b) The Demo folder (if you want to run the demo)
     c) The MISC folder (contains distance csv and the unimod database)
     d) The BoneDB folder, which contains the curated ClassiCOL collagen fasta files
     e) Download the requirements.txt file to install all additional packages
   **Put all these folders in the ClassiCOL_version_x_x_x folder downloaded from GitHub**

<img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/ClassiCOL-files-layout.PNG"/>

3. Open Anaconda command Prompt and navigate to the location of the folder to where you downloaded the ClassiCOL folders.
4. Install the required packages using `pip install -r requirements.txt`.

### Usage
Use the following command to start the algorithm with the demo data:
```sh
$ python ClassiCOL.py -d path_to_the_script -l path_to_folder_containing_your_search_results -s MASCOT -t Mammalia -c number_of_CPUs
```

You can use the arguments as follows:
  - `-l` folder location containing your personal Mascot \*.csv, MaxQuant \*.txt, or Manual \*.csv output files. In case you want to test the algorithm a MASCOT output file is provided in the Demo folder. Accessable by using `-l Demo`
```sh
$ python ClassiCOL.py -d path_to_the_script -l Demo -s MASCOT
```
  - `-s MASCOT`, `MaxQuant` or `Manual` (specify the search engine used)
  - `-t` (optional) you can restrict the taxonomy by specifying it, e.g., Pecora or for species: Bos_taurus or both: Homo_sapiens/Canis
  - `-m` specify the fixed modification used during protein extraction, e.g., C,45.987721 or multiple with C,45.98/M,...
  - `-f` (optional) location of the folder containing a custom database in fasta format
  - `-d` the directory to where the ClassiCOL algorithm is located on your computer
  - `-c` The number of CPUs you want to use (default = 3 less than available on your computer)
### A dummy example
1. **Input files:**
   - MASCOT.csv:
     Download your results directly from MASCOT in csv format
   - MaxQuant.txt
     Use the output datafile containing peptides and locational data from MaxQuant in txt format
   - Manual.csv:
     A manual csv can be made and used as input. This file should include a sequence and (if present) the modification/s with localisation. N-term location =0, first amino acid has location 1, and C-term uses -1 as location number e.g.:
```csv
seq,modifications
GAAGLPGPK,6|Oxidation
GFSGLDGAK,
AGPPGPPGPAGK,3|Oxidation|9|Oxidation
```

2. **Batch searches:**
   - MASCOT: Place all MASCOT csv files in the same folder. The algoithm will automatically analyse all files in this folder
   - MaxQuant: As with MASCOT, you can place all files in the same folder. Additionally if 1 output file contains multiple experiments, the algorithm will automatically recognise this and analyse each experiment individually
   - Manual: Same as MASCOT

3. **The ClassiCOL output:**
ClassiCOL will put all results in the folder 'ClassiCOL_outputs', here each experiment will get its own folder for easy access. This will contain the heatmap, sunburst plot, sunburst plot with species missingness, rescored_barplot, rescored_lineplot, temporary csv output files and the final csv output file. For batch searches there will be a summary output file in the ClassiCOL_output folder.

4. **Interpretation of the results:**
ClassiCOL will provide an estimation of taxonomy based on the available sequences in the ClassiCOL database and peptides from your search engine. It is always up to the user to interpret what these results mean!

- **The Heatmap**:
  The heatmap shows the path the algorithm will take given the NCBI taxonomy (y axis) and how the proteins relate to each other (x axis). The colours show abundance of peptides assigned to each protein after isoBLAST.
  
  <img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/heatmap_example_html.png" width="1500" height="1500" />
  
- **The sunburst**:
  This figure shows an interactive overview of the output of your ClassiCOL search. A colour scheme is used to highlight the most likely classification/s (the more yellow = the more likely). By hovering over the sunburst plot you can see the number of attributed peptides and the number of isoBLASTed peptides. You can zoom in by clicking on the sunburst plot, and zoom out by clicking on the center node (or by refreshing).

<img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/sunburst_example_html.png"/>
  
- **The sunburst with missingness**:
  This plot shows exactly the same results as the sunburst plot, however now it includes all known species in NCBI that were not present during the ClassiCOL analysis. Only branches neighboring the main branch are shown up to the Order level. e.g. attached to the Family node, every missing genus (no represenative in the database used) will be shown.

  <img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/sunburst_with_missingness_example_html.png" />
  
- **The temporary output csv**:
  This csv is generated after the initial classification. The species/taxa are ranked to likelihood and peptides-proteins are shown that were used during the classification.
  
- **Rescored barplot**:
  For each classification, the top result is rescored. This rescoring is based on uniqueness within the top scoring group of species, meaning that all peptides shared amongst these species will be excluded. The overlap that has uniqueness is shown in this barplot.

  <img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/barplot_example_html.png"/>
  
- **Rescored lineplot**:
  This lineplot shows how the scoring changes amongst top scoring candidates. When a dropoff is noticed after rescoring, lower-scoring candidates can be considered as discardable. When no drop-off is noticable, the sample can be comprised of a physical and/or genetic mixture.

  <img src="https://github.com/EngelsI/ClassiCOL/blob/main/240405_tarandus_1_1_p/lineplot_example_html.png" width="400" height="500" />
  
- **The final output csv**:
  This is an easy-to-navigate output after rescoring. This includes peptide-protein information and classification information.
  
- **The batch summary csv**:
  This is a minimal information file that gives an overview of the top results alongside some metadata from the batch search.


**WARNING_1:** Depending on the number of unique peptides in your sample and the number of species you want to consider, the isoBLAST calculations could take a while (about 2min for +/- 1000 unique peptides per species). An overnight search is recommended. Batch searches will go much quicker towards the end.

**WARNING_2:** The algorithm can use a substantial amount of the available CPU and memory. When not enough is free, there is a chance the algorithm will go into error.
