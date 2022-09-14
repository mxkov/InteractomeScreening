# Interactome Screening
This repository contains the source code and raw results for the interactome screening project.

### Prerequisites
- Make sure that the working directory is set to wherever this Readme is located.
- All packages used are listed in ```scripts_common/full-session-info.R```, while the versions are specified in ```sessionInfo.txt```.
- Additionally, [oncoboxlib](https://gitlab.com/oncobox/oncoboxlib) needs to be installed (version 1.2.3).

### Data
The repository already contains the pre-calculated results of the analysis, both raw and summarized. To summarize and analyze the results, all you need to do is unpack ```pathways/arr_hgnc.zip```. In order to run the entire pipeline from scratch, however, you will need to do the following:

- Download the TCGA data to ```data/``` using the manifest files from its subdirectories. All files listed in a manifest should be downloaded to the same directory where this manifest is located. Additionally, files from the .tar.gz archives from ```data/TCGA-*/Clinical/``` and ```data/TCGA-*/Biospecimen/``` should be extracted to these respective directories.
- Extract the files from ```pathways/arr_hgnc.zip``` to ```pathways/```.
- Download classical pathway databases from [OncoboxPD](https://open.oncobox.com): start a session, then switch from ```Folders``` to ```Databases``` in the menu at the top of the page. Database names and versions used in this project are specified in ```pathways/databases/versions.txt```. Place the database directories under ```pathways/databases/```.

<sup>Hopefully I haven't missed anything.</sup>

### How to Run
To proceed to analyzing the pre-calculated results, skip steps 1-2 (and possibly step 3).

1. In ```scripts_common/input-locos.R```, specify which tumor localizations you want to process. For a full list of localizations, see ```names(gs)``` from ```scripts_common/locos-info.R```.
2. Run ```scripts_common/runs/runall_1.R```, then ```scripts_common/runs/runall_2.R```. The first script will prepare normalized expression matrices, PALs, formatted clinical data, etc. for the tumor localizations specified, and the second script will run the survival analysis. The two stages were split into two scripts because the second one uses way less RAM and I usually restart RStudio between the two. To do everything at once, simply toggle the three ```run_surv*``` arguments in ```scripts_common/runs/runall_1.R```.
3. Summarize the raw results using ```scripts_common/analyze-results/summarize.R```. A file named ```summary.csv``` will appear in ```locos/```. You might want to back up the existing ```summary.csv``` first.
4. For further analysis and visualization options, see other scripts in ```scripts_common/analyze-results/```. 

### How Everything Works
- Due to the convoluted ways in which some clinical data are organized in TCGA, their structure needs to be specified manually for each project. This is done in ```locos/*/scripts/get-args.R```.
- Each localization can be processed manually using ```locos/*/scripts/manual.R```. See the comments in ```scripts_common/master.R``` to find out which arguments you may need to specify there.
- Both ```scripts_common/runs/runall_1,2.R``` and ```locos/*/scripts/manual.R``` call the ```doLoco()``` function from ```scripts_common/master.R```. ```doLoco()``` is the big function that pre-processes data and calculates PALs, checks the data, and calls survival analysis functions from ```scripts_common/analysis/analysis-surv.R```. Pretty much every step performed by ```doLoco()``` can be skipped, which allows for more flexibility (e.g. you won't have to pre-process data all over again if the downstream analysis goes wrong).
- The pre-processed data and raw results are written to ```locos/*/results_TCGA/``` by default. The raw results of my own -- the ones actually used for the report -- are already present in this repository, along with some small files needed to analyze them. Most of the other files were too large.
- ```scripts_common/analyze-results/``` contains scripts for summarizing the raw results, building plots, etc.
- ```scripts_common/analysis/analysis-diff.R``` contains functions for analyzing the performance of PALs and gene expression levels as diagnostic biomarkers, i.e. to discriminate between histological types. This kind of analysis can also be called by ```doLoco()``` for a small number of tumor localizations, but it was not included in the report.
- The scripts in ```scripts_common/analysis/tools/``` contain small handy functions that I had to use more than once.

### Known Issues
- Normal tissue controls are not yet implemented for classical PALs, as I did not get to try them out yet.

### How to Cite
- When using classical pathways from [OncoboxPD](https://open.oncobox.com), please cite:

  Marianna A. Zolotovskaia, Victor S. Tkachev, Anastasia A. Guryanova, Alexander M. Simonov, Mikhail M. Raevskiy, Victor V. Efimov, Ye Wang, Marina I. Sekacheva, Andrew V. Garazha, Nicolas M. Borisov, Denis V. Kuzmin, Maxim I. Sorokin, Anton A. Buzdin (2022). OncoboxPD: human 51 672 molecular pathways database with tools for activity calculating and visualization. *Computational and Structural Biotechnology Journal*, *20*, 2280-2291. ISSN 2001-0370, https://doi.org/10.1016/j.csbj.2022.05.006.

### Contacts
Feel free to ask questions, report issues, give feedback, etc. to kovalenko.ma@phystech.edu.

