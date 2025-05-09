## KZ Augur and Nextstrain Pipeline

* Updated 2023-09-20 to streamlit, see previous version for command line

#### Install
- install docker https://docs.docker.com/engine/install/ubuntu/
    - Follow postinstallation instructions
    - `sudo groupadd docker`
    - `sudo usermod -aG docker $USER`
    - Restart your terminal and type in `groups`. If docker is not there restart your computer and check again.
- install miniconda https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
- Git clone this repo
    - `cd` to where you want this folder
    - `git clone git@pandora.mriglobal.org:cprice/kz_augur.git`
- Set up your environment
    - `cd kz_augur`
    - `conda env create -f KZ_augur.yml`
    - `conda activate KZ_augur`
    - `pip install streamlit`

#### How to run
- In a terminal, move to the kz_augur folder
- Type in `streamlit run app.py`
- A browser will pop up with your streamlit app
- Add at least 3 fastq files to one of the two references (TBEV or CCHF)
- Click on "Nextstrain"
- Select either TBEV or CCHF from the dropdown
- Uncheck any strains you do not want in your nextflow analysis, if any
- Click Submit
- Nextstrain will open in a new window. The first run will take a few minutes as it builds the docker container. Just reload the page once it's done.

#### Alternative Conda Creation
The conda environment has been difficult to reconstruct. Following these directions after setting up Docker has been successful:
- conda create -c conda-forge --name kz nodejs=16
- conda activate kz
- conda install -c conda-forge -c bioconda augur
- curl -fsSL --proto '=https' https://nextstrain.org/cli/installer/linux | bash
- conda config --add channels bioconda
- conda config --add channels conda-forge
- conda install -c bioconda samtools
- conda install bcftools
- conda install minimap2
- conda install samtools
- conda install scikit-learn
- conda install sourmash
- conda install umap-learn
- conda install seaborn
- conda install bokeh
- conda install holoviews
- conda install scikit-image
- conda install matplotlib
- pip install streamlit

#### Process Description
The pipeline begins with the mapping of the fastq sequences to the selected reference, either Crimean Congo Hemorrhagic Fever (CCHF) or Tick Borne Encephalitis Virus (TBEV) using Minimap2. The output is a .SAM file, which is then converted and sorted into a .BAM file using Samtools. Then, Bcftools performs variant calling on the sorted .BAM file to create a compressed vcf.gz file. Tabix is then used to create an index file for the vcf file. The final component of this initial step is a consensus sequence generation from the VCF file using Bcftools, aligning the variants back to the reference genome to create the consensus .FASTA file.

If the reference is specified as CCHF, the method only retains the shortest “S” segment of the genome as the sole focus of analysis. This process also integrates the uploaded sequence data alongside user specified meta data into an existing metadata table, which ensures each unique entry has a unqiue identifier for future runs. 

After uploading the files for analysis, the next step is to run nexstrain. Before running the Augur/Nexstrain pipeline a multiple sequence alignment is generated using Augur to the specified reference genome. This is one of the main inputs to the Augur/Nexstrain pipeline. In the pipeline, the first step is to constructed the phylogenetic tree using IQ-TREE via the augur tree function. Next is the tree refinement step, wherein the intial tree is refined using the augur refine function. This step adjusts branch lengths and resolves potential conflicts in the metadata to ensure a most accurate depiction of evolutionary relationships. The next step is the augur ancestral function, which takes the refined tree and interpolates ancestral sequences between steps, filling in gaps of how sequences may have evolved where data is absent. Mutations are then annotated using the augur translate function, which transforms the ancestral sequences into amino acid representations to identify mutations. Finally, augur traits is called as the final function which uses metadata features to correlated phylogenetic traits with metadata features. The output is then exported to a .json format compatible with Auspice, a visualization tool. All generated data by the augur pipeline and the metadata is then able to be rednered using a nextstrain view command. 

Another method of viewing the uploaded sequence data is through running a sourmash embedding. The function starts by aggregating any input labels with the genbank labeled data for the given CCHF or TBEV input. Each sequence has a MinHash sketch created, and the MinHashes are aggregated. A similarity matrix is created using sourmash compare_all_pairs, which performs a pairwise comparison between all the aggregated sketches. The output is then visualized as an altair dynamic plot, which can be used to explore similarity between the data.

All data can be exported as a final step into a .zip archive
