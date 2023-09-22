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
    - `git clone https://github.com/mriglobal/kz_analysis.git`
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
