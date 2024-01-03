![HOGVAX Logo](HOGVAX_logo.png)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10055170.svg)](https://doi.org/10.5281/zenodo.10055170)

# Abstract
Vaccination is the key component to overcoming the global COVID-19 pandemic. Peptide vaccines present a safe and cost-efficient alternative to traditional vaccines. They are rapidly produced and adapted for new viral variants. The vaccine's efficacy essentially depends on two components: the peptides included in the vaccine and the ability of major histocompatibility complex (MHC) molecules to bind and present these peptides to cells of the immune system. Due to the high diversity of MHC alleles and their diverging specificities in binding peptides, choosing a set of peptides that maximizes population coverage is a challenging task. Further, peptide vaccines are limited in their size allowing only for a small set of peptides to be included. Thus, they might fail to immunize a large part of the human population or protect against upcoming viral variants. Here, we present HOGVAX, a combinatorial optimization approach to select peptides that maximize population coverage. Furthermore, we exploit overlaps between peptide sequences to include a large number of peptides in a limited space and thereby also cover rare MHC alleles. We model this task as a theoretical problem, which we call the *Maximal Scoring k-Superstring Problem*. Additionally, HOGVAX is able to consider haplotype frequencies to take linkage disequilibrium between MHC loci into account. Our vaccine formulations contain significantly more peptides compared to vaccine sequences built from concatenated peptides. We predicted over 98% population coverage for our vaccine candidates of MHC class I and II based on single-allele and haplotype frequencies. Moreover, we predicted high numbers of per-individual presented peptides leading to a robust immunity in the face of new virus variants.

# Please Cite
Sara C. Schulte, Alexander T. Dilthey, Gunnar W. Klau, *HOGVAX: Exploiting epitope overlaps to maximize population coverage in vaccine design with application to SARS-CoV-2*, Cell Systems, Volume 14, Issue 12, 2023, Pages 1122-1130.e3, ISSN 2405-4712, https://doi.org/10.1016/j.cels.2023.11.001.

# Execute HOGVAX
HOGVAX uses the Gurobi solver for which a license is required. Further information can be found [here](https://www.gurobi.com/academia/academic-program-and-licenses/). Please also make sure that Conda is properly installed. For more information, please follow this [link](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

## Data
HLA haplotype and allele frequencies, binding affinity predictions, and input peptides are available on [Zenodo](https://doi.org/10.5281/zenodo.10055170) in the OptiVax_Data.zip folder. Click [here](https://zenodo.org/records/10450477/files/OptiVax_Data.zip?download=1) and the download will start immediately (~2GB). Please unzip the data folder after downloading. 

## Preprocessing
The provided data is already in correct format for HOGVAX. If you use your own data, we provide you with a jupyter notebook file that you can use to convert your data into the format required by HOGVAX. You may need to modify the code for your purposes. Each step in the jupyter notebook file comes along with a short explanation.

## Execution
Open a shell in the [HOGVAX](HOGVAX/) folder and create the conda environment from the `yaml` file for proper execution of HOGVAX.

```shell
conda env create -f hogvax_env.yaml
conda activate hogvax_env
```

To execute HOGVAX, call the python script with the necessary arguments, see the [example arguments](HOGVAX/example_arguments.txt).

```shell
python hogvax.py <arguments>
```
The `--help` argument gives the following list of arguments.

```angular2html
  -h, --help            show this help message and exit
  --k K, -k K           Maximal length of vaccine sequence
  --populations POPULATIONS [POPULATIONS ...], -pop POPULATIONS [POPULATIONS ...]
                        Target population(s). Default "World"
  --peptides PEPTIDES, -pep PEPTIDES
                        Preprocessed peptide file with every peptide in a new
                        line.
  --allele-frequencies F_DATA, -af F_DATA
                        (Normalized) allele frequency file.
  --ba-threshold BA_THRESHOLD, -t BA_THRESHOLD
                        Binding affinities are converted to binary data, where 
                        everything >= BA_THRESHOLD is set to 1.
  --binding-affinities BA_MATRIX, -ba BA_MATRIX
                        Binding affinity file for input peptides and alleles.
  --required_epitopes REQUIRED_EPITOPES, -epi REQUIRED_EPITOPES
                        File of peptides you want to be present in vaccine
  --min-hits MIN_HITS, -mh MIN_HITS
                        Minimum number of hits for an allele to be covered
  --maximize-peptides   Maximize number of peptides in the vaccine in a second
                        optimization
  --embedding-length EMBEDDING_LENGTH
                        Set length of embedding if used
  --embedded-peptides EMBEDDED_PEPTIDES
                        File containing embedded peptides
  --embedded-epitope_features EMBEDDED_EPITOPE_FEATURES
                        Path to embedded epitope features
  --outdir OUTDIR, -o OUTDIR
                        Output directory
  --verbose [LOGGING_ENABLED], -v [LOGGING_ENABLED]
```
