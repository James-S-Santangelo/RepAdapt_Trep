## RepAdapt ANGSD Snakemake pipeline

### Description of repository

#### Outputs from this pipeline

This pipeline analyses 1,992 white clover (_Trifolium repens_) samples from 
26 cities around the world. From each city, we sampled between 5 and 10 plants 
(mean = 8.66) from each of 1 to 5 urban and rural populations (mean = 4.42). 
This pipeline includes the following analysis steps:

1. Estimate allele frequencies at SNPs across all 1,992 samples
2. For each SNP in (1), re-estimate allele counts and frequencies in each population
of each city independently. 
3. For population from every city, estimate pi, Waterson's theta, and Tajima's D in 
20 Kb windows across the genome.
4. For each city, estimate pairwise urban-rural Hudson's Fst in 20 Kb windows across
the genome. Done for all pairwise urban-rural population combinations. 

#### General pipeline overview

See the [Rulegraph](./workflow/rulegraph.pdf) for the dependency graph of rules
included in this workflow. This is helpful for knowing the order in which rules are
run in case you want to transform the pipeline into shell scripts rather than
using Snakemake.

Note this pipeline is identical to that on the
[master](https://github.com/RepAdapt/ANGSD_Snakemake_Pipeline/tree/master?tab=readme-ov-file)
branch, but will automatically rename files according to the RepAdapt prefix
you have chosen. The prefix simply needs to be specified in the
[config](./config/server.yaml) file.

This repository contains an example Snakemake pipeline for using ANGSD
to analyze low coverage data, defined in RepAdapt as < 3X. The pipeline
performs the following actions:

1. [sweeps](./workflow/rules/angsd_sweeps_fst_thetas.smk): Estimates Site
   Allele Frequency (SAF) likelihood files for each population. Uses these to
   estimate 2DSFS for all pairwise combinations of populations, and 1DSFS for
   each population. The 1D amd 2D SFSs are then used to estimate per-site and
   windowed thetas (e.g., Pi, Waterson's) and per-site and windowed Fst,
   respectively. All these actions are parallelized across chromosomes.
2. [GEA](./workflow/rules/angsd_gea_allele_frequencies.smk): Infers SNPs across
   all samples, and then infers base counts and allele frequencies at each of
   these SNPs in each population. These actions are parallelized across
   chromosomes. 

The pipeline uses `Conda`, `Snakemake`, and `Apptainer` (formerly
`Singularity`) for workflow management and reproducibility. All Snakefiles and
directories are well-documented, but here is a brief overview of the pipeline
and directories in this repository:

A detailed description of required input files, etc. can be found below.

#### Overview of directories

- [config](./config): Snakemake configuration files for different clusters.
- [resources](./resources): Text files used in pipeline (e.g., sample
  information, chromosomes, etc.)
- [workflow](./workflow): Main Snakemake workflow with rules, environments,
  scripts, notebooks, and cluster profiles for running the pipeline on Compute
  Canada SLURM-based clusters.


### Using the pipeline

This pipeline assumes that reads have already been aligned to a reference
genome and that you have the resulting BAM files. These can be generated using
scripts 01 to 05b from the [RepAdapt mpileup SNP calling
pipeline](https://github.com/pbattlay/RepAdapt/tree/main/snp_calling_pipeline/mpileup_pipeline).

Finally, [./config/server.yaml](./config/server.yaml) should be updated to reflect
the paths to files on your machine (e.g., reference genome, etc.).

This pipeline requires `Conda` and `Singularity`:

- A minimal installation of `Conda` (i.e., Miniconda) can be installed by
  following the instructions for your platform
  [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- Installation of `Apptainer` (formerly `Singularity`) requires Admin
  privileges, but using `Apptainer` to run pre-created containers does not.
  Installation instructions can be found
  [here](https://apptainer.org/docs/admin/main/installation.html). All
  containers used in this pipeline are avalaible in [this public
  reposity](https://cloud.sylabs.io/library/james-s-santangelo), though they
  will be automatically pulled and executed by the pipeline. 

Assuming `Conda` is installed, the this repository's `Conda` environment can be
replicated by running the following command:

`conda env create -f environment.yaml -n <env_name>`

This will create a `Conda` environment named _<env\_name>_ containing a minimal
set of dependencies required to run the pipeline (e.g., Python 3.12 and
Snakemake 8.30.0).

After activating the environment (`conda activate <env_name>`), the pipeline can
be executed from the [workflow](./workflow) directory by running a command that
looks something like:

`snakemake --configfile ../config/serverl.yaml --profile profile/ --sdm conda apptaienr -j <cores>`

Here, `<cores>` is the number of cores available for executing parallel processes. 

Note that the YAML configfile at [config](./config/server.yaml) will
likely need to be modified to accommodate the paths to files on you server.

For execution on a SLURM cluster, the pipeline can be executed by running:

`snakemake --configfile ../config/server.yaml --profile profile/ --workflow-profile slurm/ --sdm conda apptainer`

Note that the YAML configfile at [slurm](./workflow/slurm/config.yaml) will
likely need to be modified to accommodate the paths on your cluster.

