# sv-path-predictor
A Snakemake pipeline to predict structural variant pathogenicity and outputs a compound score calculated from three different predictors. It uses [CADD-SV](https://github.com/kircherlab/CADD-SV), [TADA](https://github.com/jakob-he/TADA), and [StrVCTVRE](https://github.com/andrewSharo/StrVCTVRE/).

# Pipeline overview
Dag is to be inserted here.

# Getting started
1. Install Snakemake version >= 6.7.0
2. Clone the repo and `cd sv-path-predictor/`
3. Prepare the inputs: `config/manifest.tab` and `config/config.yaml`
   1. cadd-sv specific- you must have these softlinked to the directory
      1. ``ln -s /path/to/softwares/pipelines/CADD-SV/annotations/ .``
      2. ``ln -s /path/to/softwares/pipelines/CADD-SV/models/ .``
4. If the dry run `-np` looks good, start your analysis
    ```console
    snakemake -s sv_path.smk --keep-going --printshellcmds --use-envmodules --use-conda
    ```
