# spectracular-pub
Selected stable versions of spectracular

## Setup of Python environment.
You can do it in three ways:

1. Use the `requirements.txt` file to build/modify your own env.

2. Use the following one-liner to create a seperate conda environment.
`conda create -n spectracular -c conda-forge -c astropy 'tqdm>=4.31' 'specutils==0.2.2' 'scipy>=1.1' 'python>=3.6' 'numpy>=1.15' 'matplotlib>=3.0.2' 'astropy>=3.1.1' 'pyyaml>=3.13'`

3. Use the `environment.yml` file with conda to skip the resolving step.
`conda env create -f environment.yml`

## Get additional data – Coelho models – and free space.
If you wish to use Spectracular without its fitting feature, skip this step.

We require you to download and unpack a large archive of high-res models:
`http://specmodels.iag.usp.br/fits_search/compress/s_coelho14_highres.tgz`

Spectracular also requires another ~5 GB to process these models in its own way.

## Fill in the spectracular/config.py file.

Provide full paths (ended with `/` character) to:
* a directory containing the unpacked `s_coelho14_highres/` catalog;
* a directory where Spectracular can put its processing results (another ~5 GB);
* the provided SVO filter profiles;
* additional SVO filter profiles, if you wish to use other than provided.

## Try both examples!
Or one of them, if you skipped the step no. 2.

First, run the `example-XSHOOTER-V.py` script, which uses Spectracular's core code stucture (without its fitting features).
Next, try `example-modelling.py`. The first run might take additional 5–10 minutes to complete (and take these ~5 GB of your free space mentioned above).

## Write your own scripts, provide other data.
Have fun.
