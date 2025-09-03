# Organelle Measurements

## Preparation

1. Required softwares:
    - YeaZ-GUI
        - To use without GUI, download BF model of YeaZ and move it under `organelle_measure/unet/`
    - ilastik
    - Python
        - nd2reader
        - numpy
        - pandas
        - scipy
        - scikit-image
        - scikit-learn
        - matplotlib
        - seaborn
        - plotly (optional)
2. Download the project
    ```bash
    git clone git@github.com:smukherjilab/organelle-measure.git
    cd organelle-measure
    ```
3. Install the package required by the scripts
    - YeaZ requires `tensorflow=1.15`, which requires `python=3.6`, so I installed it in a separate environment:
        ```bash
        mkvirtualenv yeaz -p 3.6
        pip install -r reqirements_yeaz.txt
        ```
    - the other parts can be installed in an environment with newer python 
        ```bash
        mkvirtualenv organelle    # if you want to run the project in an isolated environment 
        pip install -r requirements.txt
        pip install -e .
        ```
4. Start your IDE _in the proper directory_

    You can use any IDE you like. I use vscode, and it needs to be started in `organelle-measure/` folder for relative paths to work.
    ```bash
    cd path/to/organelle-measure
    code . & 
    ```

## File Structure

- `organelle_measure/`
    - reusable scripts used by notebooks and scripts.
- `notebooks/`: playground to try new methods
    - `notebook1.py`
    - `notebook2.py`
    - ... 
- `scripts/`: scripts to batch process the images. Often uses `batch_apply()` function in `organelle_measure/`
    - `postprocess_vacuole_allround1data.py`
    - ...
- `test/`: files need by `notebooks/`
- `images/`
    - `raw/`: read-only microscopy images, ignored by git and synced by another software
        - `EYrainbow_glucose_largerBF`
        - `EYrainbowWhi5Up_betaEstrodiol`
        - `EYrainbow_leucine_large`
        - `EYrainbow_1nmpp1_1st`
        - `EYrainbow_rapamycin_1stTry`
        - `EYrainbow_rapamycin_CheckBistability`
            - `camera-{before/after}_EYrainbow_{experiment}-{condition}_field-{f}.nd2`
            - ...
            - `spectral-{blue/green/yellow/red}_EYrainbow_{experiment}-{condition}_field-{f}.nd2`
            - ...
            - `unmixed-{blue/green/yellow/red}_EYrainbow_{experiment}-{condition}_field-{f}.nd2`
            - ...
    - `cell/`
        - `EYrainbow_glucose_largerBF`
        - `EYrainbowWhi5Up_betaEstrodiol`
        - `EYrainbow_leucine_large`
        - `EYrainbow_1nmpp1_1st`
        - `EYrainbow_rapamycin_1stTry`
        - `EYrainbow_rapamycin_CheckBistability`
            - `binCell_EYrainbow_{experiment}-{condition}_field-{f}.tif`
            - ...
    - `preprocessed/`
        - `EYrainbow_glucose_largerBF`
        - `EYrainbowWhi5Up_betaEstrodiol`
        - `EYrainbow_leucine_large`
        - `EYrainbow_1nmpp1_1st`
        - `EYrainbow_rapamycin_1stTry`
        - `EYrainbow_rapamycin_CheckBistability`
            - `{organelle}_EYrainbow_{experiment}-{condition}_field-{f}.tif`
            - ...
            - `probability_{organelle}_EYrainbow_{experiment}-{condition}_field-{f}.tif`
            - ...
        - `leucine-large-blue-gaussian`
            - `probability_spectral-blue_EYrainbow_leu-0_hour-3,field-{f}.h5`
            - ...
    - `labelled/`
        - `EYrainbow_glucose_largerBF`
        - `EYrainbowWhi5Up_betaEstrodiol`
        - `EYrainbow_leucine_large`
        - `EYrainbow_1nmpp1_1st`
        - `EYrainbow_rapamycin_1stTry`
        - `EYrainbow_rapamycin_CheckBistability`
            - `label-{organelle}_EYrainbow_{experiment}-{condition}_field-{f}.tiff`
- `data/`: output of the pipeline, ignored by git and synced by another software
    - `ilastik/`: ilastik projects used to segment the images
    - `results/`
        - `EYrainbow_glucose_largerBF`
        - `EYrainbowWhi5Up_betaEstrodiol`
        - `EYrainbow_leucine_large`
        - `EYrainbow_1nmpp1_1st`
        - `EYrainbow_rapamycin_1stTry`
        - `EYrainbow_rapamycin_CheckBistability`
            - `cell_EYrainbow_{experiment}-{condition}_field-{f}.csv`
            - `{organelle}_EYrainbow_{experiment}-{condition}_field-{f}.csv`
    - `figures/`
    - `spectra/`
    - ...

## Pipeline

> **General rule:** <br> run `python ./scripts/{xxx}.py` after modifying `args` in each script. 
It should be a `pandas.DataFrame`, whose columns are the keyword arguments to the batch-applied function, while each row is an input.

1. Segment, label, and register cell masks:
    1. Every experiment other than "EYrainbow_leucine_large"
        - Script: `segment_cell.py`
        - Inputs: bright field camera images of the cell boundaries.
            - `int12(Y,X)`
            - `images/raw/{Experiment}/camera_EYrainbow_{experiment}-{condition}_field-{f}.nd2`
        - Outputs: label image of different cells, registered to the FOV of the confocal detector 
            - `uint16(Y,X)`
            - `images/cell/{Experiment}/binCell_EYrainbow_{experiment}-{condition}_field-{f}.tif`
2. Preprocess organelle images:
    1. peroxisome and vacuole:
        1. If not "EYrainbow_leucine_large"
            - Script: `preprocess_blue.py`
            - Inputs: Unmixed ND2 image of peroxisomes and vacuoles, as 2 channels in the same image
                - `int12(2,Z,Y,X)`, 
                - `images/raw/{Experiment}/unmixed-blue_EYrainbow_{experiment}-{condition}_field-{f}.nd2`
            - Outputs: TIF labelled images of organelles, each organelle are saves in a separated image.
                - `uint16(Z,Y,X)`
                - `images/preprocessed/{Experiment}/{organelle}_EYrainbow_{experiment}-{condition}_field-{f}.tif`
    2. ER
        - Script: `preprocess_green.py`
        - Inputs: ND2 confocal image of ER, single channel
            - `int12(Z,Y,X)`, 
            - `images/raw/{Experiment}/spectral-green_EYrainbow_{experiment}-{condition}_field-{f}.nd2`
        - Outputs: TIF labelled image of ER, single channel
            - `uint16(Z,Y,X)` 
            - `images/preprocessed/{Experiment}/{organelle}_EYrainbow_{experiment}-{condition}_field-{f}.tif`
    3. Golgi, mitochondrion, lipid droplet
        - Script: `preprocess_yellowNred.py`
        - Inputs: ND2 raw images of golgi, unmixed ND2 image mitochondrion and lipid droplets, as 2 channels in the same image
            - `int12(Z,Y,X)`, 
            - `images/raw/{Experiment}/spectral-green_EYrainbow_{experiment}-{condition}_field-{f}.nd2`
        - Outputs: TIF labelled images of the organelles, each organelle is a singel channel z-stack image.
            - `uint16(Z,Y,X)` 
            - `images/preprocessed/{Experiment}/{organelle}_EYrainbow_{experiment}-{condition}_field-{f}.tif`
3. Segment organelle images with `ilastik`, then export to simple segmentation and probability.
    - ilastik project files can be found at `data/ilastik`
    - Outputs:
        - Probability:         `float(2,Z,Y,X)`
4. Postprocess organelle images
    1. peroxisome, Golgi, lipid droplet:
        - Script: `postprocess_globular.py`
        - Inputs:
            - organelle reference image:
                - ilastik probability image
                - `float(Z,Y,X)`
                - TIFF file type
            - cell image
                - TIFF label image
                - `uint16(Y,X)`
                - `binCell_EYrainbow_{experiment}-{condition}_field-{f}.tiff`                    
        - Outputs: 
            - TIFF label image
            - `uint16(Z,Y,X)`
            - `label-vacuole_EYrainbow_{experiment}-{condition}_field-{f}.tiff`    
    2. vacuole:
        - Script: `postprocess_vacuole.py`
        - Inputs: 
            - vacuole image:
                - ilastik probability image
                - `float(2,Z,Y,X)`
                - `probability_ER_EYrainbow_{experiment}-{condition}_field-{f}.h5`
            - cell image
                - TIFF label image
                - `uint16(Y,X)`
                - `binCell_EYrainbow_{experiment}-{condition}_field-{f}.tiff`
        - Outputs: 
            - TIFF label image
            - `uint16(Z,Y,X)`
            - `label-vacuole_EYrainbow_{experiment}-{condition}_field-{f}.tiff`    
    3. ER:
        - Script: `postprocess_ER.py`
        - Inputs: 
            - ilastik probability image
            - `float(2,Z,Y,X)`
            - `probability_ER_EYrainbow_{experiment}-{condition}_field-{f}.h5`
        - Outputs: 
            - TIFF label image
            - `uint16(Z,Y,X)`
            - `label-ER_EYrainbow_{experiment}-{condition}_field-{f}.tiff` 
    4. mitochondrion:
        - Script: `postprocess_mito.py`
        - Inputs: 
            - ilastik probability image
            - `float(2,Z,Y,X)`
            - `probability_mito_EYrainbow_{experiment}-{condition}_field-{f}.h5` 
        - Outputs: 
            - TIFF label image
            - `uint16(Z,Y,X)`
            - `label-mito_EYrainbow_{experiment}-{condition}_field-{f}.tiff` 
5. Measure cell
    - Script: `measure_cell.py`
    - Inputs: 
        - TIFF label image
        - `uint16(Y,X)`
        - Segmented Cell Image from segment_cell.py (Step 1)                   
    - Outputs: 
        - `data/results/cell_EYrainbow_{experiment}-{condition}_field-{f}.csv`
6. Measure organelle properties
    - Script: `measure_organelle.py`
    - Inputs: 
        - TIFF label image
        - `uint16(Z,Y,X)`
        - `label-mito_EYrainbow_{experiment}-{condition}_field-{f}.tiff` 
    - Outputs: 
        - `data/results/{organelle}_EYrainbow_{experiment}-{condition}_field-{f}.csv`
7. Data analysis and visualization
    - Script: `csv2figures.py`
    - Inputs: `data/results`
    - Outputs: `data/`

