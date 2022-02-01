# SE_SL_2022
Instructions for use of code to implement data analysis and circuit results as published in Seth W. Egger and Stephen G. Lisberger. "Neural structure of a sensory decoder for motor control." *Nature Communications* (2022)

We have included our full set of analysis, behavioral modeling, and circuit simulation tools for the user. As such, we provide some general use instructions. However, our main objective here is to provide the user with code to recreate the results from Egger and Lisberger (2022), and the associated scripts/functions are illustrative for running the code in general.

To recreate the results, the user should download the following data sets from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5889167.svg)](https://doi.org/10.5281/zenodo.5889167)
, and place the unzipped data and in the current MATLAB path:
1. Behavioral*
    1. Reggie_MultiSizePursuit.mat
    2. Xtra_MultiSizePursuit.mat

2. Biomimetic circuit**
    1. circuitN.zip

*.mat files contain trial-by-trial data extracted from the original data files. .mat files contain all the data necessary to run behavioral analysis. Raw data files are as recorded on the day of each experiment. Both are contained in monkeyR_data.zip and monkeyX_data.zip.

System requirements: MATLAB R2020b
Tested on: Ubuntu 20.04; OSX 10.15.7

## Recreation of results
### Figure 1
1. Use the script `GainNoiseModelBehavior`
2. For monkey R
    1. Run the code under "Reggie" to set up analysis parameters (line 1)
    2. Run code under "Main analysis" to plot results of behavioral analysis (line 13)
3. For monkey X
    1. Run the code under "Xtra" (line 56)
    2. Run the code under "Main analysis" (line 68)
4. For bootstrap analysis, change "Ncv" to the number of desired bootstraps. (Warning, bootstrapping is implemented serially. Large values of Ncv will take a long time to complete)

### Figure 2
1. Use the script `GainNoiseModelBehavior` as for Figure 1.

### Figure 3
For behavioral results:
1. Using the script `GainNoiseModelBehavior`, run code for monkeys R and X as in Figure 1
2. Run code under "Gain noise w vs standard model w" (line 110)

For theoretical results:
1. At the MATLAB command line, run

```
TestGainNoiseFitting('sigma',SIGMA)
```

with SIGMA being the standard deviation of the noise in the gain process (e.g. 0, 0.1, or 0.5 for the left, middle and right graphs of panel b, respectively).

### Figure 4
The circuit as instantiated in Figure 4 can be plot by running the following from the MATLAB command line:

```
plotTuningProperties('N1280_g*log2shat_gainNoiseOn_20210602.mat')
```

### Figure 5
The circuit as instatiated in panels b and e can be plot by running the following from the command line:

```
plotTuningProperties('N1280_g*log2shat_gainNoiseOff_20210602.mat')
```
The circuit as instatiated in panels c and f can be plot by:

```
plotTuningProperties('N1280_g*log2shat_gainNoiseOn_20210602.mat')
```

### Supplementary Figure 1
1. Use the script `GainNoiseModelBehavior`
2. For monkey R
    1. Run the code under "Reggie" to set up analysis parameters (line 1)
    2. Run code under "Main analysis w/ eccentricity threshold" to plot results of behavioral analysis (line 21)

3. For monkey X
    1. Run the code under "Xtra" (line 56)
    2. Run the code under "Main analysis w/ eccentricity threshold" (line 76)

### Supplementary Figure 2
1. Use the script `GainNoiseModelBehavior`
2. For monkey R
    1. Run the code under "Reggie" to set up analysis parameters (line 1)
    2. Run code under "Perform core analysis in multiple analysis time windows" (line 29)

3. For monkey X
    1. Run the code under "Xtra" (line 56)
    2. Run the code under "Perform core analysis in multiple analysis time windows" (line 84)

### Supplementary Figure 3
1. At the MATLAB command line, run

    ```
    modelparams.variant='BLS';
    modelparams.smin = 0.1;
    modelparams.smax = 20;
    modelparams.b = 0;
    modelparams.method.type='quad';
    modelparams.method.dx=1;
    TestGainNoiseFitting('modelparams',modelparams,'sigma',SIGMA,'gains',[0.3 0.6,0.9])
    ```

  with SIGMA being the standard deviation of the noise in the gain process (e.g. 0.01 or 0 for panels c and d, respectively).

### Supplementary Figure 4
1. Use the script `GainNoiseModelBehavior`
2. For monkey R
    1. Run the code under "Reggie" to set up analysis parameters (line 1)
    2. Run code under "Main analysis" to plot results of behavioral analysis (line 13)

3. For monkey X
    1. Run the code under "Xtra" (line 56)
    2. Run the code under "Main analysis" (line 68)

### Supplementary Figure 5
Due to the combined size of the 729 different instantiations of the biomimetic circuit used in this figure, it is not possible for us to share the data associated with the original results. To recreate the results requires iterating the model for each combination of parameter values. We have provided a wrapper function that will set up the list of parameterizations used and run one iteration of parameters selected from this list. Here, we provide the necessary code to run one such iteration. It is *highly* recommended that the user deploys the code on a computational cluster such that each iteration can be run in parallel.

1. Set up the save options at the MATLAB command line using:

    ```
    saveOpts.On = true;
    saveOpts.location = [DESTINATION_DIRECTORY NAME DATE ID '_0'];
    saveOpts.Figs = false;
    ```

where `DESTINATION_DIRECTORY` is the desired directory where the results of each iteration of the circuit will be saved, `NAME` specifies a global file name, and `ID` is unique identifier (ID) that can later be used to recover the data associated with this parameter sweep (we recommend inserting the date of creation using `ID = datestr(now,'yyyymmdd')`).

2. Pass the save options to the wrapper function that will build the desired list and execute instance `LISTNUM` of that list using

    ```
    gainNoiseNeuralModelParameterSweeps_cluster(LISTNUM,'decoderAlgorithm','g*log2shat','N',1280,'saveOpts',saveOpts);
    ```

This will execute the biomimetic circuit as in the main paper, but with parameters taken from the `LISTNUM`th entry of the list of 729 different parameter combinations. See the General Use section below for details on how to specify the parameterization list. The output will save a file at `saveOpts.location` with `_0` replaced by `_LISTNUM`. To generate all the results, one must run this code for each parameterization of the list (e.g. 1 through 729).

3. Analyze the results using

    ```
    gainNoiseNeuralModelParameterSweep_analysis(PARAMETER_SWEEP_DIR,'calcNeuronPursuitCorrelation',true,'dataDate',ID,'gainNoise',0.4)
    ```

where `PARAMETER_SWEEP_DIR` is the directory that contains the results of the parameter sweep analysis from step 2 and `ID` is the unique identifier specified in step 1.

### Supplementary Figure 6
1. At the MATLAB command line, run

    ```
    CircuitN_vs_behavioralCorrelation('baseName','_g*log2shat_gainNoiseOff_20210602.mat','baseNameNoise','_g*log2shat_gainNoiseOn_20210602.mat')
    ```

### Supplementary Table 1
1. 1. Use the script `GainNoiseModelBehavior`
2. For monkey R
    1. In the code under "Reggie" (line 1), change set the option 'saveTable' to `saveTable` defined as

        ```
        saveTable.On = true
        saveTable.directory = DESIRED_SAVE_LOCATION      
        ```      
      2. Run the code under "Reggie" to set up analysis parameters (line 1)
      3. Run code under "Main analysis" to plot results of behavioral analysis (line 13)
      4. Results are saved as a LaTeX table in "R_tableYYYYMMDD.tex" in DESIRED_SAVE_LOCATION with YYYY being the year, MM being the month, and DD being the day of the month code was executed.

3. For monkey X
    1. In the code under "Xtra" (line 64), change set `saveTable` to

        ```
        saveTable.On = true
        saveTable.directory = DESIRED_SAVE_LOCATION
        ```

    2. Run the code under "Xtra" (line 56)
    3. Run the code under "Main analysis" (line 68)
    4. Results are saved as a LaTeX table in "X_tableYYYYMMDD.tex" in DESIRED_SAVE_LOCATION with YYYY being the year, MM being the month, and DD being the day of the month code was executed.

## General use
### Behavioral analysis
*Extraction of raw data*

Also provided is the raw data as collected using Maestro (Lisberger Lab software; https://sites.google.com/a/srscicomp.com/maestro/home). To extract the data from these files, open MATLAB and run the following code:

  ```
  d = SetupSmoothPursuitProject(sname,'MultiSizePursuit',directory)
  ```

with `sname` corresponding to the name of the monkey that one wishes to recover the data for and `directory` corresponding to that contains all the relevant data and code.

*Note, the code expects the project to have the following directory structure:

```
directory
└───analysis
    └───extraction
└───data
    └───sname1
        |   sname1_MultiSizePursuit.mat
        └───raw_data_folders
    └───sname2
        |   sname2_MultiSizePursuit.mat
        └───raw_data_folders

```

*Note, if directory/data/sname/MultiSizePursuit.mat exists, `SetupSmoothPursuitProject` will load that data and then check if more data needs to be extracted. Data collected on a given day should be contained in one `raw_data_folder` and data collected on different days should be sorted into different folders.

*Analysis of data*

Data and analysis and model fitting are all performed by the function `gainNoiseMultiSize`. A variety of options are available to the user. Please see the execution of analysis in `GainNoiseModelBehavior` for example inputs to `gainNoiseMultiSize`

### Behavioral modeling
*Simulation of candidate behavioral models*

As part of our exploration of simple behavioral models that might be able to explain the observed data, we developed a tool for simulating behavior. To verify our model fits, we also developed methods for testing the fit of candidate models to the simulated behavior. The function `TestGainNoiseFitting` does most of the heavy lifting. Below we will give a couple examples of how to use the function.

To simulate the simple gain noise model, from the command line in MATLAB simply type:

  ```
  TestGainNoiseFitting('sigma',SIGMA,'gains',GAINS)
  ```

This will simulate the model with the mean gain levels specified in GAINS and the standard deviation of gain noise specified by SIGMA.

To simulate a different model, one must specify the model parameters. See instructions for reproducing Supplementary Figure 3 for an example using the BLS estimator. See the code itself for other possible model simulations.

To fit a model to the simulated data, one must specify the fit parameters. For example, to simulate the BLS model and fit the BLS model to the data one should specify the model parameters as:

  ```
  modelparams.variant='BLS';
  modelparams.method.type = 'quad';
  modelparams.method.dx=1;
  modelparams.smin = 0.1;
  modelparams.smax = 20;
  modelparams.b = 0;
  ```

and the fit parameters as:

  ```
  fitparams.variant = 'BLS';
  fitparams.smin = 0.1;
  fitparams.smax = 20;
  ```

and then run:

  ```
  TestGainNoiseFitting('modelparams',modelparams,'gains',GAINS,'sigma',SIGMA,'fitparams',fitparams);
  ```

*Note: the code is designed to allow for the simulated model to be different than the fit model for testing of model identifiability.

### Biomimetic circuit modeling
The code for recreating the biomimetic circuit results from the paper used above load results from runs of the model previously. Below we provide brief instructions to implement the circuit, mostly by way of demonstrating how to implement the circuit used in the main text of Egger and Lisberger (2022).

*Circuit basics*

The basic idea is to generate a population of MT neurons and their responses to a given set of target speeds and a target size. `DirSpeedSizeLocMT` is the function which generates the MT population response. The population response is then decoded by the putative downstream circuit. `DecodeMT` is the function that performs decoding. A full 'experiment' with different speeds and target sizes is wrapped in the function `NeuralModel`.

To set up an experiment, the first step is to define the tuning properties of the population. To set up the range of preferred directions, direction tuning amplitudes, and direction tuning widths, set:

  ```
  thetaTuning.range = [-180,180,1800];
  thetaTuning.amplitudeRange = [20,200,1000];
  thetaTuning.widthRange = [20,90,1000];
  ```

Each trio specifies the `[MIN,MAX,NUMBER]` of values to sample from. For example, `thetaTuning.range = [-180,180,1800]` creates 1800 different preferred directions, linearly spaced between -180 and 180 deg. Similarly, to set up the range of preferred log_2(speeds), speed tuning amplitudes, and speed tuning widths, set:

  ```
  speedTuning.range = [-1,8,1000];
  speedTuning.amplitudeRange = [1,20,1000];
  speedTuning.widthRange = [0.64,2.8,1000];
  ```

Then supply the ranges to `NeuralModel` like this:

  ```
  NeuralModel('thetas',0,'speeds',4:4:20,'gainNoise',0.4,...
          'theta',thetaTuning,'speed',speedTuning,'N',1280)
  ```

This will simulate the response of a population of N = 1280 MT neurons to each of the 3 sizes used in the paper (2, 6, and 20 deg), moving in the 0 deg direction at 4, 8, 12, 16, and 20 deg/s.

*Iterate neuron number*

In the paper, we present analysis of the effect of the number of MT neurons on the ciruit's behavior. We have provided a wrapper function which will do so. Simply enter the following into the MATLAB command line:

  ```
  saveOpts.On = true;
  saveOpts.Figs = false;
  saveOpts.locationBase = 'DESIRED_SAVE_LOCATION';
  gainNoiseNeuralModelsweepN('saveOpts',saveOpts)
  ```

This is simulate populations of 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, and 10240 MT neurons, as implemented in the main paper, with and without gain noise. The output of each number of neurons and gain noise condition with be saved to the directory `DESIRED_SAVE_LOCATION` (see code for details). The corresponding results can then be analyzed using `CircuitN_vs_behavioralCorrelation` as for Supplementary Figure 6 above, with the correct base file names.

*Iterate MT neuron parameterization*

In the paper, we also iterate 729 different parameterizations of the threshold nonlinearities and surround properties of the model neurons. Because the computational cost is prohibitive on a single computer, we recommend the user implements the code in a parallel fashion on a computational cluster. We provide a wrapper function that will set up a grid of parameter realizations and then run one instance from this list, instance `simi`.

  ```
  [ws,sigGs,Gs] = gainNoiseNeuralModelParameterSweeps_cluster(simi,'surround_weights',SURROUND_WEIGHTS,'thresholds',THRESHOLDS,'exponentials',EXPONENTIALS)
  ```

This will run the circuit model as in the main paper, but with parameters chosen from the `simi`th entry of a list made up of every possible combination of parameters in `SURROUND_WEIGHTS`, `THRESHOLDS`, and `EXPONENTIALS`. After the user has run this for every possible `simi`, the results of every parameterization can then be analyzed as for Supplementary Figure 5 using `gainNoiseNeuralModelParameterSweep_analysis`.
