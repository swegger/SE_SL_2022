# SE_SL_2022
Instructions for use of code to implement data analysis and circuit results as published in Seth W. Egger and Stephen G. Lisberger. "Neural structure of a sensory decoder for motor control." *Nature Communications* (2022)

We have included our full set of analysis, behavioral modeling, and circuit simulation tools for the user. As such, we provide some general use instructions. However, our main objective here is to provide the user with code to recreate the results from Egger and Lisberger (2022), and the associated scripts/functions are illustrative for running the code in general.

To recreate the results, the user needs the follwing data sets, unzipped and in the current MATLAB directory:
1. Behavioral
  a) Reggie_MultiSizePursuit.mat
  b) Xtra_MultiSizePursuit.mat
2. Biomimetic circuit
  a) parameterSweep.zip
  b) circuitN.zip

System requirements: MATLAB R2020b
Tested on: Ubuntu 20.04; OSX 10.15.7

## Recreation of results
### Figure 1
1. Use the script "GainNoiseModelBehavior.m"
2. For monkey R
  1. Run the code under "Reggie" to set up analysis parameters (line 1)
  2. Run code under "Main analysis" to plot results of behavioral analysis (line 17)
3. For monkey X
  1. Run the code under "Xtra" (line 64)
  2. Run the code under "Main analysis" (line 76)
4. For bootstrap analysis, change "Ncv" to the number of desired bootstraps. (Warning, bootstrapping is implemented serially. Large values of Ncv will take a long time to complete)

### Figure 2
1. Use the script "GainNoiseModelBehavior" as for Figure 1.

### Figure 3
For behavioral results:
1. Using the script "GainNoiseModelBehavior.m," run code for monkeys R and X as in Figure 1
2. Run code under "Gain noise w vs standard model w" (line 120)

For theoretical results:
1. At the MATLAB command line, run

    ```
    TestGainNoiseFitting('sigma',SIGMA)
    ```

with SIGMA being the standard deviation of the noise in the gain process (e.g. 0, 0.1, or 0.5 for the left, middle and right graphs of panel b, respectively).

### Figure 4

### Figure 5

### Supplementary Figure 1
1. Use the script "GainNoiseModelBehavior.m"
2. For monkey R
  1. Run the code under "Reggie" to set up analysis parameters (line 1)
  2. Run code under "Main analysis w/ eccentricity threshold" to plot results of behavioral analysis (line 25)
3. For monkey X
  1. Run the code under "Xtra" (line 64)
  2. Run the code under "Main analysis w/ eccentricity threshold" (line 84)

### Supplementary Figure 2
1. Use the script "GainNoiseModelBehavior.m"
2. For monkey R
  1. Run the code under "Reggie" to set up analysis parameters (line 1)
  2. Run code under "Perform core analysis in multiple analysis time windows" (line 33)
3. For monkey X
  1. Run the code under "Xtra" (line 64)
  2. Run the code under "Perform core analysis in multiple analysis time windows" (line 92)

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
1. Use the script "GainNoiseModelBehavior.m"
2. For monkey R
  1. Run the code under "Reggie" to set up analysis parameters (line 1)
  2. Run code under "Main analysis" to plot results of behavioral analysis (line 17)
3. For monkey X
  1. Run the code under "Xtra" (line 64)
  2. Run the code under "Main analysis" (line 76)

### Supplementary Figure 5

### Supplementary Figure 6
1. At the MATLAB command line, run

    ```
    CircuitN_vs_behavioralCorrelation('baseName','_g*log2shat_gainNoiseOff_20210602.mat','baseNameNoise','_g*log2shat_gainNoiseOn_20210602.mat')
    ```

### Supplementary Table 1
1. 1. Use the script "GainNoiseModelBehavior.m"
2. For monkey R
  1. In the code under "Reggie" (line 1), change set the option 'saveTable' to `saveTable` defined as

      ```
      saveTable.On = true
      saveTable.directory = DESIRED_SAVE_LOCATION      
      ```      
  2. Run the code under "Reggie" to set up analysis parameters (line 1)
  3. Run code under "Main analysis" to plot results of behavioral analysis (line 17)
  4. Results are saved as a LaTeX table in "R_tableYYYYMMDD.tex" in DESIRED_SAVE_LOCATION with YYYY being the year, MM being the month, and DD being the day of the month code was executed.
3. For monkey X
  1. In the code under "Xtra" (line 64), change set `saveTable` to

      ```
      saveTable.On = true
      saveTable.directory = DESIRED_SAVE_LOCATION
      ```

  2. Run the code under "Xtra" (line 64)
  3. Run the code under "Main analysis" (line 76)
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

*Note, if directory/data/sname1/MultiSizePursuit.mat exists, SetupSmoothPursuitProject will load that data and then check if more data needs to be extracted.

*Analysis of data*
Data and analysis and model fitting are all performed by the function gainNoiseMultiSize. A variety of options are available to the user. Please see the execution of analysis in GainNoiseModelBehavior for example inputs to gainNoiseMultiSize

### Behavioral modeling
*Simulation of candidate behavioral models*
As part of our exploration of simple behavioral models that might be able to explain the observed data, we developed a tool for simulating behavior. To verify our model fits, we also developed methods for testing the fit of candidate models to the simulated behavior. The function TestGainNoiseFitting does most of the heavy lifting. Below we will give a couple examples of how to use the function.

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
