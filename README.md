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
    `TestGainNoiseFitting('sigma',SIGMA)`
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
    `modelparams.variant='BLS';
    modelparams.method.type='quad';
    modelparams.method.dx=1;
    TestGainNoiseFitting('modelparams',modelparams,'sigma',SIGMA,'gains',[0.3 0.6,0.9])`
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
    `CircuitN_vs_behavioralCorrelation('baseName','_g*log2shat_gainNoiseOff_20210602.mat','baseNameNoise','_g*log2shat_gainNoiseOn_20210602.mat')`

### Supplementary Table 1
1. 1. Use the script "GainNoiseModelBehavior.m"
2. For monkey R
  1. In the code under "Reggie" (line 1), change set `saveTable` to `true`
  2. Run the code under "Reggie" to set up analysis parameters (line 1)
  3. Run code under "Main analysis" to plot results of behavioral analysis (line 17)
  4. Results are saved as a LaTeX table in "R_tableYYYYMMDD.tex" in current working directory with YYYY being the year, MM being the month, and DD being the day of the month code was executed.
3. For monkey X
  1. In the code under "Xtra" (line 64), change set `saveTable` to `true`
  2. Run the code under "Xtra" (line 64)
  3. Run the code under "Main analysis" (line 76)
  4. Results are saved as a LaTeX table in "X_tableYYYYMMDD.tex" in current working directory with YYYY being the year, MM being the month, and DD being the day of the month code was executed.

## General use
### Behavioral analysis

### Behavioral modeling

### Biomimetic circuit modeling
