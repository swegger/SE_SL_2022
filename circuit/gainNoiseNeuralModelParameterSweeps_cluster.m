function [ws,sigGs,Gs] = gainNoiseNeuralModelParameterSweeps_cluster(simi,varargin)
%% gainNoiseNeuralModelParameterSweeps
%
%   [ws,sigGs,Gs] = gainNoiseNeuralModelParameterSweeps()
%
%   Sweeps through a default set of parameters and measures the resulting w
%   and sigG of circuit model output for comparison to behavior.
%
%%

%% Defaults
surround_weights_default = linspace(0,0.1,9);
thresholds_default = linspace(0,0.5,9);
exponentials_default = linspace(0.01,1,9);
gainNoises_default = 0;

saveOpts_default.On = true;
saveOpts_default.location = ...
    ['model_' datestr(now,'yyyymmdd') '_0'];
saveOpts_default.Figs = false;

thetas_default = 0;
speeds_default = 4:4:20;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'simi')
addParameter(Parser,'thetas',thetas_default)
addParameter(Parser,'speeds',speeds_default)
addParameter(Parser,'surround_weights',surround_weights_default)
addParameter(Parser,'thresholds',thresholds_default)
addParameter(Parser,'exponentials',exponentials_default)
addParameter(Parser,'gainNoises',gainNoises_default)
addParameter(Parser,'decoderAlgorithm','g*log2shat')
addParameter(Parser,'N',1280)
addParameter(Parser,'saveOpts',saveOpts_default)

parse(Parser,simi,varargin{:})

simi = Parser.Results.simi;
thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
surround_weights = Parser.Results.surround_weights;
thresholds = Parser.Results.thresholds;
exponentials = Parser.Results.exponentials;
gainNoises = Parser.Results.gainNoises;
decoderAlgorithm = Parser.Results.decoderAlgorithm;
N = Parser.Results.N;
saveOpts = Parser.Results.saveOpts;

%% Set up parameter combinations
[SWs,Ts,Es,GNs] = ndgrid(surround_weights,thresholds,exponentials,gainNoises);

params = [SWs(:),Ts(:),Es(:),GNs(:)];

%% For each parameter combination, run simulation
sizeProps.minEccentricity = 1;
sizeProps.maxEccentricity = 30;
disp(['Sweep ' num2str(simi) ' of ' num2str(size(params,1))])
sizeProps.surround_weight = params(simi,1);
sizeProps.exponential = params(simi,3);
sizeProps.threshold = params(simi,2);

thetaTuning.range = [-180,180,1800];
thetaTuning.amplitudeRange = [20,200,1000];
thetaTuning.widthRange = [20,90,1000];

speedTuning.range = [-1,8,1000];
speedTuning.amplitudeRange = [1,20,1000];
speedTuning.widthRange = [0.64,2.8,1000];

% With gain noise
saveOpts.location = [saveOpts.location(1:end-1) num2str(simi)];
[~, ~, ~, ~, ~, ws(simi), sigGs(simi), Gs(simi,:)] = NeuralModel(...
    'thetas',thetas,'speeds',speeds,'decoderAlgorithm',decoderAlgorithm,...
    'theta',thetaTuning,'speed',speedTuning,'N',N,...
    'gainNoise',params(simi,4),'sizeProps',sizeProps,...
    'plotMT',false,'plotDecoding',false,'plotResults',false,...
    'saveOpts',saveOpts);
    
