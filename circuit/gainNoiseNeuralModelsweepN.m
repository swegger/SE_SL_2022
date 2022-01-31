function gainNoiseNeuralModelsweepN(varargin)
%% gainNoiseNeuralModelsweepN
%
%   Implements neural circuit model in different ways to explore the
%   effects of different parameterizations on expression of gain noise
%   effect.
%
%%

%% Defaults
Cov_default.sigf = 0.55;
Cov_default.thetaLengthConstant = 0.4;
Cov_default.speedLengthConstant = 0.3;
Cov_default.separationLengthConstant = 0.3;
Cov_default.alpha = 0;
Cov_default.diffAlpha = 0;

saveOpts_default.On = false;
saveOpts_default.Figs = false;
saveOpts_default.locationBase = '';

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'thetas',0)
addParameter(Parser,'speeds',4:4:20)
addParameter(Parser,'mymakeaxisflg',false)
addParameter(Parser,'Cov',Cov_default)
addParameter(Parser,'saveOpts',saveOpts_default)
addParameter(Parser,'N',[20 40 80 160 320 640 1280 2560 5120 10240])

parse(Parser,varargin{:})

thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
mymakeaxisflg = Parser.Results.mymakeaxisflg;
Cov = Parser.Results.Cov;
saveOpts = Parser.Results.saveOpts;
N = Parser.Results.N;

%% g*log2shat w/ dir tuning from -180 to 180
Cov.diffAlpha = 0;
Cov.separationLengthConstant = 0.3;

thetaTuning.range = [-180,180,1800];
thetaTuning.amplitudeRange = [20,200,1000];
thetaTuning.widthRange = [20,90,1000];

speedTuning.range = [-1,8,1000];
speedTuning.amplitudeRange = [1,20,1000];
speedTuning.widthRange = [0.64,2.8,1000];

for Ni = 1:length(N)
    disp(Ni)
    % Without gain noise
    saveOpts.location = [saveOpts.locationBase 'N' num2str(N(Ni)) '_g*log2shat_gainNoiseOff_' datestr(now,'yyyymmdd')];
    NeuralModel('thetas',thetas,'speeds',speeds,'gainNoise',0,'Cov',Cov,...
        'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','g*log2shat',...
        'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts,'N',N(Ni),...
        'plotMT',false,'plotDecoding',false,'plotResults',false)
    
    % With gain noise
    saveOpts.location = [saveOpts.locationBase 'N' num2str(N(Ni)) '_g*log2shat_gainNoiseOn_' datestr(now,'yyyymmdd')];
    NeuralModel('thetas',thetas,'speeds',speeds,'gainNoise',0.4,'Cov',Cov,...
        'theta',thetaTuning,'speed',speedTuning,'decoderAlgorithm','g*log2shat',...
        'mymakeaxisflg',mymakeaxisflg,'saveOpts',saveOpts,'N',N(Ni),...
        'plotMT',false,'plotDecoding',false,'plotResults',false)

end