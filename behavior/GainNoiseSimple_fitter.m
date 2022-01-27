function [ws, sigma, G, ll] = GainNoiseSimple_fitter(s,m,c,varargin)
%% GainNoise_fitter
%
%   [ws, sigma, wp, G, ll] = GainNoise_fitter(s,m)
%       Fits the parameters of the full Gain Noise model:
%           m = = (G+etaG)*(s+etaS) + etaM
%       with etaS distributed as a scalar Gaussian with mean zero and
%       variance ws^2*s^2, etaM distributed as a scalar Gaussian with mean
%       zero and variance wp^2*(G*s)^2, and etaG distributed as a Gaussian
%       with mean zero and variance sigma^2
%
%
%%

%% Defaults
nslices = 100;
sm_default = linspace(0.001,60,nslices);
Gprime_default = linspace(0.01,1.5,nslices);
%sm_default = linspace(0.01,80,5000);

options_default = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'s')
addRequired(Parser,'m')
addRequired(Parser,'c')
addParameter(Parser,'f',@simpleGainNoise)
addParameter(Parser,'sm',sm_default)
addParameter(Parser,'Gprime',Gprime_default)
addParameter(Parser,'params0',[0.1,0.1,1])
addParameter(Parser,'lb',[0, 0, 0])
addParameter(Parser,'ub',[Inf, Inf, Inf])
addParameter(Parser,'options',options_default)

parse(Parser,s,m,c,varargin{:})

s = Parser.Results.s;
m = Parser.Results.m;
c = Parser.Results.c;
f = Parser.Results.f;
sm = Parser.Results.sm;
Gprime = Parser.Results.Gprime;
params0 = Parser.Results.params0;
lb = Parser.Results.lb;
ub = Parser.Results.ub;
options = Parser.Results.options;


%% Fit the data

[params, ll, EXITFLAG,OUTPUT] = fmincon(@(p)(minimizer(p,s,m,c,sm,Gprime)),...
    params0,[],[],[],[],lb,ub,[],options);

ws = params(1);
sigma = params(2);
G = params(3:end);

%% Functions

function out = minimizer(p,s,m,c,sm,Gprime)
    gains = p(3:end);
    g = gains(c);
    out = -LogLikelihood_data_take_gainNoiseSimple(s,m,g,p(1),p(2),...
        'sm',sm,'Gprime',Gprime);
    
function out = simpleGainNoise(s,G)
    out = G.*s;