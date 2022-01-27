function [ws, w0, wp, ll] = FlexWsConG_fitter(s,m,c,varargin)
%% GainNoise_fitter
%
%   [ws, sigma, wp, G, ll] = FlexConG_fitter(s,m)
%       Fits the parameters of the FlexConG model:
%           m = = G*(s+etaS) + etaM
%       with etaS distributed as a scalar Gaussian with mean zero and
%       variance ws^2*s^2, etaM distributed as a scalar Gaussian with mean
%       zero and variance wp^2*(G*s)^2, and G set as
%           G = w_0^2/(w_0^2 + w_s^2);
%
%
%%

%% Defaults
sm_default = linspace(0.001,60,100);

options_default = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'s')
addRequired(Parser,'m')
addRequired(Parser,'c')
addParameter(Parser,'f',@simpleGainNoise)
addParameter(Parser,'sm',sm_default)
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
params0 = Parser.Results.params0;
lb = Parser.Results.lb;
ub = Parser.Results.ub;
options = Parser.Results.options;


%% Fit the data

[params, ll, EXITFLAG,OUTPUT] = fmincon(@(p)(minimizer(p,s,m,c,sm)),...
    params0,[],[],[],[],lb,ub,[],options);

wp = params(1);
w0 = params(2);
ws = params(3:end);

%% Functions

function out = minimizer(p,s,m,c,sm)
    out = -LogLikelihood_data_take_FlexWsConG(s,m,c,p(2),p(3:end),p(1),...
        'sm',sm);
    
function out = simpleGainNoise(s,G)
    out = G.*s;