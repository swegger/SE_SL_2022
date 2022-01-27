function [ws, sigma, G, b, ll] = BLS_G_fitter(s,m,c,varargin)
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
nslices = 80;
sm_default = linspace(2,40,nslices);
Gprime_default = linspace(0.01,1,nslices);

options_default = optimset('Display','iter');

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'s')
addRequired(Parser,'m')
addRequired(Parser,'c')
addParameter(Parser,'sm',sm_default)
addParameter(Parser,'Gprime',Gprime_default)
addParameter(Parser,'params0',[0.1,0.1,0,1])
addParameter(Parser,'lb',[0, 0, -Inf, 0])
addParameter(Parser,'ub',[Inf, Inf, Inf, Inf])
addParameter(Parser,'options',options_default)
addParameter(Parser,'smin',0.1)
addParameter(Parser,'smax',20)
addParameter(Parser,'dx',1)

parse(Parser,s,m,c,varargin{:})

s = Parser.Results.s;
m = Parser.Results.m;
c = Parser.Results.c;
sm = Parser.Results.sm;
Gprime = Parser.Results.Gprime;
params0 = Parser.Results.params0;
lb = Parser.Results.lb;
ub = Parser.Results.ub;
options = Parser.Results.options;
smin = Parser.Results.smin;
smax = Parser.Results.smax;
dx = Parser.Results.dx;


%% Fit the data

[params, ll, EXITFLAG,OUTPUT] = fmincon(@(p)(minimizer(p,s,m,c,sm,Gprime,smin,smax,dx)),...
    params0,[],[],[],[],lb,ub,[],options);

ws = params(1);
sigma = params(2);
b = params(3);
G = params(4:end);

%% Functions

function out = minimizer(p,s,m,c,sm,Gprime,smin,smax,dx)
    gains = p(4:end);
    g = gains(c);
    out = -LogLikelihood_data_take_BLS_G(s,m,g,p(1),p(2),p(3),...
        'sm',sm,'Gprime',Gprime,'smin',smin,'smax',smax,'dx',dx);