function ll = LogLikelihood_data_take_FlexWsConG(s,m,c,w0,ws,wp,varargin)
%% LogLikelihood_data_take_gainNoise
%
%   ll = LogLikelihood_data_take_gainNoise(s,m,G,w,wp,sigma)
%       Computes the log likelihood of the data given the gain noise model
%       and its parameters (G, w, wp, sigma)
%
%
%%

%% Defaults
nslices = 200;
% sm_default = logspace(log10(0.01),log10(80),nslices);
sm_default = linspace(0.001,60,nslices);

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'s')
addRequired(Parser,'m')
addRequired(Parser,'c')
addRequired(Parser,'w0')
addRequired(Parser,'ws')
addRequired(Parser,'wp')
addParameter(Parser,'f',@simpleGainNoise)
addParameter(Parser,'sm',sm_default)

parse(Parser,s,m,c,w0,ws,wp,varargin{:})

s = Parser.Results.s;
m = Parser.Results.m;
c = Parser.Results.c;
w0 = Parser.Results.w0;
ws = Parser.Results.ws;
wp = Parser.Results.wp;
f = Parser.Results.f;
sm = Parser.Results.sm;

%% Calculate log likelihoods

l = length(sm);

% Set up Simpson's nodes
w = ones(1,l);
h = (sm(end)-sm(1))/l;
w(1,2:2:l-1) = 4;
w(1,3:2:l-1) = 2;
w(1,:) = w(1,:)*h/3;
W = w(1,:)';

% Set up variables
WS = ws(c);
G = w0^2./(w0^2 + WS.^2);
SM = repmat(sm(:),[1, length(s)]);
S = repmat(s,[size(SM,1),1]);
GP = repmat(G,[size(SM,1),1]);
WS = repmat(WS,[size(SM,1),1]);
M = repmat(m,[size(SM,1),1]);

% Calculate probability distributions
pS = p_sm_take_s(SM,S,WS);
pM = p_m_take_Gprime_sm(M,SM,GP,f,wp);

% Perform integration
integrand = pS.*pM;
likelihood = W'*integrand;

% Find log likelihood across pairs of m and s
ll = sum(log(likelihood),2);

%% Functions
    
function out = p_m_take_Gprime_sm(m,sm,Gprime,f,wp)
    out = 1./sqrt(2*pi*wp^2*f(sm,Gprime).^2) .* exp( -(m - f(sm,Gprime)).^2./(2*wp^2*f(sm,Gprime).^2) );
    
function out = p_sm_take_s(sm,s,ws)
    out = 1./sqrt(2*pi*ws.^2.*s.^2) .* exp( -(sm-s).^2./(2*ws.^2.*s.^2) );
    
function out = simpleGainNoise(sm,Gprime)
    out = Gprime.*sm;
    