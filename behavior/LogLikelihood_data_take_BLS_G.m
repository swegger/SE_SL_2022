function ll = LogLikelihood_data_take_BLS_G(s,m,G,ws,sigma,b,varargin)
%% LogLikelihood_data_take_gainNoise
%
%   ll = LogLikelihood_data_take_gainNoise(s,m,G,w,wp,sigma)
%       Computes the log likelihood of the data given the gain noise model
%       and its parameters (G, w, wp, sigma)
%
%
%%

%% Defaults
nslices = 100;
sm_default = linspace(0.001,60,nslices);
Gprime_default = linspace(0.01,1.5,nslices);

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'s')
addRequired(Parser,'m')
addRequired(Parser,'G')
addRequired(Parser,'ws')
addRequired(Parser,'sigma')
addRequired(Parser,'b')
addParameter(Parser,'sm',sm_default)
addParameter(Parser,'Gprime',Gprime_default)
addParameter(Parser,'smin',0.1)
addParameter(Parser,'smax',20)
addParameter(Parser,'dx',1)

parse(Parser,s,m,G,ws,sigma,b,varargin{:})

s = Parser.Results.s;
m = Parser.Results.m;
G = Parser.Results.G;
ws = Parser.Results.ws;
sigma = Parser.Results.sigma;
b = Parser.Results.b;
sm = Parser.Results.sm;
Gprime = Parser.Results.Gprime;
smin = Parser.Results.smin;
smax = Parser.Results.smax;
dx = Parser.Results.dx;

%% Calculate log likelihoods on G/p line
% 
% % Set up variables
% SM = repmat(sm(:),[1, length(s)]);
% S = repmat(s,[size(SM,1),1]);
% M = repmat(m,[size(SM,1),1]);
% GP = M./SM;
% Gtemp = repmat(G,[size(GP,1),1]);
% 
% % Calculate probability distributions
% pS = p_sm_take_s(SM,S,ws);
% pG = p_Gprime_take_G(GP,Gtemp,sigma);
% 
% % Perform integration
% integrand = pS.*pG;
% likelihood = sum(integrand,1);
% 
% %figure;
% %histogram(log(likelihood),linspace(-50,0,50))
% 
% % Find log likelihood across pairs of m and s
% ll = sum(log(likelihood),2);

%% Calculate log likelihoods
l = length(sm);

% Set up Simpson's nodes
w = ones(2,l);
h(1) = (sm(end)-sm(1))/l;
h(2) = (Gprime(end)-Gprime(1))/l;
w(1,2:2:l-1) = 4;
w(1,3:2:l-1) = 2;
w(2,2:2:l-1) = 4;
w(2,3:2:l-1) = 2;
w(1,:) = w(1,:)*h(1)/3;
w(2,:) = w(2,:)*h(2)/3;
W = w(1,:)'*w(2,:);
W = W(:);

% Set up variables
[SM,GP] = meshgrid(sm,Gprime);
SM = repmat(SM(:),[1, length(s)]);
GP = repmat(GP(:),[1, length(s)]);

S = repmat(s,[size(SM,1),1]);
Gtemp = repmat(G,[size(GP,1),1]);
M = repmat(m,[size(SM,1),1]);

% Calculate probability distributions
integrand = p_sm_take_s(SM,S,ws);
integrand = integrand.*p_Gprime_take_G(GP,Gtemp,sigma);
integrand = integrand.*p_m_take_Gprime_sm(M,SM,GP,0.01,ws,b,smin,smax,dx);

% Perform integration
% integrand = pS.*pG.*pM;
likelihood = W'*integrand;

% Find log likelihood across pairs of m and s
ll = sum(log(likelihood),2);

%% Functions

function out = p_Gprime_take_G(Gprime,G,sigma)
    out = 1./sqrt(2*pi*sigma.^2) .* exp( -(Gprime-G).^2./(2*sigma.^2) );
            
function out = p_m_take_Gprime_sm(m,sm,Gprime,wp,ws,b,smin,smax,dx)
    %out = double(m == f(sm,Gprime));
    out = 1./sqrt(2*pi*wp^2*BLSlocal(sm,Gprime,ws,smin,smax,dx).^2) .* ...
        exp( -(m - BLSlocal(sm,Gprime,ws,smin,smax,dx)-b).^2./(2*wp^2*BLSlocal(sm,Gprime,ws,smin,smax,dx).^2) );
    
function out = p_sm_take_s(sm,s,ws)
    out = 1./sqrt(2*pi*ws.^2*s.^2) .* exp( -(sm-s).^2./(2*ws^2*s.^2) );
    
function out = BLSlocal(sm,Gprime,ws,smin,smax,dx)
    method.type = 'quad';
    method.dx = dx;
    out = Gprime.*reshape(ScalarBayesEstimators(...
        reshape(sm,[numel(sm),1]),ws,smin,smax,'method',method),...
        size(Gprime));
        