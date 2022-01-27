function ll = LogLikelihood_data_take_gainNoise(s,m,G,ws,wp,sigma,varargin)
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
addRequired(Parser,'wp')
addRequired(Parser,'sigma')
addParameter(Parser,'f',@simpleGainNoise)
addParameter(Parser,'sm',sm_default)
addParameter(Parser,'Gprime',Gprime_default)

parse(Parser,s,m,G,ws,wp,sigma,varargin{:})

s = Parser.Results.s;
m = Parser.Results.m;
G = Parser.Results.G;
ws = Parser.Results.ws;
wp = Parser.Results.wp;
sigma = Parser.Results.sigma;
f = Parser.Results.f;
sm = Parser.Results.sm;
Gprime = Parser.Results.Gprime;

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
pS = p_sm_take_s(SM,S,ws);
pG = p_Gprime_take_G(GP,Gtemp,sigma);
pM = p_m_take_Gprime_sm(M,SM,GP,f,wp);

% Perform integration
integrand = pS.*pG.*pM;
likelihood = W'*integrand;

% Find log likelihood across pairs of m and s
ll = sum(log(likelihood),2);

%% Calculate log likelihoods by focusing on density near intersection of s, m, and G
% for triali = 1:length(s)
%     centroid = ([m(triali)/s(triali); s(triali)] + [G(triali); m(triali)/G(triali)])/2;
% %     ghat = (m(triali)/s(triali) + G(triali))/2;
%     if centroid(2)-5*(wp+ws)*centroid(2) < 0.001
%         sm = linspace(0.001,...
%             centroid(2)+5*(wp+ws)*centroid(2), 200);
%     else
%         sm = linspace(centroid(2)-5*(wp+ws)*centroid(2),...
%             centroid(2)+5*(wp+ws)*centroid(2), 200);
%     end
%     Gprime = linspace(centroid(1)-5*sigma, centroid(1)+5*sigma, 200);
%     
%     % Set up Simpson's nodes
%     l = length(sm);
%     w = ones(2,l);
%     h(1) = (sm(end)-sm(1))/l;
%     h(2) = (Gprime(end)-Gprime(1))/l;
%     w(1,2:2:l-1) = 4;
%     w(1,3:2:l-1) = 2;
%     w(2,2:2:l-1) = 4;
%     w(2,3:2:l-1) = 2;
%     w(1,:) = w(1,:)*h(1)/3;
%     w(2,:) = w(2,:)*h(2)/3;
%     W = w(1,:)'*w(2,:);
%     W = W(:);
%     
%     % Set up variables
%     [SM,GP] = meshgrid(sm,Gprime);
%     SM = SM(:);
%     GP = GP(:);
%     
%     S = repmat(s(triali),[size(SM,1),1]);
%     Gtemp = repmat(G(triali),[size(GP,1),1]);
%     M = repmat(m(triali),[size(SM,1),1]);
%     
%     % Calculate probability distributions
%     pS = p_sm_take_s(SM,S,ws);
%     pG = p_Gprime_take_G(GP,Gtemp,sigma);
%     pM = p_m_take_Gprime_sm(M,SM,GP,f,wp);
%     
%     % Perform integration
%     integrand = pS.*pG.*pM;
%     likelihood(triali) = W'*integrand;
% end
% ll = sum(log(likelihood),2);


%% Functions

function out = p_Gprime_take_G(Gprime,G,sigma)
    out = 1./sqrt(2*pi*sigma.^2) .* exp( -(Gprime-G).^2./(2*sigma.^2) );
    
function out = p_m_take_Gprime_sm(m,sm,Gprime,f,wp)
    out = 1./sqrt(2*pi*wp^2*f(sm,Gprime).^2) .* exp( -(m - f(sm,Gprime)).^2./(2*wp^2*f(sm,Gprime).^2) );
    
function out = p_sm_take_s(sm,s,ws)
    out = 1./sqrt(2*pi*ws.^2*s.^2) .* exp( -(sm-s).^2./(2*ws^2*s.^2) );
    
function out = simpleGainNoise(sm,Gprime)
    out = Gprime.*sm;
    