function [shat, ghat] = inferShatG(e,g,sig,s,w,varargin)
%% inferShatG
%
%   [shat, ghat] = inferShatG(e,g,sig,s,w)
%       Infers the measured stimulus speed (shat) and the noisy gain (ghat)
%       based on the observed eye speed (e), mean gain (g), estimate of
%       gain noise (sig), actual stimulus speed (s), and the estimated 
%       Weber fraction for speed perception (w).
%
%       The algorithm assumes that the observed eye speed results from
%           e = (g + eta_g)(s + eta_s)
%       Where eta_g is noise drawn from a zero mean Gaussian with variance
%       sig^2 and eta_s is drawn from a zero mean Gaussian with variance
%       w^2s^2 (see Egger and Lisberger, 2020). Therefore,
%           ghat = g + eta_g,
%           shat = s + eta_s,
%       and
%           e = ghat*shat.
%
%       Algorithm makes use of the assumption that motor noise is
%       negligable to search the likelihhood function, p(ghat|g)*p(shat|s)
%       along only those values that are consistent with e = ghat*shat
%
%       By default, shat and ghat are found solving for where the
%       derivative of the log likelihood function is zero and the second
%       derivative, evaluated at shat, is less than zero. The dervivative
%       is
%           dll/dshat = shat^2/(w^2s^2) + shat^3/(w^2s) - egshat/sig^2 +
%                       e^2/sig^2
%       and the second derivative is
%           d^2ll/dshat^2 = -4shat^3/(w^2s^2) + 3shat^2/(w^2s) - eg/sig^2
%       The roots of dll/dshat can be found by Matlab's roots function. The
%       solution is the real root where d^2ll/dshat^2 < 0. Multiple
%       solutions are possible for extremely unlikely values of e; in this
%       case, the value of shat closest to e is chosen as the solution.
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'e')
addRequired(Parser,'g')
addRequired(Parser,'sig')
addRequired(Parser,'s')
addRequired(Parser,'w')
addParameter(Parser,'algorithm','roots')

parse(Parser,e,g,sig,s,w,varargin{:})

e = Parser.Results.e;
g = Parser.Results.g;
sig = Parser.Results.sig;
s = Parser.Results.s;
w = Parser.Results.w;
algorithm = Parser.Results.algorithm;

%% Infer G and shat from eye speed
switch algorithm
    case {'minsearch'}
        minfun = @(p)( -logProbG(e./p,g,sig)-logProbS(p,s,w) );
        shat = fminsearch(minfun,e);
        ghat = e./shat;

    case {'roots'}
        % Find the maxmum log likelihood based on where the derivative is
        % zero
        possibleShat = roots([-1/(w^2*s^2) 1/(w^2*s) 0 -e*g/sig^2 e.^2/sig^2]);
        for ri = 1:length(possibleShat)
            realRoot(ri) = isreal(possibleShat(ri));
            tempRoot = real(possibleShat(ri));
            maximaRoot(ri) = -4*tempRoot^3/(w^2*s^2) + 3*tempRoot^2/(w^2*s) - e*g/sig^2 < 0;
        end
        if sum(realRoot & maximaRoot) == 1
            shat = real(possibleShat(realRoot & maximaRoot));
            ghat = e./shat;
        else
            warning(['Multiple solutions possible, taking solution closest to e = ' num2str(e) '.'])
            possibleShat = possibleShat(realRoot & maximaRoot);
            shat = possibleShat(min(abs(possibleShat - e)) == possibleShat - e);
            ghat = e./shat;
        end
        
end

%% Functions

%% probG
function p_g = logProbG(ghat,g,sig)
    p_g = -log(sqrt(2*pi*sig.^2)) - (ghat-g).^2./(2*sig^2);
    
%% probS
function p_s = logProbS(sm,s,w)
    p_s = -log(sqrt(2*pi*w^2*s^2)) - (sm-s).^2./(2*w^2*s^2);