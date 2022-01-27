function [e, likelihood] = ScalarBayesEstimators(m,wm,xmin,xmax,varargin)
%% ScalarBayesEstimators
%
%   e = ScalarBayesEstimators(m,wm,xmin,xmax)
%       Computes the BLS estimates (e) of x, given a set of measurements
%       (m). Assumes uniform prior over x between xmin and xmax
%
%   e = ScalarBayesEstimators(m,wm,xmin,xmax,'method',method_opts)
%       Computes the BLS estimate using the the method specifiied by
%       method_opts.
%           method_opts.type = 'integral'
%               Use integral functions, default
%           method_opts.type = 'trapz'
%               Approximates the integral using trapz using the step
%               size method_opts.dx
%           method_opts.type = 'quad'
%               Approximates the integral using Simpson's quadrature using
%               the step size method_opts.dx
%
%   e = ScalarBayesEstimators(...,'estimator',estimator_opts)
%       Computes the estimate specified by estimator_opts.type
%           estimator_opts.type = 'BLS'
%               Bayes least-squares; default.
%           estimator_opts.type = 'ObsAct'
%               Observer-Actor model. User must supply the weber fraction
%               on production in the estimator_opts.wp
%
%   e = ScalarBayesEstimators(m,wm,xmin,xmax,'prior',prior_opts)
%       Computes the BLS estimate using a prior specified by the structure
%       prior_opts.
%           TODO: prior types
%
%   e = ScalarBayesEstimators(...,'estimator',estimator_opts)
%       with
%           estimator_opts.type = 'weightedMean'
%       and
%           estimator_opts.weights = weights
%       Computs the BLS estimate using a weighted average of the
%       measurements according a set of weights.
%           TODO: integral and trapz support
%
%   e = ScalarBayesEstimators(...,'estimator',estimator_opts)
%       with
%           estimator_opts.type = 'sequential'
%       and
%           estimator_opts.transitionFunction = function handle for
%           transition probabilities
%       and
%           estimator_opts.p = vector of parameters for transition function
%       Computs the estimate using a sequential message passing
%       algorithm. In brief, updates the joint probability distribution of
%       the measurement and sample according to
%           p(x(t),m(1),m(2),...,m(t)) = p(m(t)|x(t))...
%               *int/(p(x(t-1),m(1),m(2),...,m(t-1))*p(x(t)|x(t-1))dx(t-1))
%       Algorithm is equivalent to the foward-backward message passing
%       algorithm (see Bishop, 2006; chapter 13), without the backward
%       passes. Estimate is then calculated as the expected value of x(t). 
%           TODO: integral and trapz support
%
%   NOTE: Methods requiring mmx.mex are mostly obsolete.
%
%%

%% Defaults
method_opts.type = 'integral';
estimator_opts.type = 'BLS';
prior_opts.type = 'uniform';

%% Parse inputs
inputs = inputParser;
addRequired(inputs,'m');
addRequired(inputs,'wm');
addRequired(inputs,'xmin');
addRequired(inputs,'xmax');
addParameter(inputs,'method',method_opts,@isstruct);
addParameter(inputs,'estimator',estimator_opts,@isstruct);
addParameter(inputs,'prior',prior_opts,@isstruct);

parse(inputs,m,wm,xmin,xmax,varargin{:})

m = inputs.Results.m;
wm = inputs.Results.wm;
xmin = inputs.Results.xmin;
xmax = inputs.Results.xmax;
method = inputs.Results.method;
estimator = inputs.Results.estimator;
prior = inputs.Results.prior;       %% TODO generalize to aribitrary prior


if isfield(estimator,'wy')
    wy = estimator.wy;
else
    wy = 0;
end

if ~isfield(estimator,'ObsAct')
    estimator.ObsAct = 0;
end

%% Compute the estimate
switch estimator.type
    case {'BLS','BLSbiasedLapse','BLS_wm_wp_sigp'}
        switch method.type
            case 'integral'
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'trapz'
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3));
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4)/(1+wy.^2)^estimator.ObsAct;
                e = permute(e,[2 1]);
                
            case 'quad'
                switch prior.type
                    case 'uniform'
                        % Number of measurements
                        N = size(m,2);
                        
                        % Create x-vector
                        dx = method.dx;
                        x = xmin:dx:xmax;
                        
                        % Create Simpson's nodes
                        l = length(x);
                        h = (xmax - xmin)/l;
                        w = ones(1,l);
                        w(2:2:l-1) = 4;
                        w(3:2:l-1) = 2;
                        w = w*h/3;
                        
                        % Reshape measurements for processing
                        M = permute(m,[2 3 1]);
                        M = repmat(M,[1,1,1,l]);
                        x = reshape(x,[1 1 1 l]);
                        X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                        
                        % Generate estimate
                        w = reshape(w,[1 1 1 l]);
                        w = repmat(w,[1 1 size(m,1) 1]);
                        likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                        %                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                        e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                        e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                        
                    case 'Gaussian'
                        % Parameterize prior
                        mu = prior.mu;
                        sig = prior.sig;
                        
                        % Number of measurements
                        N = size(m,2);
                        
                        % Create x-vector
                        dx = method.dx;
                        x = mu-3*sig:dx:mu+3*sig;
                        
                        % Create Simpson's nodes
                        l = length(x);
                        h = (xmax - xmin)/l;
                        w = ones(1,l);
                        w(2:2:l-1) = 4;
                        w(3:2:l-1) = 2;
                        w = w*h/3;
                        
                        % Reshape measurements for processing
                        M = permute(m,[2 3 1]);
                        M = repmat(M,[1,1,1,l]);
                        x = reshape(x,[1 1 1 l]);
                        X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                        P = 1/sqrt(2*pi*sig^2) * exp( -(X(1,:,:,:)-mu).^2/(2*sig^2));
                        
                        % Generate estimate
                        w = reshape(w,[1 1 1 l]);
                        w = repmat(w,[1 1 size(m,1) 1]);
                        likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                        %                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                        e = sum(w.*X(1,:,:,:).*likelihood.*P,4)./sum(w.*likelihood.*P,4);
                        e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                       
                    otherwise
                        error(['Prior ' prior.type ' not recognized!'])
                end
                
            case 'MonteCarlo'
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
                
            case 'MonteCarlo_batch'
                % Set up integration variables
                options.N = method.N;
                options.batch_sz = method.batch_sz;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
        end
        
    case {'ObsAct','ObsActLapse'}
%         wy = estimator.wy;
        switch method.type
            case 'integral'
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N)/(1+wy.^2);
                end
                
            case 'trapz'
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3));
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4);
                e = permute(e,[2 1])/(1+wy.^2);
                
            case 'quad'
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                
                % Create Simpson's nodes
                l = length(x);
                h = (xmax - xmin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = repmat(M,[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(m,1) 1]);
                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                e = permute(e,[3 2 1])/(1+wy.^2);
                
            case 'MonteCarlo'
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2);
                
            case 'MonteCarlo_batch'
                % Set up integration variables
                options.N = method.N;
                options.batch_sz = method.batch_sz;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2);
        end
    
    case 'MAP'
        switch method.type
            case 'fminsearch'
                for i = 1:size(m,1)
                    post = @(x)(posteriorMAP(x,m(i,:),wm,xmin,xmax));
                    e(i) = fminsearch(post,(xmax-xmin)/2+xmin)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'analytical'
                N = size(m,2);
                e = mean(m,2) .* (-1 + sqrt(1 + 4*wm^2 .* mean(m.^2,2)./mean(m,2).^2)) / (2*wm^2);
%                 if N == 1
%                     e = m.*(-1+sqrt(1+4*wm.^2))./(2*wm.^2);
%                 elseif N == 2
%                     e = sum(m,2)./(4*wm^2) - sqrt( (sum(m,2)./(4*wm^2)).^2 - sum(m.^2,2)/(2*wm^2) );
%                 else
%                     error('N > 2 for analytical solution not yet supported!')
%                 end
                e(e < xmin) = xmin;
                e(e > xmax) = xmax;
                e = e/(1+wy.^2)^estimator.ObsAct;
                
            case 'integral'
                error('Not yet supported!')
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'trapz'
                error('Not yet supported!')
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3));
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4);
                e = permute(e,[2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'quad'
                error('Not yet supported!')
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                
                % Create Simpson's nodes
                l = length(x);
                h = (xmax - xmin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = repmat(M,[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(m,1) 1]);
                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'MonteCarlo'
                error('Not yet supported!')
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
                
            case 'MonteCarlo_batch'
                error('Not yet supported!')
                % Set up integration variables
                options.N = method.N;
                options.batch_sz = method.batch_sz;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
        end
        
    case 'MLE'
        switch method.type
            case 'fminsearch'
                for i = 1:size(m,1)
                    post = @(x)(posteriorMLE(x,m(i,:),wm,xmin,xmax));
                    e(i) = fminsearch(post,(xmax-xmin)/2+xmin)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'analytical'
                N = size(m,2);
                e = mean(m,2) .* (-1 + sqrt(1 + 4*wm^2 .* mean(m.^2,2)./mean(m,2).^2)) / (2*wm^2)/(1+wy.^2)^estimator.ObsAct;
                
            case 'integral'
                error('Not yet supported!')
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'trapz'
                error('Not yet supported!')
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3));
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4);
                e = permute(e,[2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'quad'
                error('Not yet supported!')
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                
                % Create Simpson's nodes
                l = length(x);
                h = (xmax - xmin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = repmat(M,[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(m,1) 1]);
                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'MonteCarlo'
                error('Not yet supported!')
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
                
            case 'MonteCarlo_batch'
                error('Not yet supported!')
                % Set up integration variables
                options.N = method.N;
                options.batch_sz = method.batch_sz;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo_batch','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
        end
        
    case 'aveEstimates'
        weights = estimator.weights;
       % TODO
       switch method.type
            case 'integral'
                error('not yet supported!')
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'trapz'
                error('not yet supported!')
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin - round(wm*xmin):dx:xmax + round(wm*xmax);
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3))/(1+wy.^2)^estimator.ObsAct;
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4);
                e = permute(e,[2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'quad'
                for i = 1:size(m,2)
                    mtemp = m(:,i);
                    N = size(mtemp,2);
                    
                    % Create x-vector
                    dx = method.dx;
                    x = xmin:dx:xmax;
                    
                    % Create Simpson'€™s nodes
                    l = length(x);
                    h = (xmax - xmin)/l;
                    w = ones(1,l);
                    w(2:2:l-1) = 4;
                    w(3:2:l-1) = 2;
                    w = w*h/3;
                    
                    % Reshape measurements for processing
                    M = permute(mtemp,[2 3 1]);
                    M = repmat(M,[1,1,1,l]);
                    x = reshape(x,[1 1 1 l]);
                    X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                    
                    % Generate estimate
                    w = reshape(w,[1 1 1 l]);
                    w = repmat(w,[1 1 size(mtemp,1) 1]);
                    likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                    %                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                    e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                    e = permute(e,[3 2 1]);
                    E(:,i) = e/(1+wy.^2)^estimator.ObsAct;
                end
                e = sum(repmat(weights,size(m,1),1).*E,2);
                
           case 'MonteCarlo'
               error('Not yet supported!')
                % Number of measurements
                m = sum(repmat(weights,size(m,1),1).*m,2);
                N = size(m,2);
                
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
       end
    
    
    case 'aveMeasurements'
        weights = ones(1,size(m,2))/size(m,2);
       % TODO
       switch method.type
            case 'integral'
                error('not yet supported!')
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'trapz'
                error('not yet supported!')
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin - round(wm*xmin):dx:xmax + round(wm*xmax);
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3));
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4);
                e = permute(e,[2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'quad'
                % Number of measurements
                m = sum(repmat(weights,size(m,1),1).*m,2);
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                
                % Create Simpson'€™s nodes
                l = length(x);
                h = (xmax - xmin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = repmat(M,[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(m,1) 1]);
                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                
           case 'MonteCarlo'
               error('Not yet supported!')
                % Number of measurements
                m = sum(repmat(weights,size(m,1),1).*m,2);
                N = size(m,2);
                
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
       end
       
    case 'weightedMean'
        weights = estimator.weights;
       % TODO
       switch method.type
            case 'integral'
                error('not yet supported!')
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N)/(1+wy.^2)^estimator.ObsAct;
                end
                
            case 'trapz'
                error('not yet supported!')
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin - round(wm*xmin):dx:xmax + round(wm*xmax);
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3));
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4);
                e = permute(e,[2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'quad'
                % Number of measurements
                m = sum(repmat(weights,size(m,1),1).*m,2);
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                
                % Create Simpson'€™s nodes
                l = length(x);
                h = (xmax - xmin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = repmat(M,[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(m,1) 1]);
                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                
           case 'MonteCarlo'
               error('Not yet supported!')
                % Number of measurements
                m = sum(repmat(weights,size(m,1),1).*m,2);
                N = size(m,2);
                
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
       end

    case 'weightedMeanOptimal'
        weights = estimator.weights;
       % TODO
       switch method.type
            case 'integral'
                error('not yet supported!')
                N = size(m,2);
                fBLS = @(m,xmin,xmax,wm,N)(integral(@(x)(x.*(1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true)./integral(@(x)((1./(sqrt(2*pi)*wm*x)).^N .* exp( -(m-x)'*(m-x) ./ (2*wm.^2.*x.^2) )),xmin,xmax,'ArrayValued',true));
                for i = 1:size(m,1)
                    e(i) = fBLS(m(i,:),xmin,xmax,wm,N);
                end
                
            case 'trapz'
                error('not yet supported!')
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin - round(wm*xmin):dx:xmax + round(wm*xmax);
                x = reshape(x,[1 1 1 length(x)]);
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = reshape(m,1,1,1,length(x));
                x = repmat(x,size(m,1),1,size(m,3));
                
                % Generate estimate
                likelihood = (1./(sqrt(2*pi)*wm*x)).^N .* exp( -(mmx('mult',permute(x-M,[2 1 3 4]),x-M))./(2*wm.^2.*x.^2) );
                e = trapz(x.*likelihood,4)./trapz(likelihood,4);
                e = permute(e,[2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'quad'
                % Number of measurements
                m = sum(repmat(weights,size(m,1),1).*m,2);
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                
                % Create Simpson'€™s nodes
                l = length(x);
                h = (xmax - xmin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                % Reshape measurements for processing
                M = permute(m,[2 3 1]);
                M = repmat(M,[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(m,1) 1]);
                likelihood = ( (1./sqrt(2*pi*sum(weights.^2))/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2).*sum(weights.^2) ) );
%                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                
           case 'MonteCarlo'
               error('Not yet supported!')
                % Number of measurements
                m = sum(repmat(weights,size(m,1),1).*m,2);
                N = size(m,2);
                
                % Set up integration variables
                options.N = method.N;
                numeratorFun = @(x,m)(MonteCarloIntegrand_numerator(x,m,wm));       % Numerator of BLS funciton
                denominatorFun = @(x,m)(MonteCarloIntegrand_denominator(x,m,wm));    % Denominator of BLS function
                
                % Find the numerator and denominator of the BLS function
                numerator = ndintegrate(numeratorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                denominator = ndintegrate(denominatorFun,[xmin xmax],'method','MonteCarlo','options',options,'ExtraVariables',m);
                
                e = numerator./denominator/(1+wy.^2)^estimator.ObsAct;
       end
       
    case 'sequential'
        switch method.type
            case 'integral'
                % TODO
                error('Integral method not yet supported for sequential estimator')
                
            case 'quad'
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin-3*wm*xmin:dx:xmax+3*wm*xmax;
                
                % Generate matrices for all combinations of x and m
                l = length(x);
                M = repmat(m(:,1)',[l 1]);
                X = repmat(reshape(x,numel(x),1),[1 size(m,1)]);
                
                % Prior
                alpha = zeros(l,size(m,1));
                alpha(x >= xmin-1 & x <= xmax,:) = 1/(xmax-xmin);
%                alpha = repmat(1/(xmax-xmin),[l size(m,1)]);
                
                % Generate transition matrix
                [X1, X2] = meshgrid(x);
                deterministicModel = estimator.deterministicModel;
                probabilisticModel = estimator.probabilisticModel;
                boundry = estimator.boundry;
                T = TransitionFunction2(X1,X2,deterministicModel,probabilisticModel,boundry);
                
                % Create Simpson's nodes
                l = length(x);
                h = (xmax+3*wm*xmax - xmin-3*wm*xmin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                W = repmat(w,[l 1]);
                %W = repmat(permute(w,[1 3 2]),[l size(m,1)]);
                
                % Initialize distribution
                alpha = alpha.*(1./(sqrt(2*pi)*wm*X)) .* exp( -(M-X).^2 ./ (2*wm.^2.*X.^2) );       % joint probability of the first measurement and each value of x
                A(:,:,1) = alpha;
                 
                % Perform update of joint probability for each step, using the appropriate transition model
                for i = 2:N
                    M = repmat(m(:,i)',[l 1]);
                    for j = 1:size(m,1)
                        alpha_(:,j) = sum(repmat(alpha(:,j)',[size(T,1) 1]).*T.*W,2);
                    end
                    %alpha_ = alpha;
                    A_(:,:,i) = alpha_;
                    %alpha_ = permute(sum(mmx('mult',mmx('mult',repmat(permute(alpha,[1 3 2]),[1 size(T,2) 1]),T),W),1),[2 3 1]);
                    %alpha_ = sum(W.*repmat(permute(alpha(:,:),[3 2 1]),[size(T,3) 1 1]).*T,3);          % propogated joint probability before measurement update
                    alpha = (1./(sqrt(2*pi)*wm*X)) .* exp( -(M-X).^2 ./ (2*wm.^2.*X.^2) ) .* alpha_; % propogated joint probability after measurement update
                    A(:,:,i) = alpha;
                    
                    % Find the normalization constant
                    Z = repmat(sum(repmat(w',[1 size(m,1)]).*alpha),[l 1]);
                end
                
                % Find the normalization constant
                Z = repmat(sum(repmat(w',[1 size(m,1)]).*alpha),[l 1]);
%                Z = sum(repmat(w',[1 size(m,1)]).*alpha);
                
                % Determine estimate
%                e = permute(sum(repmat(w',[1 size(m,1)]).*X.*alpha)./Z,[2 1]);
                e = permute(sum(repmat(w',[1 size(m,1)]).*X.*alpha./Z),[2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'trapz'
                % TODO
                error('trapizoidal method not yet supported for sequential estimator')
                
            case 'MonteCarlo'
                % TODO
                error('Monte Carlo method not yet supported for sequential estimator')
                
        end
        
    case 'NoisyMem'
        switch method.type
            case 'integral'
                % TODO
                
            case 'trapz'
                % TODO
                
            case 'quad'
                % Number of measurements
                N = size(m,2);
                
                % Create x-vector
                dx = method.dx;
                x = xmin:dx:xmax;
                
                % Create Simpson'€™s nodes
                l = length(x);
                h = (xmax - xmin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                % Reshape measurements for processing
                M = repmat(permute(m(:,1),[2 3 1]),[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(m,1) 1]);
                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .* exp( -((X-M).^2)./(2*wm.^2.*X(1,:,:,:).^2) ) );
                for i = 2:N
                    M = repmat(permute(m(:,2),[2 3 1]),[1,1,1,l]);
                    likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .* exp( -((X-M).^2)./(2*wm.^2.*X(1,:,:,:).^2) ) );
                end
                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                
            case 'MonteCarlo'
                % TODO
        end
        
    case 'SubOptMemBias'
        switch method.type
            case 'integral'
                % TODO
            case 'trapz'
                % TODO
                
            case 'analytical'
                % TODO
                
            case 'quad'
                % Generate MLE estimate
                N = size(m,2);
                wm_drift = estimator.wm_drift;
                w_int = estimator.w_int;
                
                if N == 1
                    m = mean(m,2) .* ...
                        (-1 + sqrt(1 + 4*wm^2 .* mean(m.^2,2)./mean(m,2).^2)) ...
                        / (2*wm^2);
                    
                    % Create x-vector
                    dx = method.dx;
                    x = xmin:dx:xmax;
                    
                    % Create Simpson'€™s nodes
                    l = length(x);
                    h = (xmax - xmin)/l;
                    w = ones(1,l);
                    w(2:2:l-1) = 4;
                    w(3:2:l-1) = 2;
                    w = w*h/3;
                    
                    % Reshape measurements for processing
                    M = permute(m,[2 3 1]);
                    M = repmat(M,[1,1,1,l]);
                    x = reshape(x,[1 1 1 l]);
                    X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                    
                    % Generate estimate
                    w = reshape(w,[1 1 1 l]);
                    w = repmat(w,[1 1 size(m,1) 1]);
                    likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .*...
                        exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                    e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                    e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                    
                elseif N == 2
                    a1 = wm_drift^2/(wm_drift.^2+wm.^2);
                    a2 = wm^2/(wm_drift.^2+wm.^2);
                    w12 = wm_drift*wm/sqrt(wm_drift.^2+wm.^2);             % A "better" formulation of w12 is to multiply this result by sqrt(2) so that it represents the "effective" wm for one measurement
                    weights = [a1 a2];
                    mbar = sum(repmat(weights,size(m,1),1).*m,2)/2;        % Division by 2 is right, but a "better" formulation is to multiply w12 by sqrt(2) so w12 represents the "effective" wm
                    m2bar = sum(repmat(weights,size(m,1),1).*m.^2,2)/2;    % Division by 2 is right, but a "better" formulation is to multiply w12 by sqrt(2) so w12 represents the "effective" wm
                    m = mbar .* ...
                        (-1 + sqrt(1 + 4*w12.^2 .* m2bar./mbar.^2))...
                        ./ (2*w12.^2);
                
                    
                    % Create x-vector
                    dx = method.dx;
                    x = xmin:dx:xmax;
                    
                    % Create Simpson'€™s nodes
                    l = length(x);
                    h = (xmax - xmin)/l;
                    w = ones(1,l);
                    w(2:2:l-1) = 4;
                    w(3:2:l-1) = 2;
                    w = w*h/3;
                    
                    % Reshape measurements for processing
                    M = permute(m,[2 3 1]);
                    M = repmat(M,[1,1,1,l]);
                    x = reshape(x,[1 1 1 l]);
                    X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                    
                    % Generate estimate
                    w = reshape(w,[1 1 1 l]);
                    w = repmat(w,[1 1 size(m,1) 1]);
                    likelihood = ( (1./sqrt(2*pi)/w_int/X(1,:,:,:)) .*...
                        exp( -(sum((X-M).^2,1))./(2*w_int.^2.*X(1,:,:,:).^2) ) );
                    
                    % For some very small values of w_int, the weighted sum
                    % is zero to numerical precision; replacing with
                    % realmin in the denominator fixes this, but returns
                    % bad estimates for extremely unlikely measurement
                    % combinations (e.g. m2 >> m1)
                    denominator = sum(w.*likelihood,4);
                    if any(denominator(:) == 0)
                        warning(['For some very small values of w_int, the weighted sum' ...
                            'is zero to numerical precision; replacing with '...
                            'realmin in the denominator fixes this, but returns '...
                            'bad estimates for extremely unlikely measurement '...
                            'combinations (e.g. m2 >> m1)'])
                        denominator(denominator == 0) = realmin;
                    end
                    e = sum(w.*X(1,:,:,:).*likelihood,4)./denominator;
                    e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                    
                else
                    error('N > 2 not supported for supOptMemBias model!')
                end
                
            case 'MonteCarlo'
                % TODO
                
            case 'MonteCarlo_batch'
                % TODO
        end
                
    case 'NestedModel'
        switch method.type
            case 'integral'
                % TODO
            case 'trapz'
                % TODO
                
            case 'analytical'
                % TODO
                
            case 'quad'
                % Generate MLE estimate
                N = size(m,2);
                wm_drift = estimator.wm_drift;
                w_int = estimator.w_int;
                alpha = estimator.alpha;
                
                if N == 1
                    m = mean(m,2) .* ...
                        (-1 + sqrt(1 + 4*wm^2 .* mean(m.^2,2)./mean(m,2).^2)) ...
                        / (2*wm^2);
                    
                    % Create x-vector
                    dx = method.dx;
                    x = xmin:dx:xmax;
                    
                    % Create Simpson'€™s nodes
                    l = length(x);
                    h = (xmax - xmin)/l;
                    w = ones(1,l);
                    w(2:2:l-1) = 4;
                    w(3:2:l-1) = 2;
                    w = w*h/3;
                    
                    % Reshape measurements for processing
                    M = permute(m,[2 3 1]);
                    M = repmat(M,[1,1,1,l]);
                    x = reshape(x,[1 1 1 l]);
                    X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                    
                    % Generate estimate
                    w = reshape(w,[1 1 1 l]);
                    w = repmat(w,[1 1 size(m,1) 1]);
                    likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .*...
                        exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                    
                    % For some very small values of w_int, the weighted sum
                    % is zero to numerical precision; replacing with
                    % realmin in the denominator fixes this, but returns
                    % bad estimates for extremely unlikely measurement
                    % combinations (e.g. m2 >> m1)
                    denominator = sum(w.*likelihood,4);
                    if any(denominator(:) == 0)
                        warning(['For some very small values of w_int, the weighted sum' ...
                            'is zero to numerical precision; replacing with '...
                            'realmin in the denominator fixes this, but returns '...
                            'bad estimates for extremely unlikely measurement '...
                            'combinations (e.g. m2 >> m1)'])
                        denominator(denominator == 0) = realmin;
                    end
                    e = sum(w.*X(1,:,:,:).*likelihood,4)./denominator;
                    e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                    
                elseif N == 2
                    a1 = wm^2/(wm_drift.^2+wm.^2);
                    a2 = wm_drift^2/(wm_drift.^2+wm.^2);
                    w12 = wm_drift*wm/sqrt(wm_drift.^2+wm.^2) * sqrt(2);
                    weights = [a1 a2];
                    mbar = sum(repmat(weights,size(m,1),1).*m,2);
                    m2bar = sum(repmat(weights,size(m,1),1).*m.^2,2);
                    m = mbar .* ...
                        ((-1 + sqrt(1 + 4*w12.^2 .* m2bar./mbar.^2))...
                        ./ (2*w12.^2)).^alpha;
                    
                    
                    % Create x-vector
                    dx = method.dx;
                    x = xmin:dx:xmax;
                    
                    % Create Simpson'€™s nodes
                    l = length(x);
                    h = (xmax - xmin)/l;
                    w = ones(1,l);
                    w(2:2:l-1) = 4;
                    w(3:2:l-1) = 2;
                    w = w*h/3;
                    
                    % Reshape measurements for processing
                    M = permute(m,[2 3 1]);
                    M = repmat(M,[1,1,1,l]);
                    x = reshape(x,[1 1 1 l]);
                    X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                    
                    % Generate estimate
                    w = reshape(w,[1 1 1 l]);
                    w = repmat(w,[1 1 size(m,1) 1]);
                    likelihood = ( (1./sqrt(2*pi)/w_int/X(1,:,:,:)) .*...
                        exp( -(sum((X-M).^2,1))./(2*w_int.^2.*X(1,:,:,:).^2) ) );
                    e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                    e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                else
                    error('N > 2 not supported for supOptMemBias model!')
                end
                
            case 'MonteCarlo'
                % TODO
                
            case 'MonteCarlo_batch'
                % TODO
        end
        
    case {'BLS_wm1wm2'}
        switch method.type
            case 'integral'
                % TODO
                
            case 'trapz'
                % TODO
                
            case 'quad'
                switch prior.type
                    case 'uniform'
                        N = size(m,2);
                        wm_drift = estimator.wm_drift;
                        if N == 1
                            wms = wm;
                        elseif N == 2
                            wms = [wm_drift wm];
                        else
                            error('BLS_wm1wm2 estimator not supported for N > 2!')
                        end
                        
                        % Create x-vector
                        dx = method.dx;
                        x = xmin:dx:xmax;
                        
                        % Create Simpson's nodes
                        l = length(x);
                        h = (xmax - xmin)/l;
                        w = ones(1,l);
                        w(2:2:l-1) = 4;
                        w(3:2:l-1) = 2;
                        w = w*h/3;
                        
                        % Reshape measurements for processing
                        M = permute(m,[2 3 1]);
                        M = repmat(M,[1,1,1,l]);
                        x = reshape(x,[1 1 1 l]);
                        X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                        
                        % Generate estimate
                        w = reshape(w,[1 1 1 l]);
                        w = repmat(w,[1 1 size(m,1) 1]);
                        likelihood = likelihoods(M,X,wms);
                        e = sum(w.*X(1,:,:,:).*likelihood,4)./sum(w.*likelihood,4);
                        e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                        
                        
                    case 'Gaussian'
                        N = size(m,2);
                        wm_drift = estimator.wm_drift;
                        
                        % Parameterize prior
                        mu = prior.mu;
                        sig = prior.sig;
                        
                        if N == 1
                            % Create x-vector
                            dx = method.dx;
                            x = xmin:dx:xmax;
                            
                            % Create Simpson's nodes
                            l = length(x);
                            h = (xmax - xmin)/l;
                            w = ones(1,l);
                            w(2:2:l-1) = 4;
                            w(3:2:l-1) = 2;
                            w = w*h/3;
                            
                            % Reshape measurements for processing
                            M = permute(m,[2 3 1]);
                            M = repmat(M,[1,1,1,l]);
                            x = reshape(x,[1 1 1 l]);
                            X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                            
                            % Reshape measurements for processing
                            M = permute(m,[2 3 1]);
                            M = repmat(M,[1,1,1,l]);
                            x = reshape(x,[1 1 1 l]);
                            X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                            P = 1/sqrt(2*pi*sig^2) * exp( -(X(1,:,:,:)-mu).^2/(2*sig^2));
                            
                            % Generate estimate
                            w = reshape(w,[1 1 1 l]);
                            w = repmat(w,[1 1 size(m,1) 1]);
                            likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                            
                            e = sum(w.*X(1,:,:,:).*likelihood.*P,4)./sum(w.*likelihood.*P,4);
                            e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                            
                        elseif N == 2                         
                            
                            % Create x-vector
                            dx = method.dx;
                            x = xmin:dx:xmax;
                            
                            % Create Simpson'€™s nodes
                            l = length(x);
                            h = (xmax - xmin)/l;
                            w = ones(1,l);
                            w(2:2:l-1) = 4;
                            w(3:2:l-1) = 2;
                            w = w*h/3;
                            
                            % Reshape measurements for processing
                            M = permute(m,[2 3 1]);
                            M = repmat(M,[1,1,1,l]);
                            x = reshape(x,[1 1 1 l]);
                            X = repmat(x,[size(M,1) 1 size(M,3) 1]);
                            P = 1/sqrt(2*pi*sig^2) * exp( -(X(1,:,:,:)-mu).^2/(2*sig^2));
                            
                            % Generate estimate
                            w = reshape(w,[1 1 1 l]);
                            w = repmat(w,[1 1 size(m,1) 1]);
                            l1 = ( (1./sqrt(2*pi)/wm_drift/X(1,:,:,:)) .* ...
                                exp( -(X(1,:,:,:)-M(1,:,:,:)).^2 ./...
                                (2*wm_drift.^2.*X(1,:,:,:).^2) ) );
                            l2 = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .* ...
                                exp( -(X(1,:,:,:)-M(1,:,:,:)).^2 ./...
                                (2*wm.^2.*X(1,:,:,:).^2) ) );
                            likelihood = l1.*l2;
                            e = sum(w.*X(1,:,:,:).*likelihood.*P,4)./sum(w.*likelihood.*P,4);
                            e = permute(e,[3 2 1])/(1+wy.^2)^estimator.ObsAct;
                        else
                            error('N > 2 not supported for supOptMemBias model!')
                        end
                        
                    otherwise
                        error(['Prior ' prior.type ' not recognized!'])
                end
            case 'MonteCarlo'
                % TODO
                
            case 'MonteCarlo_batch'
                % TODO
        end
        
    case {'EKF'}
        switch method.type
            case 'integral'
                % TODO
                
            case 'trapz'
                % TODO
                
            case 'quad'
                estimatorBLS = estimator;
                 % Not fully generalized, but optimal for 1 measurement
%                 estimatorBLS.type = 'BLS';
%                 e = ScalarBayesEstimators(m(:,1),wm,xmin,xmax,'method',method,...
%                     'estimator',estimatorBLS);
%                 wmi = wm;
%                 if size(m,2) > 1
%                     for mi = 2:size(m,2)
%                         K = wmi^2/(wmi^2 + wm^2);           % Update gain
%                         wmi = wmi*wm/sqrt(wmi^2 + wm^2);    % Update weber fraction
%                         err = m(:,mi) - e;                % Calculate errors
%                         f_e = ScalarBayesEstimators(err + (xmin+xmax)/2,wm,xmin,xmax,'method',method,...
%                             'estimator',estimatorBLS) - (xmin+xmax)/2;          % Apply f_BLS nonlinearity to errors (centered on prior mean)
%                         e = e + K*f_e;
%                     end
%                 end
                
                % Fully generalized
                estimatorBLS.type = 'BLS';
                e0 = (xmin+xmax)/2;
                wmi = wm;
                for mi = 1:size(m,2)
                    if mi == 1
                        err = m(:,mi) - e0;
                        f_e = ScalarBayesEstimators(err+e0,wm,xmin,xmax,'method',method,...
                            'estimator',estimatorBLS) - e0;
                        e = e0 + f_e;
                    else
                        K = wmi^2/(wmi^2 + wm^2);           % Update gain
                        wmi = wmi*wm/sqrt(wmi^2 + wm^2);    % Update weber fraction
                        err = m(:,mi) - e;                % Calculate errors
                        f_e = ScalarBayesEstimators(err + e0,wm,xmin,xmax,'method',method,...
                            'estimator',estimatorBLS) - e0;          % Apply f_BLS nonlinearity to errors (centered on prior mean)
                        e = e + K*f_e;
                    end
                end
                
            case 'MonteCarlo'
                % TODO
                
            case 'MonteCarlo_batch'
                % TODO
        end
        
    case {'EKFk1k2'}
        switch method.type
            case 'integral'
                % TODO
                
            case 'trapz'
                % TODO
                
            case 'quad'
                estimatorBLS = estimator;
                 % Not fully generalized, but optimal for 1 measurement
%                 estimatorBLS.type = 'BLS';
%                 e = ScalarBayesEstimators(m(:,1),wm,xmin,xmax,'method',method,...
%                     'estimator',estimatorBLS);
%                 wmi = wm;
%                 if size(m,2) > 1
%                     for mi = 2:size(m,2)
%                         K = wmi^2/(wmi^2 + wm^2);           % Update gain
%                         wmi = wmi*wm/sqrt(wmi^2 + wm^2);    % Update weber fraction
%                         err = m(:,mi) - e;                % Calculate errors
%                         f_e = ScalarBayesEstimators(err + (xmin+xmax)/2,wm,xmin,xmax,'method',method,...
%                             'estimator',estimatorBLS) - (xmin+xmax)/2;          % Apply f_BLS nonlinearity to errors (centered on prior mean)
%                         e = e + K*f_e;
%                     end
%                 end
                
                % Fully generalized
                estimatorBLS.type = 'BLS';
                e0 = (xmin+xmax)/2;
                wmi = wm;
                for mi = 1:size(m,2)
                    if mi == 1
                        err = m(:,mi) - e0;
                        f_e = ScalarBayesEstimators(err+e0,wm,xmin,xmax,'method',method,...
                            'estimator',estimatorBLS) - e0;
                        e = e0 + estimator.k(mi)*f_e;
                    else
                        err = m(:,mi) - e;                % Calculate errors
                        f_e = ScalarBayesEstimators(err + e0,wm,xmin,xmax,'method',method,...
                            'estimator',estimatorBLS) - e0;          % Apply f_BLS nonlinearity to errors (centered on prior mean)
                        e = e + estimator.k(mi)*f_e;
                    end
                end
                
            case 'MonteCarlo'
                % TODO
                
            case 'MonteCarlo_batch'
                % TODO
        end
        
    case {'EKFf'}
        switch method.type
            case 'integral'
                % TODO
                
            case 'trapz'
                % TODO
                
            case 'quad'
                estimatorBLS.wy = estimator.wy;
                estimatorBLS.ObsAct = estimator.ObsAct;
                estimatorBLS.type = 'BLS';
                e = ScalarBayesEstimators(m(:,1),wm,xmin,xmax,'method',method,...
                    'estimator',estimatorBLS);
                wmi = wm;
                if size(m,2) > 1
                    for mi = 2:size(m,2)
                        K = wmi^2/(wmi^2 + wm^2);           % Update gain
                        wmi = wmi*wm/sqrt(wmi^2 + wm^2);    % Update weber fraction
                        err = m(:,mi) - e;                % Calculate errors
                        f_e = estimator.f(err);          % Apply f_BLS nonlinearity to errors (centered on prior mean)
                        e = e + K*f_e;
                    end
                end
                
            case 'MonteCarlo'
                % TODO
                
            case 'MonteCarlo_batch'
                % TODO
        end
end

%% Functions

function out = MonteCarloIntegrand_numerator(x,m,wm)

% Reshape measurements for processing
N = size(m,2);
l = length(x);
M = permute(m,[2 3 1]);
M = repmat(M,[1,1,1,l]);
x = reshape(x,[1 1 1 l]);
X = repmat(x,[size(M,1) 1 size(M,3) 1]);

% compute likelihood
likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
out = squeeze(permute(X(1,:,:,:).*likelihood,[4 3 2 1]));


function out = MonteCarloIntegrand_denominator(x,m,wm)

% Reshape measurements for processing
N = size(m,2);
l = length(x);
M = permute(m,[2 3 1]);
M = repmat(M,[1,1,1,l]);
x = reshape(x,[1 1 1 l]);
X = repmat(x,[size(M,1) 1 size(M,3) 1]);

% compute likelihood
likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(mmx('mult',permute(X-M,[2 1 3 4]),X-M))./(2*wm.^2.*X(1,:,:,:).^2) ) );
out = squeeze(permute(likelihood,[4 3 2 1]));

function out = posteriorMAP(x,m,wm,xmin,xmax)

N = size(m,2);
out = -(1./sqrt(2*pi)./wm./x).^N .* exp( -(sum((x-m).^2))./(2*wm.^2.*x.^2) );
out(x < xmin | x > xmax) = 0;

function out = posteriorMLE(x,m,wm,xmin,xmax)

N = size(m,2);
out = -(1./sqrt(2*pi)./wm./x).^N .* exp( -(sum((x-m).^2))./(2*wm.^2.*x.^2) );


function l = likelihoods(M,X,wms)

ls = nan(size(X));
for wmi = 1:length(wms)
    wm = wms(wmi);
    ls(wmi,:,:,:) = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .* ...
        exp( -(X(1,:,:,:)-M(wmi,:,:,:)).^2 ./...
        (2*wm.^2.*X(1,:,:,:).^2) ) );
end
l = prod(ls,1);
