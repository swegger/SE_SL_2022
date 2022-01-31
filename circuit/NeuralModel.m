function [n, M, rNN, e, tuning, w, sigG, G] = NeuralModel(varargin)
%% NeuralModel_v1
%
%
%
%%

%% Defaults
thetas_default = linspace(-90,90,20);
speeds_default = 2.^(linspace(-1,8,20)); %cumprod(2*ones(1,6));

theta_default.range = [-90,90,1800];
theta_default.Amp = 10;
theta_default.sig = 45;

speed_default.range = [-1,8,1000];
speed_default.Amp = 10;
speed_default.sig = 1.64;
speed_default.d = 0.1;

Cov_default.sigf = 0.55;
Cov_default.thetaLengthConstant = 0.4;
Cov_default.speedLengthConstant = 0.3;
Cov_default.separationLengthConstant = 0.3;
Cov_default.alpha = 0;
Cov_default.diffAlpha = 0;

sizeProps_default.minEccentricity = 1;
sizeProps_default.maxEccentricity = 30;

saveOpts_default.On = false;
saveOpts_default.Figs = false;

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'sizes',[2,6,20])
addParameter(Parser,'thetas',thetas_default)
addParameter(Parser,'speeds',speeds_default)
addParameter(Parser,'theta',theta_default)
addParameter(Parser,'speed',speed_default)
addParameter(Parser,'Cov',Cov_default)
addParameter(Parser,'n0',1)
addParameter(Parser,'epsilon',0.05)
addParameter(Parser,'gainNoise',0)
addParameter(Parser,'motorNoise',0)
addParameter(Parser,'sizeProps',sizeProps_default)
addParameter(Parser,'N',NaN)
addParameter(Parser,'decoderAlgorithm','g*log2shat')
addParameter(Parser,'mymakeaxisflg',true)
addParameter(Parser,'plotMT',true)
addParameter(Parser,'plotDecoding',true)
addParameter(Parser,'plotResults',true)
addParameter(Parser,'saveOpts',saveOpts_default)

parse(Parser,varargin{:})

sizes = Parser.Results.sizes;
thetas = Parser.Results.thetas;
speeds = Parser.Results.speeds;
theta = Parser.Results.theta;
speed = Parser.Results.speed;
Cov = Parser.Results.Cov;
n0 = Parser.Results.n0;
epsilon = Parser.Results.epsilon;
gainNoise = Parser.Results.gainNoise;
motorNoise = Parser.Results.motorNoise;
sizeProps = Parser.Results.sizeProps;
N = Parser.Results.N;
decoderAlgorithm = Parser.Results.decoderAlgorithm;
mymakeaxisflg = Parser.Results.mymakeaxisflg;
plotMT = Parser.Results.plotMT;
plotDecoding = Parser.Results.plotDecoding;
plotResults = Parser.Results.plotResults;
saveOpts = Parser.Results.saveOpts;

%% Generate population size
% fun=@(x1) (1.14*x1.^-0.76); % Albright & Desimone 87 Exp Brain Res (mm/deg)
fun=@(x1) (6*x1.^-0.9); % Erickson cell density function (mm/deg)
%N = ceil(20*(integral(fun,sizeProps.minEccentricity,sizeProps.maxEccentricity)+9));  % By Pack Lab. Should probably multiply by 2!!!!!
if isnan(N)
    N = ceil( 20* ( 6*( 10*sizeProps.maxEccentricity.^(1/10) - 10*sizeProps.minEccentricity.^(1/10) ) + 9 ) ); % Definite integral of cell density function; should probably multiply by 2!!!!!
end

%% Initial simulation w/ maximum number of neurons
% Generate population tuning parameters
tuning = tuningFunctions(N,theta,speed,Cov,n0,sizeProps);

% Simulate MT and then decode
[nTemp, MTemp, rNNTemp, ~, ~] = DirSpeedSizeLocMT(thetas,speeds,max(sizes),'trialN',400,'tuning',tuning,'plotflg',false);

[Ds,Ss] = meshgrid(thetas,speeds);
s = cat(3,Ds',Ss');
switch decoderAlgorithm
    case {'bioRxiv2020'}
        gainTemp = gainFunction(nTemp(1,:,:,:),tuning,0);
        normalizer = [mean(gainTemp(:))/mean(log2(speeds)), 0]; % For -90 90 theta range
        
    case {'simpleSumVA^2'}
        [e_temp, gain_temp] = DecodeMT(nTemp,tuning,s,...
            'gainNoise',0,'epsilon',epsilon,'b',[1 0],'decoderAlgorithm',decoderAlgorithm,...
            'mymakeaxisflg',false,'plotflg',false,...
            'motorNoise',motorNoise);
        normalizer = [sqrt(sum(sum(e_temp(1,:,:,2),3))/numel(e_temp(1,:,:,2)))/10, 0];
        
    case {'gainConstant'}
        [e_temp, gain_temp] = DecodeMT(nTemp,tuning,s,...
            'gainNoise',0,'epsilon',epsilon,'b',[1 0],'decoderAlgorithm',decoderAlgorithm,...
            'mymakeaxisflg',false,'plotflg',false,...
            'motorNoise',motorNoise);
        normalizer = [sum(sum(log2(e_temp(1,:,:,2)),3))/numel(e_temp(1,:,:,2))/log2(10), 0];
        
    otherwise
        [e_temp, gain_temp] = DecodeMT(nTemp,tuning,s,...
            'gainNoise',0,'epsilon',epsilon,'b',[1 0],'decoderAlgorithm',decoderAlgorithm,...
            'mymakeaxisflg',false,'plotflg',false,...
            'motorNoise',motorNoise);
        normalizer = [sum(sum(e_temp(1,:,:,2),3))/numel(e_temp(1,:,:,2))/10, 0];
end
    
%% Run MT simulations
for szi = 1:length(sizes)

    [n{szi}, M{szi}, rNN{szi}, ~, ~] = DirSpeedSizeLocMT(thetas,speeds,sizes(szi),'trialN',400,'tuning',tuning,'plotflg',plotMT,'mymakeaxisflg',mymakeaxisflg);

end

%% Run decoding simulations

% Decoder properties
% [Ds,Ss] = meshgrid(thetas,speeds);
% s = cat(3,Ds',Ss');

for szi = 1:length(sizes)
    [e{szi}, gain{szi}, figureHandles, Rs{szi}, PVals{szi}] = DecodeMT(n{szi},tuning,s,...
        'gainNoise',gainNoise,'epsilon',epsilon,'b',normalizer,'decoderAlgorithm',decoderAlgorithm,...
        'mymakeaxisflg',false,'plotflg',plotDecoding,...
        'motorNoise',motorNoise);

    if saveOpts.Figs
        savefig(figureHandles.neuon_est_corr,[saveOpts.location '_MT_pursuit_corr_' num2str(sizes(szi)) ])
    end

    eBar = mean(e{szi},3);
    eVar = var(e{szi},1,3);

    VeM(:,:,szi) = squeeze(eBar(:,:,:,2));
    VeVAR(:,:,szi) = squeeze(eVar(:,:,:,2));

%     figure(h)
%     subplot(1,length(N),szi)
%     plot(speeds,squeeze(e{szi}(1,:,:,2)),'ko')
%     hold on
%     plotUnity;
end

%% Fit gainSDN model across sizes
OPTIONS = optimset('Display','off');
[w,sigG] = fit_gainSDN(speeds,VeM,VeVAR,0.1,gainNoise,OPTIONS);


%% Estimate gain{szi}~ain{szi}
for szi = 1:length(sizes)
    ys = e{szi}(1,:,:,2);
    xs = repmat(speeds,[size(ys,1),1,size(ys,3)]);
    betas(:,szi) = regress(ys(:),[xs(:) ones(size(xs(:)))]);
    G(:,szi) = betas(1,szi);
end

%% Plotting
if plotResults
    szcolors = [0 1 0;...
        1 0 0;...
        0   0   0];

    h = figure('Name','Target v Eye speed','Position',[26 366 621 387]);
    for szi = 1:length(sizes)
        szind = length(sizes)-szi+1;
        plot(speeds,squeeze(mean(e{szind}(1,:,:,2),3)),'o-',...
            'Color',szcolors(szind,:),'MarkerFaceColor',szcolors(szind,:))
        hold on
    end
    axis square
    xlabel('Target speed (deg/s)')
    ylabel('Eye speed (deg/s)')
    plotUnity;
    if mymakeaxisflg
        mymakeaxis(gca,'xticks',[0,10,20],'yticks',[0 10 20]);
    end

    if saveOpts.Figs
        savefig(h,[saveOpts.location '_targVeye'])
    end

    h = figure;
    Ncolors = colormap('lines');


    Ncolors = [szcolors; Ncolors];

    colors = projectColorMaps('speeds','samples',1:length(speeds),...
        'sampleDepth',length(speeds));

    % figure
    for di = 1:length(thetas)
        subplot(1,length(thetas),di)
        for szi = 1:length(sizes)
            szind = length(sizes)-szi+1;
            plot(permute(VeM(di,:,szind),[2,3,1]),permute(VeVAR(di,:,szind),[2,3,1]),...
                'o','Color',Ncolors(szind,:),'MarkerFaceColor',Ncolors(szind,:),'MarkerSize',10)
            hold on
            x = linspace(0,max(speeds),100);
            plot(betas(1,szind)*x+betas(2,szind),gainSDN(betas(1,szind)*x+betas(2,szind),x,w,sigG),'-','Color',Ncolors(szind,:))
            xlabel('Mean eye speed (deg/s)')
            ylabel('Eye speed variance (deg/s)^2')
        end

        %     ylim([0 3.5])
        axis square
        if mymakeaxisflg
            mymakeaxis(gca);
        end
    end

    if saveOpts.Figs
        savefig(h,[saveOpts.location 'muVvar'])
    end

    %% Tuning
    h = figure('Name','MT RF centers','Position',[713 316 1530 420]);
    subplot(1,2,1)
    plot(tuning.size.x,tuning.size.y,'k.')
    hold on
    for szi = 1:length(sizes)
        [x,y] = ellipse(sizes(szi)/2,sizes(szi)/2,0,0,pi/360);
        plot(x,y,'r')
    end
    axis([-1.2*sizeProps.maxEccentricity 1.2*sizeProps.maxEccentricity -1.2*sizeProps.maxEccentricity 1.2*sizeProps.maxEccentricity])
    axis square
    ind1 = find(sqrt(tuning.size.x.^2+tuning.size.y.^2) > 9 & sqrt(tuning.size.x.^2+tuning.size.y.^2) < 15,1);
    [x,y] = ellipse(tuning.size.radius(ind1),tuning.size.radius(ind1),tuning.size.x(ind1),tuning.size.y(ind1),pi/360);
    plot(x,y,'k')
    ind2 = find(sqrt(tuning.size.x.^2+tuning.size.y.^2) < 6 & sqrt(tuning.size.x.^2+tuning.size.y.^2) > 3 & tuning.speed.pref > 10 ...
        & tuning.theta.pref < 30 & tuning.theta.pref > -30 & tuning.theta.Amp > 50,1);
    [x,y] = ellipse(tuning.size.radius(ind2),tuning.size.radius(ind2),tuning.size.x(ind2),tuning.size.y(ind2),pi/360);
    plot(x,y,'-','Color',[0,174,239]/255)
    plotVertical(0);
    plotHorizontal(0);
    xlabel('Horizontal position (deg)')
    ylabel('Vertical position (deg)')
    if mymakeaxisflg
        mymakeaxis(gca);
    end

    subplot(1,2,2)
    for szi = 1:length(sizes)
        plot(sqrt(tuning.size.x.^2+tuning.size.y.^2),squeeze(mean(M{szi},2)),'o')
        hold on
        plot(sqrt(tuning.size.x(ind1).^2+tuning.size.y(ind1).^2),squeeze(mean(M{szi}(:,:,ind1),2)),...
            'ko','MarkerFaceColor',[0 0 0])
        plot(sqrt(tuning.size.x(ind2).^2+tuning.size.y(ind2).^2),squeeze(mean(M{szi}(:,:,ind2),2)),...
            'ko','MarkerFaceColor',[0 0 0])
        plotVertical(sizes(szi)/2);
    end

    xlabel('Eccentricity (deg)')
    ylabel('Mean response')
    if mymakeaxisflg
        mymakeaxis(gca);
    end

    ExInd = ind2; %find(tuning.speed.pref > 15,1);
    figure('Name','Popultion','Position',[440 31 559 767])
%     f = @(x,p)(tuning.theta.Amp .* exp( -(x-p).^2 ./ tuning.theta.sig.^2 /2 ));
    x = linspace(tuning.theta.range(1),tuning.theta.range(2),200);
    randN = 50;
    if randN > N
        randN = N;
    end
    randInds = randsample(N,randN);
    subplot(4,1,1)
    for ni = 1:randN
        if length(tuning.theta.Amp) > 1
            f = tuning.theta.Amp(randInds(ni)) .* ...
                exp( -(x-tuning.theta.pref(randInds(ni))).^2 ./ ...
                tuning.theta.sig(randInds(ni)).^2 /2 );
        else
            f = tuning.theta.Amp .* ...
                exp( -(x-tuning.theta.pref(randInds(ni))).^2 ./ ...
                tuning.theta.sig.^2 /2 );
        end
        plot(x,f,'Color',[0.6 0.6 0.6])
%         plot(x,f(x,tuning.theta.pref(randInds(ni))),'Color',[0.6 0.6 0.6])
        hold on
    end
    if length(tuning.theta.Amp) > 1
        f = tuning.theta.Amp(ExInd) .* ...
            exp( -(x-tuning.theta.pref(ExInd)).^2 ./ ...
            tuning.theta.sig(ExInd).^2 /2 );
    else
        f = tuning.theta.Amp .* ...
            exp( -(x-tuning.theta.pref(ExInd)).^2 ./ ...
            tuning.theta.sig.^2 /2 );
    end
    if ~isempty(f)
    plot(x,f,'k','LineWidth',2)
    end
%     plot(x,f(x,tuning.theta.pref(ExInd)),'k','LineWidth',2)
    xlabel('Direction (deg)')
    ylabel('Spikes/s')
    axis tight

    subplot(4,1,2)
    x = logspace(log10(2^(tuning.speed.range(1))),log10(2^(tuning.speed.range(2))),200);
%     x = linspace(2^(tuning.speed.range(1)),2^(tuning.speed.range(2)),200);
%     f = @(x,p)(tuning.speed.Amp .* exp( -(log2(x./p)).^2 ./ tuning.speed.sig ));
    for ni = 1:randN
        if length(tuning.theta.Amp) > 1
            f = tuning.theta.Amp(randInds(ni)) .* ...
                exp( -(log2(x./tuning.speed.pref(randInds(ni)))).^2 ./ ...
                tuning.speed.sig(randInds(ni)) );
        else
            f = tuning.theta.Amp .* ...
                exp( -(log2(x./tuning.speed.pref(randInds(ni)))).^2 ./ ...
                tuning.speed.sig );
        end
        plot(x,f,'Color',[0.6 0.6 0.6])
%         plot(x,f(x,tuning.speed.pref(randInds(ni))),'Color',[0.6 0.6 0.6])
        hold on
    end
    if length(tuning.theta.Amp) > 1
        f = tuning.theta.Amp(ExInd) .* ...
            exp( -(log2(x./tuning.speed.pref(ExInd))).^2 ./ ...
            tuning.speed.sig(ExInd) );
    else
        f = tuning.theta.Amp .* ...
            exp( -(log2(x./tuning.speed.pref(ExInd))).^2 ./ ...
            tuning.speed.sig );
    end
    if ~isempty(f)    
    plot(x,f,'k','LineWidth',2);
    end
%     plot(x,f(x,tuning.speed.pref(ExInd)),'k','LineWidth',2);
    set(gca,'XScale','log')
    xlabel('log(speed)')
    ylabel('Spikes/s')
    axis tight

    subplot(4,1,[3,4])
    scatter(tuning.speed.pref,tuning.theta.pref,40,squeeze(n{3}(1,4,20,:)),'filled')
    hold on
    scatter(tuning.speed.pref(ExInd),tuning.theta.pref(ExInd),200,squeeze(n{3}(1,4,20,(ExInd))),'filled')
    hax = gca;
    hax.XScale = 'log';
    hax.TickDir = 'out';
    cax = myrgb(64,[0,0,0]/255,[255,0,0]/255);
    colormap(cax)
    colorbar
    axis tight
    plotVertical(speeds(4));
    axis square
    xlabel('log(Pref speed)')
    ylabel('Pref direction')

    %% Gain
    h = figure('Name','Gain vs speed','Position',[26 366 621 387]);
    for szi = 1:length(sizes)
        plot(speeds,1e4*squeeze(mean(gain{szi}(1,:,:),3))/normalizer(1),'o-',...
            'Color',szcolors(szi,:),'MarkerFaceColor',szcolors(szi,:))
        hold on
%         plot(speeds,squeeze(gain{szi}(1,:,randsample(size(gain{szi},3),100))),'o',...
%             'Color',szcolors(szi,:),'MarkerFaceColor',szcolors(szi,:))
%         hold on
    end
    axis square
    xlabel('Target speed (deg/s)')
    ylabel('Gain')
    if mymakeaxisflg
        mymakeaxis(gca,'xticks',[0,10,20]);
    end

    %% Example correlation
    zE = (squeeze(e{3}(1,4,:,2))-mean(squeeze(e{3}(1,4,:,2))))/std(squeeze(e{3}(1,4,:,2)));
    zN = (squeeze(n{3}(1,4,:,ind2))-mean(squeeze(n{3}(1,4,:,ind2))))/std(squeeze(n{3}(1,4,:,ind2)));

    figure('Name','Neuron estimate correlation')
    scatter(zN,zE,80,[0 0 0],'filled')
    hold on
    axis equal
    plotUnity;
    axis square
    plotHorizontal(0);
    plotVertical(0);
    xlabel('z-score (neuron)')
    ylabel('z-score (pursuit)')

    mymakeaxis(gca)
end

%% Save
if saveOpts.On
    save(saveOpts.location,'-v7.3')
end

%% Functions

%% sampleEccentricities
function eccentricities = sampleEccentricities(N,minEccentricity,maxEccentricity)
    foveaN = 9*20; % 9*20 neurons from 0 to 1 deg
    u = rand(N-foveaN,1);
    eccentricities = ((maxEccentricity^0.1-1)*u+1).^10;
    foveaMin = 0.25;
    eccentricities = [(1-foveaMin)*rand(foveaN,1)+foveaMin; ...
        eccentricities];

    eccentricities = randsample(eccentricities,N);

%% tuningFunctions
function tuning = tuningFunctions(N,theta,speed,Cov,n0,sizeProps)

    thetas = linspace(theta.range(1),theta.range(2),theta.range(3));
    speeds = 2.^(linspace(speed.range(1),speed.range(2),speed.range(3)));

    A(:,1) = randsample(thetas,N,true);
    A(:,2) = randsample(speeds,N,true);
    B = sortrows(A,[2,1]);

    tuning.theta = theta;
    tuning.theta.pref = B(:,1);
    if isfield(theta,'widthRange')
        dirwidths = linspace(theta.widthRange(1),theta.widthRange(2),theta.widthRange(3));
        if length(dirwidths) == 1
            dirwidths = dirwidths*ones(N,1);
        end
        tuning.theta.sig(:,1) = randsample(dirwidths,N,true);
    end
    if isfield(theta,'amplitudeRange')
        amps = linspace(theta.amplitudeRange(1),theta.amplitudeRange(2),theta.amplitudeRange(3));
        if length(amps) == 1
            amps = amps*ones(N,1);
        end
        tuning.theta.Amp(:,1) = randsample(amps,N,true);
    end

    tuning.speed = speed;
    tuning.speed.pref = B(:,2);
    if isfield(speed,'widthRange')
        speedwidths = linspace(speed.widthRange(1),speed.widthRange(2),speed.widthRange(3));
        if length(speedwidths) == 1
            speedwidths = speedwidths*ones(N,1);
        end
        tuning.speed.sig(:,1) = randsample(speedwidths,N,true);
    end
%     if isfield(speed,'amplitudeRange')
%         amps = linspace(speed.amplitudeRange(1),speed.amplitudeRange(2),speed.amplitudeRange(3));
%         if length(amps) == 1
%             amps = amps*ones(N,1);
%         end
%         tuning.speed.Amp(:,1) = randsample(amps,N,true);
%     end

    eccentricities = sampleEccentricities(N,sizeProps.minEccentricity,sizeProps.maxEccentricity);
    rhos = 2*pi*rand(N,1);
    tuning.size.x = eccentricities.*cos(rhos);
    tuning.size.y = eccentricities.*sin(rhos);
    tuning.size.radius = (0.69*eccentricities+1)/sqrt(pi); % Raduis approixmation from Ferrera and Lisberger (1997)
    if isfield(sizeProps,'threshold')
        tuning.size.threshold = sizeProps.threshold.*ones(N,1);
    else
        tuning.size.threshold = zeros(N,1);
    end
    if isfield(sizeProps,'exponential')
        tuning.size.exponential = sizeProps.exponential.*ones(N,1);
    else
        tuning.size.exponential = ones(N,1);
    end
    if isfield(sizeProps,'surround_weight')
        tuning.size.surround_weight = sizeProps.surround_weight.*ones(N,1);
    else
        tuning.size.surround_weight = ones(N,1);%zeros(N,1);
    end

    tuning.n0 = n0;

    tuning.Cov = Cov;

%% fit_gainSDN
function [w, sigG] = fit_gainSDN(speeds,VeM,VeVAR,w0,sigG0,OPTIONS)
    %%
    minimizer = @(p)( sum( sum( (VeVAR - gainSDN(VeM,repmat(speeds(:)',[size(VeM,1),1,size(VeM,3)]),p(1),p(2))).^2 ) ) );
    p = fmincon(minimizer,[w0 sigG0],[],[],[],[],[0,0],[Inf,Inf],[],OPTIONS);
    w = p(1);
    sigG = p(2);

%% gainSDN
function v = gainSDN(ve,vs,w,sigG)
    %%
    v = w.^2.*ve.^2 + (sigG.^2 + sigG.^2.*w.^2).*vs.^2;

%% gainFunction
function out = gainFunction(n,tuning,gainNoise)
    %%
    for thetai = 1:size(n,1)
        for speedi = 1:size(n,2)
            out(thetai,speedi,:) = log2(tuning.speed.pref)' ...
                * permute(n(thetai,speedi,:,:),[4,3,1,2]) + ...
                gainNoise*randn(size(n,3),1)';

            % SDN
            %         out(thetai,speedi,:) = log2(tuning.speed.pref)' ...
            %             * permute(n(thetai,speedi,:,:),[4,3,1,2]);
            %         out(thetai,speedi,:) = out(thetai,speedi,:) + ...
            %             out(thetai,speedi,:)*gainNoise.*randn(1,1,size(n,3));
        end
    end
