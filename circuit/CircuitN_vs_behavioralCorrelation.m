function CircuitN_vs_behavioralCorrelation(varargin)
%% CircuitN_vs_behavioralCorrelation
%
%   Performs analysis of MT-behavior correlation and shift in signal 
%   dependent noise for biomimetic model of different populatin sizes.
%
%%

%% Defaults
% N = logspace(1,4,4);
% N = [20 40 80 160 320 640 1280 2560 5120];
% 
% baseName = '_vanilla_gainNoiseOff_20200504.mat';
% baseNameNoise = '_vanilla_gainNoiseOn_20200504.mat';
noiseInteractionSim_default.On = false;
OPTIONS_default = optimset('Display','off');

%% Parse inputs
Parser = inputParser;

addParameter(Parser,'N',[20 40 80 160 320 640 1280 2560 5120 10240])
addParameter(Parser,'baseName','_vanilla_gainNoiseOff_20200504.mat')
addParameter(Parser,'baseNameNoise','_vanilla_gainNoiseOn_20200504.mat')
addParameter(Parser,'noiseInteractionSim',noiseInteractionSim_default)
addParameter(Parser,'quants',linspace(0.1,1,10))
addParameter(Parser,'OPTIONS',OPTIONS_default)

parse(Parser,varargin{:})

N = Parser.Results.N;
baseName = Parser.Results.baseName;
baseNameNoise = Parser.Results.baseNameNoise;
noiseInteractionSim = Parser.Results.noiseInteractionSim;
quants = Parser.Results.quants;
OPTIONS = Parser.Results.OPTIONS;

%% Iterate for each population size

for Ni = 1:length(N)
    load(['N' num2str(N(Ni)) baseName],'Rs','w','sigG','VeM','VeVAR','tuning','s')
    for szi = 1:3
        dX = repmat(s(:,:,1),[1,1,size(Rs{szi},3),size(Rs{szi},4)]) - ...
            repmat(permute(tuning.theta.pref,[3,2,1]),[size(s,1),size(s,2),1,2]);
        quantsR(Ni,szi,:) = quantile(Rs{szi}(dX < 45 & dX > -45 | dX > 135 | dX < -135),quants);
        quantsR2(Ni,szi,:) = quantile(Rs{szi}(dX < 45 & dX > -45 | dX > 135 | dX < -135).^2,quants);
        maxR(Ni,szi) = max(Rs{szi}(dX < 45 & dX > -45 | dX > 135 | dX < -135));
        rvalsModel{Ni,szi} = Rs{szi}(dX < 45 & dX > -45 | dX > 135 | dX < -135);
    end
%     subplot(length(N),1,Ni)
%     histogram(tuning.theta.Amp.*tuning.speed.Amp,linspace(0,400,200))
%     histogram(Rs{2}(:),linspace(-1,1,200))
%     hold on
%     plotVertical(maxR(Ni,2));
%     plotVertical(0.3378);
    wfit(Ni) = w;
    sigFit(Ni) = sigG;
    
    for szi = 1:3;%[1,3]
        w_standard(Ni,szi) = fit_standardSDN(...
            VeM(:,:,szi),VeVAR(:,:,szi),0.1,OPTIONS); 
    end
    
    varSummary(Ni,:) = mean(VeVAR,2);
    
    load(['N' num2str(N(Ni)) baseNameNoise],'Rs','w','sigG','VeM','VeVAR','gainNoise','tuning')
    for szi = 1:3
        dX = repmat(s(:,:,1),[1,1,size(Rs{szi},3),size(Rs{szi},4)]) - ...
            repmat(permute(tuning.theta.pref,[3,2,1]),[size(s,1),size(s,2),1,2]);
        quantsRnoise(Ni,szi,:) = quantile(Rs{szi}(dX < 45 & dX > -45 | dX > 135 | dX < -135),quants);
        quantsR2noise(Ni,szi,:) = quantile(Rs{szi}(dX < 45 & dX > -45 | dX > 135 | dX < -135).^2,quants);
        maxRnoise(Ni,szi) = max(Rs{szi}(dX < 45 & dX > -45 | dX > 135 | dX < -135));
        
        rvalsModelnoise{Ni,szi} = Rs{szi}(dX < 45 & dX > -45 | dX > 135 | dX < -135);
    end
    wfit_noise(Ni) = w;
    sigFit_noise(Ni) = sigG;
    
    for szi = 1:3;%[1,3]
        w_standard_noise(Ni,szi) = fit_standardSDN(...
            VeM(:,:,szi),VeVAR(:,:,szi),0.1,OPTIONS); 
    end
    
    varSummary_noise(Ni,:) = mean(VeVAR,2);
    
    if noiseInteractionSim.On
        for ni = 1:length(noiseInteractionSim.gainNoise)
            disp(['Sim ' num2str(ni) ' of ' num2str(length(N)) ' for N = ' num2str(N(ni))])
            for szi = 1:3
                load(['N' num2str(N(Ni)) baseNameNoise],'n','tuning','s','epsilon','normalizer',...
                    'decoderAlgorithm','motorNoise')
                [e, ~, ~, Rs, ~] = DecodeMT(n{szi},tuning,s,...
                    'gainNoise',noiseInteractionSim.gainNoise(ni),...
                    'epsilon',epsilon,'b',normalizer,'decoderAlgorithm',decoderAlgorithm,...
                    'mymakeaxisflg',false,'plotflg',false,...
                    'motorNoise',motorNoise);
                
                eBar = mean(e,3);
                eVar = var(e,1,3);
                
                VeM(:,:,szi) = squeeze(eBar(:,:,:,2));
                VeVAR(:,:,szi) = squeeze(eVar(:,:,:,2));
                wtemp(szi) = fit_standardSDN(...
                    VeM(:,:,szi),VeVAR(:,:,szi),0.1,OPTIONS);
                
                noiseInteractionSim.maxR(Ni,ni,szi) = max(Rs(:));
            end
            noiseInteractionSim.deltaW(Ni,ni,1) = wtemp(1) - wtemp(3);
            noiseInteractionSim.deltaW(Ni,ni,2) = wtemp(2) - wtemp(3);
            noiseInteractionSim.varSummary(Ni,ni,:) = mean(VeVAR,2);
        end
    end
    
end

%% Load in hohl et al data
hohl = load('~/Documents/Publications/GainNoise/Submissions/NatComms/Submission2/mat/hohl_et_al_speedCorr.mat');
speedCorr_data = [hohl.hohletalsig; hohl.hohletalnotsig];
hohlRquants = quantile(speedCorr_data(:,2),quants);
hohlR2quants = quantile(speedCorr_data(:,2).^2,quants);

%%
figure('Name','Maximum correlation')
for szi = 1:3
    subplot(1,3,szi)
    plot(log2(N/10),maxR(:,szi),'-o')
    hold on
    plot(log2(N/10),maxRnoise(:,szi),'-o')
    xlabel('Population size')
    ylabel('Maximum MT-behavior correlation')
    axis([1 log2(max(N)/10) 0.2 0.8])
    plotHorizontal( 0.3378 );
end
legend('\sigma = 0',['\sigma = ' num2str(gainNoise)])

figure('Name','w and \sigma')
subplot(1,2,1)
plot(N,wfit,'-o')
hold on
plot(N,wfit_noise,'-o')
xlabel('Population size')
ylabel('w')
legend('\sigma = 0',['\sigma = ' num2str(gainNoise)])

subplot(1,2,2)
plot(N,sigFit,'-o')
hold on
plot(N,sigFit_noise,'-o')
xlabel('Population size')
ylabel('\sigma')

%%
deltaW = w_standard(:,1) - w_standard(:,end);
deltaW_noise = w_standard_noise(:,1) - w_standard_noise(:,end);
figure('Name','w and \sigma')
subplot(1,2,1)
plot(log2(N/10),deltaW,'-o')
hold on
plot(log2(N/10),deltaW_noise,'-o')
plotHorizontal( 0.0698 );
plotHorizontal( 0.1034 );
ylabel('w')
legend('\sigma = 0',['\sigma = ' num2str(gainNoise)])

deltaW = w_standard(:,2) - w_standard(:,end);
deltaW_noise = w_standard_noise(:,2) - w_standard_noise(:,end);
subplot(1,2,2)
plot(log2(N/10),deltaW,'-o')
hold on
plot(log2(N/10),deltaW_noise,'-o')
%plotHorizontal( 0.0698 );
%plotHorizontal( 0.1034 );
plotHorizontal( 0.0317 );       
plotHorizontal( 0.0522 );
xlabel('Population size')
ylabel('w')
legend('\sigma = 0',['\sigma = ' num2str(gainNoise)])

%%
colors = [0 1 0; 1 0 0; 0 0 0];
figure('Name','Standard SDN with population N')
for szi = [1,3]
    plot(N,w_standard(:,szi),'o-','Color',colors(szi,:))
    hold on
end
xlabel('Population size')
ylabel('w')

%% Variance-limiting correlations

colors = [0 1 0; 1 0 0; 0 0 0];
figure('Name','Pursuit variance with population N')
subplot(1,2,1)
for szi = 1:3
    plot(log2(N/10),log2(varSummary(:,szi)),'-o','Color',colors(szi,:))
    hold on
end
axis tight
ax(:,1) = axis;

subplot(1,2,2)
for szi = 1:3
    plot(log2(N/10),log2(varSummary_noise(:,szi)),'-o','Color',colors(szi,:))
    hold on
end
axis tight
ax(:,2) = axis;

subplot(1,2,1)
axis([min(ax(1,:)),max(ax(2,:)),min(ax(3,:)),max(ax(4,:))])
xlabel('Population size')
ylabel('Mean variance')

subplot(1,2,2)
axis([min(ax(1,:)),max(ax(2,:)),min(ax(3,:)),max(ax(4,:))])
xlabel('Population size')
ylabel('Mean variance')

figure('Name','Pursuit variance with population N')
for szi = 1:3
    subplot(1,3,szi)
    plot(log2(N/10),log2(varSummary(:,szi)),'-o')
    hold on
    plot(log2(N/10),log2(varSummary_noise(:,szi)),'-o')
    ax(:,szi) = axis;
end
for szi = 1:3
    subplot(1,3,szi)
    axis([min(ax(1,:)),max(ax(2,:)),min(ax(3,:)),max(ax(4,:))])
    xlabel('Population size')
    ylabel('Mean variance')
end
legend('\sigma = 0',['\sigma = ' num2str(gainNoise)])

%% Noise/N interaction
if noiseInteractionSim.On
    colors = projectColorMaps('popN','samples',1:length(N),...
        'sampleDepth',length(N));
    figure
    for szi = 1:3
        subplot(1,3,szi)
        for Ni = 1:length(N)
            plot(noiseInteractionSim.maxR(Ni,:,szi),noiseInteractionSim.varSummary(Ni,:,szi),'o-','Color',colors(Ni,:))
            hold on
        end
        
        for Ni = 1:length(N)
            plot(noiseInteractionSim.maxR(Ni,1,szi),noiseInteractionSim.varSummary(Ni,1,szi),...
                'o','Color',colors(Ni,:),'MarkerFaceColor',colors(Ni,:),'MarkerSize',10)
        end
        axis tight
        plotVertical( 0.3378 );
        plotHorizontal( 0.7773 );
        plotHorizontal( 1.6049 );
        xlabel('Maximum correlation')
        ylabel('Variance')
    end
    
    figure
    for szi = 1:3
        subplot(2,3,szi)
        for Ni = 1:length(N)
            plot(noiseInteractionSim.maxR(Ni,:,szi),noiseInteractionSim.deltaW(Ni,:,1),'o-','Color',colors(Ni,:))
            hold on
        end
        for Ni = 1:length(N)
            plot(noiseInteractionSim.maxR(Ni,1,szi),noiseInteractionSim.deltaW(Ni,1,1),...
                'o','Color',colors(Ni,:),'MarkerFaceColor',colors(Ni,:),'MarkerSize',10)
        end
        axis tight
        plotHorizontal( 0.0698 );
        plotHorizontal( 0.1034 );
        plotVertical( 0.3378 );
        xlabel('Maximum correlation')
        ylabel('\Delta W')
        
        subplot(2,3,szi+3)
        for Ni = 1:length(N)
            plot(noiseInteractionSim.maxR(Ni,:,szi),noiseInteractionSim.deltaW(Ni,:,2),'o-','Color',colors(Ni,:))
            hold on
        end
        
        for Ni = 1:length(N)
            plot(noiseInteractionSim.maxR(Ni,1,szi),noiseInteractionSim.deltaW(Ni,1,2),...
                'o','Color',colors(Ni,:),'MarkerFaceColor',colors(Ni,:),'MarkerSize',10)
        end
        axis tight
        plotHorizontal( 0.0317 );
        plotHorizontal( 0.0522 );
        plotVertical( 0.3378 );
        xlabel('Maximum correlation')
        ylabel('\Delta W')
    end
    
    figure
    for szi = 1:3
        subplot(1,3,szi)
        for Ni = 1:length(N)
            plot3(noiseInteractionSim.maxR(Ni,:,szi),...
                noiseInteractionSim.varSummary(Ni,:,szi),...
                noiseInteractionSim.deltaW(Ni,:,2),...
                'o-','Color',colors(Ni,:))
            hold on
            
            patch([0.3378 0.3378 0.3378 0.3378 0.3378],...
                [0.7773 1.6049 1.6049 0.7773 0.7773],...
                [0.0317 0.0317 0.0522 0.0522 0.0317],...
                [0.1 0.1 0.1])
        end 
            
        for Ni = 1:length(N)
            plot3(noiseInteractionSim.maxR(Ni,1,szi),...
                noiseInteractionSim.varSummary(Ni,1,szi),...
                noiseInteractionSim.deltaW(Ni,1,2),...
                'o','Color',colors(Ni,:),'MarkerFaceColor',colors(Ni,:))
        end
        xlabel('Maximum correlation')
        ylabel('Variance')
        zlabel('\Delta w')
        grid on
        view([120 20])
    end
    
end

%% Quantile vs quantile
figure('Name','Quantiles','Position',[50 169 1370 627])
% for Ni = 1:length(N)
%     subplot(2,5,Ni)
Ni = 7;
    for szi = 1:3
        plot(hohlRquants,permute(quantsR(Ni,szi,:),[3,2,1]),'o--','Color',colors(szi,:))
        hold on
        plot(hohlRquants,permute(quantsRnoise(Ni,szi,:),[3,2,1]),'o-','Color',colors(szi,:))
        axis([-0.5 0.8 -0.5 0.8])
        plotUnity;
        axis square
        title(['N = ' num2str(N(Ni))])
        xlabel('Hohl et. al. quantiles')
        ylabel('Model quantiles')
    end
% end

%% MT-behavior distributions
for Ni = 1:length(N)
    figure('Name',['MT-behavior correlations, N = ' num2str(N(Ni))]);
    subplot(3,1,1)
    histogram(speedCorr_data(:,2))
    hold on
    plotVertical(hohlRquants(1:2:end));
    xlim([-0.5 0.5])
    
    subplot(3,1,2)
    histogram(rvalsModelnoise{Ni,2})
    hold on
    plotVertical(quantsRnoise(Ni,2,1:2:end));
    xlim([-0.5 0.5])
    
    subplot(3,1,3)
    histogram(rvalsModel{Ni,2})
    hold on
    plotVertical(quantsR(Ni,2,1:2:end));
    xlim([-0.5 0.5])
end

%% Quantile n model vs behavior
for quantInd = 1:length(quants)
    %quantInd = 3;
    figure('Name',['Model w/ and without noise at quantile ' num2str(quants(quantInd))])
    for szi = 1:3
        subplot(1,3,szi)
        semilogx(N,quantsR(:,szi,quantInd),'ko-')
        hold on
        semilogx(N,quantsRnoise(:,szi,quantInd),'ro-')
        plotHorizontal(hohlRquants(quantInd))
        xlabel('Neuron number')
        ylabel(['Correlation at quantile ' num2str(quants(quantInd))])
    end
end
%% Functions

%% gainSDN
function v = standardSDN(ve,w)
    %%
    v = w.^2.*ve.^2;
    
function w = fit_standardSDN(VeM,VeVAR,w0,OPTIONS)
    %%
    minimizer = @(p)( sum( sum( (VeVAR - standardSDN(VeM,p(1))).^2 ) ) );
    p = fmincon(minimizer,w0,[],[],[],[],[0,0],[Inf,Inf],[],OPTIONS);
    w = p(1);

