function plotTuningProperties(modelResultsFile,varargin)
%% plotTuningProperties
%
%   plotTuningProperties(modelResultsFile)
%   
%   Plots the tuning properties of the MT population modeled in
%   modelResultsFile.
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'modelResultsFile')
addParameter(Parser,'szIndex',3)
addParameter(Parser,'restrictedDirs',true)

parse(Parser,modelResultsFile,varargin{:})

modelResultsFile = Parser.Results.modelResultsFile;
szIndex = Parser.Results.szIndex;
restrictedDirs = Parser.Results.restrictedDirs;

%% Load model results from file
load(modelResultsFile)

%% Plotting
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
%if mymakeaxisflg
    mymakeaxis(gca,'xticks',[0,10,20],'yticks',[0 10 20]);
%end

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
    %if mymakeaxisflg
        mymakeaxis(gca);
    %end
end

if saveOpts.Figs
    savefig(h,[saveOpts.location 'muVvar'])
end

%% Tuning
h = figure('Name','MT RF centers','Position',[713 316 530 420]);
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
mymakeaxis(gca);

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
end
axis square
xlim([0 24])
xlabel('Target speed (deg/s)')
ylabel('Gain')
mymakeaxis(gca,'xticks',[0,10,20]);

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

%% Noise correlations
figure('Name','Noise correlations','Position',[150 157 1089 641])
h = subplot(3,2,[1 3 5]);
imagesc(1:length(tuning.theta.pref),1:length(tuning.theta.pref),rNN{szIndex} - diag(diag(rNN{szIndex})))
colormap gray
axis square
h.YDir = 'normal';
xlabel('Neuron i')
ylabel('Neuron j')
colorbar

subplot(3,2,2)
[S1,S2] = meshgrid(tuning.speed.pref);
dS = S1-S2;
mask = logical(tril(ones(size(rNN{szIndex})),-1));
alldS = abs(dS(mask));
allrNN = rNN{szIndex}(mask);
exInds = randsample(length(alldS),2000);
hS = scatter(alldS(exInds),allrNN(exInds),20,'k');
hS.MarkerFaceColor = [0 0 0];
hold on
axis([0 150 -0.4 0.8])
plotHorizontal(0);
xlabel('\Delta Speed perference')
ylabel('Noise Correlation')
mymakeaxis(gca,'yticks',[-0.25,0,0.25,0.5,0.75,1])

subplot(3,2,4)
[D1,D2] = meshgrid(tuning.theta.pref);
dD = wrapTo180(D1-D2);
mask = logical(tril(ones(size(rNN{szIndex})),-1));
alldD = abs(dD(mask));
allrNN = rNN{szIndex}(mask);
exInds = randsample(length(alldD),2000);
hS = scatter(alldD(exInds),allrNN(exInds),20,'k');
hS.MarkerFaceColor = [0 0 0];
hold on
axis([0 180 -0.4 0.8])
xlabel('\Delta Direction perference')
ylabel('Noise Correlation')
mymakeaxis(gca,'yticks',[-0.25,0,0.25,0.5,0.75,1])

subplot(3,2,6)
[X1,X2] = meshgrid(tuning.size.x);
[Y1,Y2] = meshgrid(tuning.size.y);
dX = sqrt((X1-X2).^2 + (Y1-Y2).^2);
mask = logical(tril(ones(size(rNN{szIndex})),-1));
alldX = abs(dX(mask));
allrNN = rNN{szIndex}(mask);
exInds = randsample(length(alldX),2000);
hS = scatter(alldX(exInds),allrNN(exInds),20,'k');
hS.MarkerFaceColor = [0 0 0];
hold on
axis([0 60 -0.4 0.8])
xlabel('\Delta Position')
ylabel('Noise Correlation')
mymakeaxis(gca,'yticks',[-0.25,0,0.25,0.5,0.75,1])

%% Neuron behavior correlations
% Sort by % preferred speed
speedPercent = 100*repmat(permute(tuning.speed.pref,[3,2,1]),[size(s,1),size(s,2),1,2]) ./ ...
    repmat(s(:,:,2),[1,1,size(Rs{szIndex},3),size(Rs{szIndex},4)]);

dirDiff = repmat(s(:,:,1),[1,1,size(Rs{szIndex},3),size(Rs{szIndex},4)]) - ...
    repmat(permute(tuning.theta.pref,[3,2,1]),[size(s,1),size(s,2),1,2]);

figure('Name','Neuron-estimate correlations');
dispSps = 1:size(Rs{szIndex},2);
dispDirs = 1:size(Rs{szIndex},1);
thres = 0.05;
for i = 1:2
    subplot(2,2,i)
    Dtemp = dirDiff(dispDirs,dispSps,:,i);
    if restrictedDirs
        dispDirs2 = squeeze(Dtemp(:,1,:) < 45 & Dtemp(:,1,:) > -45 | Dtemp(:,1,:) < -135 | Dtemp(:,1,:) > 135);
    else
        dispDirs2 = true(size(Rs{szIndex},3),1);
    end
    Dtemp(Dtemp < -180) = 360 + Dtemp(Dtemp < -180);
    Dtemp(Dtemp > 180) = 360 - Dtemp(Dtemp > 180);
    Dtemp = Dtemp(:,:,dispDirs2);
    Rtemp = Rs{szIndex}(dispDirs,dispSps,dispDirs2,i);
    Ptemp = PVals{szIndex}(dispDirs,dispSps,dispDirs2,i);
    %         tempInds = randsample(length(Dtemp),379);
    %         Dtemp = Dtemp(tempInds);
    %         Rtemp = Rtemp(tempInds);
    %         Ptemp = Ptemp(tempInds);
    h = scatter(abs(Dtemp(:)),Rtemp(:));
    set(h,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
    hold on
    hsig = scatter(abs(Dtemp(Ptemp < thres)),Rtemp(Ptemp < thres));
    set(hsig,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
    axis([0 180 -0.6 0.6])
    plotHorizontal(0.5);
    plotHorizontal(0);
    plotHorizontal(-0.5);
    hax = gca;
    hax.TickDir = 'out';
    xlabel('Difference from preferred direction ')
    if i == 1
        ylabel('Neuron-eye direction correlation')
    else
        ylabel('Neuron-eye speed correlation')
    end
    if mymakeaxisflg
        mymakeaxis(gca)
    end
    
end
for i = 1:2
    subplot(2,2,i+2)
    Dtemp = dirDiff(dispDirs,dispSps,:,i);
    if restrictedDirs
        dispDirs2 = squeeze(Dtemp(:,1,:) < 45 & Dtemp(:,1,:) > -45 | Dtemp(:,1,:) < -135 | Dtemp(:,1,:) > 135);
    else
        dispDirs2 = true(size(Rs{szIndex},3),1);
    end
    dispDirs2 = squeeze(Dtemp(:,1,:) < 45 & Dtemp(:,1,:) > -45 | Dtemp(:,1,:) < -135 | Dtemp(:,1,:) > 135);
    Stemp = -log2(speedPercent(dispDirs,dispSps,dispDirs2,i)/100);
    Rtemp = Rs{szIndex}(dispDirs,dispSps,dispDirs2,i);
    Ptemp = PVals{szIndex}(dispDirs,dispSps,dispDirs2,i);
    %         tempInds = randsample(length(Stemp),379);
    %         Stemp = Stemp(tempInds);
    %         Rtemp = Rtemp(tempInds);
    %         Ptemp = Ptemp(tempInds);
    h = scatter(Stemp(:),Rtemp(:));
    set(h,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
    hold on
    hsig = scatter(Stemp(Ptemp < thres),Rtemp(Ptemp < thres));
    set(hsig,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
    %         set(gca,'XScale','log')
    ax = axis;
    axis([ax(1) ax(2) -0.6 0.6])
    plotHorizontal(0.5);
    plotHorizontal(0);
    plotHorizontal(-0.5);
    axis square
    hax = gca;
    hax.TickDir = 'out';
    hax.YTick = [-0.5,0,0.5];
    hax.XTick = [1,10,100,1000];
    xlabel('Percent perferred speed')
    if i == 1
        ylabel('Neuron-eye direction correlation')
    else
        ylabel('Neuron-eye speed correlation')
    end
    if mymakeaxisflg
        mymakeaxis(gca)
    end
end

%% Functions

%% gainSDN
function v = gainSDN(ve,vs,w,sigG)
    %%
    v = w.^2.*ve.^2 + (sigG.^2 + sigG.^2.*w.^2).*vs.^2;