%% Reggie
subject = 'Reggie';
runs = 1:2;
gainNoiseFitTrainIndices = [1,3,5];
gainNoiseFitTestIndices = [2,4];
speeds = NaN;
analysisWin = [110 190];
tWin = [1 250];
Ncv = NaN;

d = load('Reggie_MultiSizePursuit.mat');

%% Main analysis
[wSDN, sigG, Gains, w_standard, wSDN_noSigG, B] = gainNoiseMultiSize(d,'runs',runs,'analysisWin',analysisWin,...
    'noSaccadeTime',0,'tWin',tWin,'plotflg',true,...
    'gainNoiseFitTrainIndices',gainNoiseFitTrainIndices,...
    'gainNoiseFitTestIndices',gainNoiseFitTestIndices,...
    'Ncv',Ncv,'alignToEye',false,...
    'speeds',speeds,'plotType','fyp');

%% Main analysis w/ eccentricity threshold
[wSDN, sigG, Gains, w_standard, wSDN_noSigG, B] = gainNoiseMultiSize(d,'runs',runs,'analysisWin',analysisWin,...
    'noSaccadeTime',0,'tWin',tWin,'plotflg',true,...
    'gainNoiseFitTrainIndices',gainNoiseFitTrainIndices,...
    'gainNoiseFitTestIndices',gainNoiseFitTestIndices,...
    'Ncv',Ncv,'eccThres',2,...
    'speeds',speeds,'plotType','fyp');

%% Perform core analysis in multiple analysis time windows
anWins(:,1) = linspace(110,190,5);
anWins(:,2) = anWins(:,1)+20;
for i = 1:size(anWins,1)
    [wSDNs(i,:), sigGs(i,:), GainsS(:,:,:,i), w_standards(:,:,i)] = gainNoiseMultiSize(d,'runs',runs,'analysisWin',anWins(i,:),...
        'noSaccadeTime',0,'tWin',tWin,'plotflg',true,...
        'gainNoiseFitTrainIndices',gainNoiseFitTrainIndices,...
        'gainNoiseFitTestIndices',gainNoiseFitTestIndices,...
        'Ncv',NaN,...
        'speeds',speeds);
end

figure
colors = [0 1 0; 1 0 0; 0 0 0];
symbols = {'^-','o-','o-'};
for di = 1:2
    subplot(2,1,di)
    for szi = 1:3
        plot(anWins(:,1)+10,squeeze(w_standards(di,szi,:)),...
            symbols{szi},'Color',colors(szi,:))
        hold on
    end
    xlabel('Time (center of bin) ms')
    ylabel('w')
end


%% Xtra
subject = 'Xtra';
runs = 1;
gainNoiseFitTrainIndices = [1,3,5];
gainNoiseFitTestIndices = [2,4];
speeds = NaN;
analysisWin = [110 190];
tWin = [1 250];
Ncv = NaN;

d = load('Xtra_MultiSizePursuit.mat');

%% Main analysis
[wSDN_x, sigG_x, Gains_x, w_standard_x, wSDN_noSigG_x, B_x] = gainNoiseMultiSize(d,'runs',runs,'analysisWin',analysisWin,...
    'noSaccadeTime',0,'tWin',tWin,'plotflg',true,...
    'gainNoiseFitTrainIndices',gainNoiseFitTrainIndices,...
    'gainNoiseFitTestIndices',gainNoiseFitTestIndices,...
    'Ncv',Ncv,'alignToEye',false,...
    'speeds',speeds,'plotType','fyp');

%% Main analysis w/ eccentricity threshold
[wSDN, sigG, Gains, w_standard, wSDN_noSigG, B] = gainNoiseMultiSize(d,'runs',runs,'analysisWin',analysisWin,...
    'noSaccadeTime',0,'tWin',tWin,'plotflg',true,...
    'gainNoiseFitTrainIndices',gainNoiseFitTrainIndices,...
    'gainNoiseFitTestIndices',gainNoiseFitTestIndices,...
    'Ncv',Ncv,'eccThres',2.5,...
    'speeds',speeds,'plotType','fyp');

%% Perform core analysis in multiple analysis time windows
anWins(:,1) = linspace(110,190,5);
anWins(:,2) = anWins(:,1)+20;
for i = 1:size(anWins,1)
    [wSDNs(i,:), sigGs(i,:), GainsS(:,:,:,i), w_standards(:,:,i)] = gainNoiseMultiSize(d,'runs',runs,'analysisWin',anWins(i,:),...
        'noSaccadeTime',0,'tWin',tWin,'plotflg',true,...
        'gainNoiseFitTrainIndices',gainNoiseFitTrainIndices,...
        'gainNoiseFitTestIndices',gainNoiseFitTestIndices,...
        'Ncv',NaN,...
        'speeds',speeds);
end

figure
colors = [0 1 0; 1 0 0; 0 0 0];
symbols = {'^-','o-','o-'};
for di = 1:2
    subplot(2,1,di)
    for szi = 1:3
        plot(anWins(:,1)+10,squeeze(w_standards(di,szi,:)),...
            symbols{szi},'Color',colors(szi,:))
        hold on
    end
    xlabel('Time (center of bin) ms')
    ylabel('w')
end

%% Gain noise w vs standard model w
figure
directions = [0 180];
ind = 1;
symbols_2 = {'+','d'; 'x','s'};
for di = 1:2
    subplot(1,2,1)
    for szi = 1:3
        w_eff = sqrt(wSDN(di).^2 + (sigG(di).^2 + sigG(di).^2*wSDN(di).^2)/B{szi}(1).^2);
        plot(w_standard(di,szi),w_eff,symbols_2{1,di},'Color',colors(szi,:))
        hold on
        
        w_eff = sqrt(wSDN_x(di).^2 + (sigG_x(di).^2 + sigG_x(di).^2*wSDN_x(di).^2)/B_x{szi}(1).^2);
        plot(w_standard_x(di,szi),w_eff,symbols_2{2,di},'Color',colors(szi,:))
    end
   axis([0 0.5 0 0.5])
    plotUnity;
    axis square
    %title(num2str(directions(di)))
    xlabel('w, flexible SDN model')
    ylabel('w, gain noise model')
    ind = ind+1;
    
    subplot(1,2,2)
    for szi = 1:3
        plot(w_standard(di,szi),wSDN_noSigG(di),symbols_2{1,di},'Color',colors(szi,:))
        hold on
        
        plot(w_standard_x(di,szi),wSDN_noSigG_x(di),symbols_2{2,di},'Color',colors(szi,:))
    end
    axis([0 0.5 0 0.5])
    plotUnity;
    axis square
    %title(num2str(directions(di)))
    xlabel('w, flexible SDN model')
    ylabel('w, rigid model')
    ind = ind+1;
end

%% Fit circuit model parameters to data
di = 1;
[gainNoise, exponent, thres, surround_weight] = ...
    gainNoiseNeuralModelParameterFit(wSDN(di),sigG(di),squeeze(Gains(1,:,di)));

%% Plot circuit model behavior
sizeProps.minEccentricity = 1;
sizeProps.maxEccentricity = 30;
sizeProps.surround_weight = surround_weight;
sizeProps.exponential = exponent;
sizeProps.threshold = thres;

% With gain noise
[~, ~, ~, ~, ~, wfit, sigGfit, Gfit] = NeuralModel_v2(...
    'thetas',0,'speeds',4:4:20,...
    'gainNoise',gainNoise,'sizeProps',sizeProps,...
    'plotMT',false,'plotDecoding',true,'plotResults',true);

%% NbackAnalysis
NbackAnalysis(d,'runs',runs,'analysisWin',analysisWin)

%% Simulation test analysis
simulationTest.On = true;
simulationTest.stdlat = [10 15 20];
simulationTest.maxlat = 50;
simulationTest.minlat = -40;
simulationTest.stdg = [0 0 0];
simulationTest.ming = 0;
[wSDNsim, sigGsim, Gainssim, w_standardsim, wSDN_noSigGsim, Bsim] = gainNoiseMultiSize(d,'runs',runs,'analysisWin',analysisWin,...
    'noSaccadeTime',0,'tWin',tWin,'plotflg',true,...
    'gainNoiseFitTrainIndices',gainNoiseFitTrainIndices,...
    'gainNoiseFitTestIndices',gainNoiseFitTestIndices,...
    'Ncv',Ncv,'simulationTest',simulationTest,...
    'speeds',speeds,'plotType','fyp');