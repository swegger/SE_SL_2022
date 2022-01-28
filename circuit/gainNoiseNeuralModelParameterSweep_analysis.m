function gainNoiseNeuralModelParameterSweep_analysis(dataDirectory,varargin)
%% gainNoiseNeuralModelParameterSweep_analysis
%
%   gainNoiseNeuralModelParameterSweep_analysis()
%   
%   Performs analysis of neural model output for default results directory
%
%%

%% Defaults
dataDate_default = '20210602';


%% Parse inputs
Parser = inputParser;

addRequired(Parser,'dataDirectory')
addParameter(Parser,'dataDate',dataDate_default)
addParameter(Parser,'calcNeuronPursuitCorrelation',false)
addParameter(Parser,'gainNoise',0.4)

parse(Parser,dataDirectory,varargin{:})

dataDirectory = Parser.Results.dataDirectory;
dataDate = Parser.Results.dataDate;
calcNeuronPursuitCorrelation = Parser.Results.calcNeuronPursuitCorrelation;
gainNoise = Parser.Results.gainNoise;

%% Find files corresponding to desired date
potentialFiles = dir([dataDirectory '/*' dataDate '*.mat']);
files = potentialFiles;%(3:end);

%% For each file, load in parameter values and results
params = nan(length(files),4);
ws = nan(length(files),1);
sigGs = nan(length(files),1);
Gs = nan(length(files),3);
w_standard = nan(length(files),3);
OPTIONS = optimset('Display','off');

for filei = 182:length(files)
    disp(['File ' num2str(filei) ' of ' num2str(length(files))])
    results = load(files(filei).name,'sizeProps','w','sigG','G','gainNoise',...
        'VeM','VeVAR','decoderAlgorithm');
    if isfield(results,'decoderAlgorithm')
        decoderAlgorithm = results.decoderAlgorithm;
    else
        decoderAlgorithm = 'bioRxiv2020';
    end
    params(filei,1) = results.sizeProps.exponential;
    params(filei,2) = results.sizeProps.threshold;
    params(filei,3) = results.sizeProps.surround_weight;
    params(filei,4) = results.gainNoise;
    ws(filei) = results.w;
    sigGs(filei) = results.sigG;
    Gs(filei,:) = results.G;
    VeM = permute(results.VeM,[3,2,1]);
    VeVAR = permute(results.VeVAR,[3,2,1]);
    
    for szi = 1:3
        w_standard(filei,szi) = fit_standardSDN(...
            VeM(szi,:),VeVAR(szi,:),0.1,OPTIONS); 
    end
    
    if calcNeuronPursuitCorrelation
        szind = 2;
        results = load(files(filei).name,'e','n','s','tuning','epsilon','normalizer');
        n = results.n;
        e = results.e{szind};
        s = results.s;
        
        % z-score
        nM = mean(n{szind},3);
        pM = mean(e,3);
        nVar = var(n{szind},[],3);
        pVar = var(e,[],3);
        
        nZ = (n{szind} - repmat(nM,[1,1,size(n{szind},3),1]))./repmat(sqrt(nVar),[1,1,size(n{szind},3),1]);
        pZ = (e - repmat(pM,[1,1,size(n{szind},3),1]))./repmat(sqrt(pVar),[1,1,size(n{szind},3),1]);
        
        % Correlate responses
        Rs = nan(size(s,1),size(s,2),size(n{szind},4),2);
        PVals = nan(size(s,1),size(s,2),size(n{szind},4),2);
        for thetai = 1:size(s,1)
            disp(thetai/size(s,1))
            for speedi = 1:size(s,2)
                for neuroni = 1:size(nZ,4)
                    [rtemp,ptemp] = corrcoef(...
                        [permute(nZ(thetai,speedi,:,neuroni),[3,4,1,2]),...
                        permute(pZ(thetai,speedi,:,:),[3,4,2,1])]);
                    
                    Rs(thetai,speedi,neuroni,:) = rtemp(1,2:3);
                    PVals(thetai,speedi,neuroni,:) = ptemp(1,2:3);
                end
            end
        end
        
        % Max corr
        maxCorr(filei,1) = max(Rs(:));
        
        % Now add noise
        for szi = 1:3
            [etemp, ~, ~, RsTemp{szi}, ~] = DecodeMT(n{szi},results.tuning,s,...
                'gainNoise',gainNoise,'epsilon',results.epsilon,'b',results.normalizer,...
                'mymakeaxisflg',false,'plotflg',false,'decoderAlgorithm',decoderAlgorithm,...
                'motorNoise',0);
    
            eBar = mean(etemp,3);
            eVar = var(etemp,1,3);
            
            VeMnoise(:,:,szi) = squeeze(eBar(:,:,:,2));
            VeVARnoise(:,:,szi) = squeeze(eVar(:,:,:,2));
            
            w_standard_noise(filei,szi) = fit_standardSDN(...
                VeMnoise(:,:,szi),VeVARnoise(:,:,szi),0.1,OPTIONS);
        end
        maxCorr_noise(filei,1) = max(RsTemp{szind}(:));
    end
end

%% Load in experimental fit
sigGobs = 0.1;

%%
figure('Name','Parameter effect on sigG')
subplot(2,3,1)
% scatter3(params(:,1),params(:,2),params(:,3),10,sigGs.^2*1e9)
% xlabel('Exp')
% ylabel('Thres')
% zlabel('SS')
% colorbar
% cax = caxis;
subplot(2,3,2)
exps = unique(params(:,1));
thres = unique(params(:,2));
sss = unique(params(:,3));
for expi = 1:length(exps)
    for thresi = 1:length(thres)
        msigG(expi,thresi) = mean(sigGs(params(:,1) == exps(expi) & params(:,2) == thres(thresi)));
        p1(expi,thresi) = exps(expi);
        p2(expi,thresi) = thres(thresi);
    end
end
scatter(p1(:),p2(:),10,msigG(:).^2*1e9)
% caxis(cax)
%scatter(params(:,1),params(:,2),10,sigGs)
xlabel('Exp')
ylabel('Thres')
subplot(2,3,4)
clear msigG p1 p2
for expi = 1:length(exps)
    for ssi = 1:length(sss)
        msigG(expi,ssi) = mean(sigGs(params(:,1) == exps(expi) & params(:,3) == sss(ssi)));
        p1(expi,ssi) = exps(expi);
        p2(expi,ssi) = sss(ssi);
    end
end
scatter(p1(:),p2(:),10,msigG(:).^2*1e9)
% caxis(cax)
%scatter(params(:,1),params(:,3),10,sigGs)
xlabel('Exp')
ylabel('SS')
subplot(2,3,5)
clear msigG p1 p2
for thresi = 1:length(thres)
    for ssi = 1:length(sss)
        msigG(thresi,ssi) = mean(sigGs(params(:,2) == thres(thresi) & params(:,3) == sss(ssi)));
        p1(thresi,ssi) = thres(thresi);
        p2(thresi,ssi) = sss(ssi);
    end
end
scatter(p1(:),p2(:),10,msigG(:).^2*1e9)
% caxis(cax)
%scatter(params(:,2),params(:,3),10,sigGs)
xlabel('Thres')
ylabel('SS')

subplot(2,3,3)
histogram((sigGs.^2)*1e9,20,'Normalization','probability')
xlabel('\sigma^2')
ylabel('Relative frequency')
% plotVertical(mean(sigGs));
% hold on
% plotVertical(sigGobs);

%%
figure('Name','Parameter effect on w')
subplot(2,3,1)
scatter3(params(:,1),params(:,2),params(:,3),10,ws)
xlabel('Exp')
ylabel('Thres')
zlabel('SS')
subplot(2,3,2)
for expi = 1:length(exps)
    for thresi = 1:length(thres)
        mw = mean(ws(params(:,1) == exps(expi) & params(:,2) == thres(thresi)));
        scatter(exps(expi),thres(thresi),10,mw)
        hold on
    end
end
%scatter(params(:,1),params(:,2),10,ws)
xlabel('Exp')
ylabel('Thres')
subplot(2,3,4)
for expi = 1:length(exps)
    for ssi = 1:length(sss)
        mw = mean(ws(params(:,1) == exps(expi) & params(:,3) == sss(ssi)));
        scatter(exps(expi),sss(ssi),10,mw)
        hold on
    end
end
%scatter(params(:,1),params(:,3),10,ws)
xlabel('Exp')
ylabel('SS')
subplot(2,3,5)
for thresi = 1:length(thres)
    for ssi = 1:length(sss)
        mw = mean(ws(params(:,2) == thres(thresi) & params(:,3) == sss(ssi)));
        p1(thresi,ssi) = thres(thresi);
        p2(thresi,ssi) = sss(ssi);
    end
end
%scatter(params(:,2),params(:,3),10,ws)
xlabel('Thres')
ylabel('SS')

subplot(2,3,3)
histogram(ws,linspace(0,0.2,100))
% plotVertical(mean(ws));
hold on
% plotVertical(sigGobs);

%%
figure('Name','Parameter effect on Gain')
subplot(1,3,1)
plot3(params(:,1),params(:,2),Gs,'o')
xlabel('Exp')
ylabel('Thres')
subplot(1,3,2)
plot3(params(:,1),params(:,3),Gs,'o')
xlabel('Exp')
ylabel('SS')
subplot(1,3,3)
plot3(params(:,2),params(:,3),Gs,'o')
xlabel('Thres')
ylabel('SS')

%% Parameter effect on maximum neuron-behavior correlation
if calcNeuronPursuitCorrelation
    figure('Name','Parameter effect on maximum neuron-behavior correlation')
    subplot(2,3,1)
    % scatter3(params(:,1),params(:,2),params(:,3),10,sigGs.^2*1e9)
    % xlabel('Exp')
    % ylabel('Thres')
    % zlabel('SS')
    % colorbar
    % cax = caxis;
    subplot(2,3,2)
    exps = unique(params(:,1));
    thres = unique(params(:,2));
    sss = unique(params(:,3));
    for expi = 1:length(exps)
        for thresi = 1:length(thres)
            mR(expi,thresi) = mean(maxCorr(params(:,1) == exps(expi) & params(:,2) == thres(thresi)));
            p1(expi,thresi) = exps(expi);
            p2(expi,thresi) = thres(thresi);
        end
    end
    scatter(p1(:),p2(:),10,mR(:))
    cax = [0.53 0.55];
    caxis(cax)
    %scatter(params(:,1),params(:,2),10,sigGs)
    xlabel('Exp')
    ylabel('Thres')
    subplot(2,3,4)
    clear msigG p1 p2
    for expi = 1:length(exps)
        for ssi = 1:length(sss)
            mR(expi,ssi) = mean(maxCorr(params(:,1) == exps(expi) & params(:,3) == sss(ssi)));
            p1(expi,ssi) = exps(expi);
            p2(expi,ssi) = sss(ssi);
        end
    end
    scatter(p1(:),p2(:),10,mR(:))
    caxis(cax)
    %scatter(params(:,1),params(:,3),10,sigGs)
    xlabel('Exp')
    ylabel('SS')
    subplot(2,3,5)
    clear msigG p1 p2
    for thresi = 1:length(thres)
        for ssi = 1:length(sss)
            mR(thresi,ssi) = mean(maxCorr(params(:,2) == thres(thresi) & params(:,3) == sss(ssi)));
            p1(thresi,ssi) = thres(thresi);
            p2(thresi,ssi) = sss(ssi);
        end
    end
    scatter(p1(:),p2(:),10,mR(:))
    caxis(cax)
    %scatter(params(:,2),params(:,3),10,sigGs)
    xlabel('Thres')
    ylabel('SS')
    
    subplot(2,3,3)
    histogram(maxCorr,linspace(0,0.6,50),'Normalization','probability')
    hold on
    histogram(maxCorr_noise,linspace(0,0.6,50),'Normalization','probability')
    plotVertical(0.3378);
    xlabel('Maximum MT-behavior correalation')
    ylabel('Relative frequency')
    legend('\sigma_g^2 = 0',['\sigma_g^2 = ' num2str(gainNoise.^2)])
    
end

%% Parameter effect on change in w_standard
figure('Name','Parameter effect delta w_standard')
subplot(2,3,1)
% scatter3(params(:,1),params(:,2),params(:,3),10,sigGs.^2*1e9)
% xlabel('Exp')
% ylabel('Thres')
% zlabel('SS')
% colorbar
% cax = caxis;
subplot(2,3,2)
exps = unique(params(:,1));
thres = unique(params(:,2));
sss = unique(params(:,3));
deltaW = w_standard(:,1) - w_standard(:,end);
deltaW_noise = w_standard_noise(:,1) - w_standard_noise(:,end);
for expi = 1:length(exps)
    for thresi = 1:length(thres)
        mdW(expi,thresi) = mean(deltaW(params(:,1) == exps(expi) & params(:,2) == thres(thresi)));
        p1(expi,thresi) = exps(expi);
        p2(expi,thresi) = thres(thresi);
    end
end
scatter(p1(:),p2(:),10,mdW(:))
xlabel('Exp')
ylabel('Thres')
subplot(2,3,4)
clear msigG p1 p2
for expi = 1:length(exps)
    for ssi = 1:length(sss)
        mdW(expi,ssi) = mean(deltaW(params(:,1) == exps(expi) & params(:,3) == sss(ssi)));
        p1(expi,ssi) = exps(expi);
        p2(expi,ssi) = sss(ssi);
    end
end
scatter(p1(:),p2(:),10,mdW(:))
xlabel('Exp')
ylabel('SS')
subplot(2,3,5)
clear msigG p1 p2
for thresi = 1:length(thres)
    for ssi = 1:length(sss)
        mdW(thresi,ssi) = mean(deltaW(params(:,2) == thres(thresi) & params(:,3) == sss(ssi)));
        p1(thresi,ssi) = thres(thresi);
        p2(thresi,ssi) = sss(ssi);
    end
end
scatter(p1(:),p2(:),10,mdW(:))
xlabel('Thres')
ylabel('SS')

subplot(2,3,3)
deltaWbins = 10*linspace(0,0.25,59);
histogram(10*deltaW,deltaWbins)
hold on
histogram(10*deltaW_noise,deltaWbins)
plotVertical(0.0698*10)
plotVertical(0.1034*10)
xlabel('w_{2 deg} - w_{20 deg}')
ylabel('Count')
legend('\sigma_G^2 = 0',['\sigma_g^2 = ' num2str(gainNoise.^2)])

%% Parameter effect on change in w_standard
figure('Name','Parameter effect delta w_standard')
subplot(2,3,1)
subplot(2,3,2)
exps = unique(params(:,1));
thres = unique(params(:,2));
sss = unique(params(:,3));
deltaW = w_standard(:,2) - w_standard(:,end);
deltaW_noise = w_standard_noise(:,2) - w_standard_noise(:,end);
for expi = 1:length(exps)
    for thresi = 1:length(thres)
        mdW(expi,thresi) = mean(deltaW(params(:,1) == exps(expi) & params(:,2) == thres(thresi)));
        p1(expi,thresi) = exps(expi);
        p2(expi,thresi) = thres(thresi);
    end
end
scatter(p1(:),p2(:),10,mdW(:))
xlabel('Exp')
ylabel('Thres')
subplot(2,3,4)
clear msigG p1 p2
for expi = 1:length(exps)
    for ssi = 1:length(sss)
        mdW(expi,ssi) = mean(deltaW(params(:,1) == exps(expi) & params(:,3) == sss(ssi)));
        p1(expi,ssi) = exps(expi);
        p2(expi,ssi) = sss(ssi);
    end
end
scatter(p1(:),p2(:),10,mdW(:))
xlabel('Exp')
ylabel('SS')
subplot(2,3,5)
clear msigG p1 p2
for thresi = 1:length(thres)
    for ssi = 1:length(sss)
        mdW(thresi,ssi) = mean(deltaW(params(:,2) == thres(thresi) & params(:,3) == sss(ssi)));
        p1(thresi,ssi) = thres(thresi);
        p2(thresi,ssi) = sss(ssi);
    end
end
scatter(p1(:),p2(:),10,mdW(:))
xlabel('Thres')
ylabel('SS')

subplot(2,3,3)
deltaWbins = 10*linspace(0,0.25,59);
histogram(10*deltaW,deltaWbins)
hold on
histogram(10*deltaW_noise,deltaWbins)
plotVertical(0.0317*10)
plotVertical(0.0522*10)
xlabel('w_{6 deg} - w_{20 deg}')
ylabel('Count')
legend('\sigma_G^2 = 0',['\sigma_g^2 = ' num2str(gainNoise.^2)])


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