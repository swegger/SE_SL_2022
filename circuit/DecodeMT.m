function [e,gain,figureHandles, Rs, PVals] = DecodeMT(n,tuning,s,varargin)
%% DecodeMT
%
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'n')
addRequired(Parser,'tuning')
addRequired(Parser,'s')
addParameter(Parser,'epsilon',250)
addParameter(Parser,'gainNoise',0.1)
addParameter(Parser,'motorNoise',0)
addParameter(Parser,'k',0.65)
addParameter(Parser,'b',[0, 0])
addParameter(Parser,'gainSDN',false)
addParameter(Parser,'restrictedDirs',true)
addParameter(Parser,'pool',NaN)
addParameter(Parser,'decoderAlgorithm','bioRxiv2020')
addParameter(Parser,'plotflg',true)
addParameter(Parser,'mymakeaxisflg',false)

parse(Parser,n,tuning,s,varargin{:})

n = Parser.Results.n;
tuning = Parser.Results.tuning;
s = Parser.Results.s;
epsilon = Parser.Results.epsilon;
gainNoise = Parser.Results.gainNoise;
motorNoise = Parser.Results.motorNoise;
k = Parser.Results.k;
b = Parser.Results.b;
gainSDN = Parser.Results.gainSDN;
restrictedDirs = Parser.Results.restrictedDirs;
pool = Parser.Results.pool;
decoderAlgorithm = Parser.Results.decoderAlgorithm;
plotflg = Parser.Results.plotflg;
mymakeaxisflg = Parser.Results.mymakeaxisflg;

if any(isnan(pool))
    pool = true(size(n,4),1);
    dpool = pool;
else
    dpool = ~pool;
end

%% Decode trial-by-trial population responses
n(n<0) = 0;
switch decoderAlgorithm
    case 'bioRxiv2020'
        numerator = vectorAverage(n,tuning,pool); 
        denominator = (epsilon + sum(n(:,:,:,dpool),4));
        gain = gainFunction(n,tuning,0)./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = repmat(gain,[1,1,1,2]).*numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = 2.^(sqrt(sqrt(sum((vA).^2,4))));
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    case 'ignoreDirection'
        for thetai = 1:size(n,1)
            for speedi = 1:size(n,2)
                numerator(thetai,speedi,:,1) = log2(tuning.speed.pref(pool))' ...
                    * permute(n(thetai,speedi,:,pool),[4,3,1,2]);
            end
        end
        numerator = repmat(numerator,[1,1,1,2]);
        denominator = (epsilon + sum(n(:,:,:,dpool),4));
        gain = gainFunction(n,tuning,0)./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = repmat(gain,[1,1,1,2]).*2.^(numerator./repmat(denominator,[1,1,1,2]));
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = vA(:,:,:,1);
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    case 'huangAndLisberger'
        numerator = vectorAverage(n,tuning,pool); 
        denominator = epsilon + sqrt(sum(huangAndLisberger(n,tuning,dpool).^2,4));
        
        gain = gainFunction(n,tuning,0)./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = repmat(gain,[1,1,1,2]).*numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = 2.^(sqrt(sqrt(sum((vA).^2,4))));
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    case 'Darlington2018'
        numerator = vectorAverage(n,tuning,pool);
        denominator = (epsilon + sum(n(:,:,:,dpool),4));
        
        gain = (gainFunction(n,tuning,0) + b(1))./(sum(n,4) + b(2));
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = repmat(gain,[1,1,1,2]).*numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = 2.^(sqrt(sqrt(sum((vA).^2,4))));
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    case 'gainConstant'
        numerator = vectorAverage(n,tuning,pool);
        denominator = epsilon + sqrt(sum(huangAndLisberger(n,tuning,dpool).^2,4));
        
        gain = 1./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = repmat(gain,[1,1,1,2]).*numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = 2.^(sqrt(sum((vA).^2,4)));
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
    
    case 'singleSqrt'
        numerator = vectorAverage(n,tuning,pool); 
        denominator = (epsilon + sum(n(:,:,:,dpool),4));
        gain = gainFunction(n,tuning,0)./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = repmat(gain,[1,1,1,2]).*numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = 2.^(sqrt(sum((vA).^2,4)));
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    case 'simpleSumVA^2'
        numerator = vectorAverage(n,tuning,pool); 
        denominator = (epsilon + sum(n(:,:,:,dpool),4));
        gain = gainFunction(n,tuning,0)./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = repmat(gain,[1,1,1,2]).*numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = sum((vA).^2,4);
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    case 'simpleSumVA'
        numerator = vectorAverage(n,tuning,pool); 
        denominator = (epsilon + sum(n(:,:,:,dpool),4));
        gain = gainFunction(n,tuning,0)./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = repmat(gain,[1,1,1,2]).*numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = sum((vA),4);
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    case 'g*log2shat'
        numerator = vectorAverage(n,tuning,pool); 
        denominator = (epsilon + sum(n(:,:,:,dpool),4));
        gain = gainFunction(n,tuning,0)./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = gain.*sqrt(sum(vA.^2,4));
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    case 'g*2^shat'
        numerator = vectorAverage(n,tuning,pool); 
        denominator = (epsilon + sum(n(:,:,:,dpool),4));
        gain = gainFunction(n,tuning,0)./b(1);
        
        if gainSDN
            gain = gain + ...
                gainNoise*gain.*randn([size(n,1),size(n,2),size(n,3)]);
        else
            gain = gain + ...
                gainNoise*randn([size(n,1),size(n,2),size(n,3)]);
        end
        
        vA = numerator./repmat(denominator,[1,1,1,2]);
        temp = atan2d(vA(:,:,:,2),vA(:,:,:,1));
        e(:,:,:,1) = temp;
        e(:,:,:,2) = gain.*2.^sqrt(sum(vA.^2,4));
        e(:,:,:,2) = e(:,:,:,2) + ...
            e(:,:,:,2).*motorNoise.*randn([size(vA,1),size(vA,2),size(vA,3)]);
        
    otherwise
        error('decoderAlgorithm not recognized!')
end

%% Correlation with MT
% z-score
nM = mean(n,3);
pM = mean(e,3);
nVar = var(n,[],3);
pVar = var(e,[],3);

nZ = (n - repmat(nM,[1,1,size(n,3),1]))./repmat(sqrt(nVar),[1,1,size(n,3),1]);
pZ = (e - repmat(pM,[1,1,size(n,3),1]))./repmat(sqrt(pVar),[1,1,size(n,3),1]);

% Correlate responses
Rs = nan(size(s,1),size(s,2),size(n,4),2);
PVals = nan(size(s,1),size(s,2),size(n,4),2);
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

% Sort by % preferred speed
speedPercent = 100*repmat(permute(tuning.speed.pref,[3,2,1]),[size(s,1),size(s,2),1,2]) ./ ...
    repmat(s(:,:,2),[1,1,size(Rs,3),size(Rs,4)]);

dirDiff = repmat(s(:,:,1),[1,1,size(Rs,3),size(Rs,4)]) - ...
    repmat(permute(tuning.theta.pref,[3,2,1]),[size(s,1),size(s,2),1,2]);

%% Plotting
if plotflg
    figureHandles.decodingPerformance = ...
        figure('Name','Decoding performance','Position',[112 378 1197 420]);
    subplot(1,3,1)
    for triali = 1:size(e,3)
        quiver(vA(:,:,triali,1),vA(:,:,triali,2),'Color',[0.6 0.6 0.6])
        hold on
    end
    quiver(cosd(s(:,:,1)).*log2(s(:,:,2)),sind(s(:,:,1)).*log2(s(:,:,2)),'k','LineWidth',2)
    axis square
    
    subplot(1,3,2)
    for triali = 1:size(e,3)
        temps = s(:,:,2);
        tempe = squeeze(e(:,:,triali,2));
        h = scatter(temps(:),tempe(:),20,[0 0 0]);
        h.MarkerFaceColor = [0 0 0];
        h.MarkerFaceAlpha = 0.01;
        h.MarkerEdgeAlpha = 0.01;
        hold on
    end
    set(gca,'XScale','log','YScale','log')
    axis square
    plotUnity;
    xlabel('Stimulus speed')
    ylabel('Speed estimate')
    if mymakeaxisflg
        mymakeaxis(gca)
    end
    
    subplot(1,3,3)
    for triali = 1:size(e,3)
        temps = s(:,:,1);
        tempe = squeeze(e(:,:,triali,1));
        h = scatter(temps(:),tempe(:),20,[0 0 0]);
        h.MarkerFaceColor = [0 0 0];
        h.MarkerFaceAlpha = 0.01;
        h.MarkerEdgeAlpha = 0.01;
        hold on
    end
    axis square
    plotUnity;
    xlabel('Stimulus direction')
    ylabel('Direction estimate')
    if mymakeaxisflg
        mymakeaxis(gca)
    end
    
    %%
    figureHandles.neuon_est_corr = ...
        figure('Name','Neuron-estimate correlations','Position',[1345 264 1176 1058]);
    dispSps = 1:size(Rs,2);
    dispDirs = 1:size(Rs,1);
    thres = 0.05;
    for i = 1:2   
        subplot(2,2,i)
        Dtemp = dirDiff(dispDirs,dispSps,:,i);
        if restrictedDirs
            dispDirs2 = squeeze(Dtemp(:,1,:) < 45 & Dtemp(:,1,:) > -45 | Dtemp(:,1,:) < -135 | Dtemp(:,1,:) > 135);
        else
            dispDirs2 = true(size(Rs,3),1);
        end
        Dtemp(Dtemp < -180) = 360 + Dtemp(Dtemp < -180);
        Dtemp(Dtemp > 180) = 360 - Dtemp(Dtemp > 180);
        Dtemp = Dtemp(:,:,dispDirs2);
        Rtemp = Rs(dispDirs,dispSps,dispDirs2,i);
        Ptemp = PVals(dispDirs,dispSps,dispDirs2,i);
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
            dispDirs2 = true(size(Rs,3),1);
        end
        dispDirs2 = squeeze(Dtemp(:,1,:) < 45 & Dtemp(:,1,:) > -45 | Dtemp(:,1,:) < -135 | Dtemp(:,1,:) > 135);
        Stemp = -log2(speedPercent(dispDirs,dispSps,dispDirs2,i)/100);
        Rtemp = Rs(dispDirs,dispSps,dispDirs2,i);
        Ptemp = PVals(dispDirs,dispSps,dispDirs2,i);
        h = scatter(Stemp(:),Rtemp(:));
        set(h,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0])
        hold on
        hsig = scatter(Stemp(Ptemp < thres),Rtemp(Ptemp < thres));
        set(hsig,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
        ax = axis;
        axis([ax(1) ax(2) -0.6 0.6])
        plotHorizontal(0.5);
        plotHorizontal(0);
        plotHorizontal(-0.5);
        axis square
        hax = gca;
        hax.TickDir = 'out';
        hax.YTick = [-0.5,0,0.5];
        hax.XTick = [-5,1,5];
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
    
else
    figureHandles = NaN;
end
    


%% Functions

%% vectorAverage
function out = vectorAverage(n,tuning,pool)
for thetai = 1:size(n,1)
    for speedi = 1:size(n,2)
        out(thetai,speedi,:,1) = cosd(tuning.theta.pref(pool))'.*log2(tuning.speed.pref(pool))' ...
            * permute(n(thetai,speedi,:,pool),[4,3,1,2]);
        out(thetai,speedi,:,2) = sind(tuning.theta.pref(pool))'.*log2(tuning.speed.pref(pool))' ...
            * permute(n(thetai,speedi,:,pool),[4,3,1,2]);
    end
end

%% gainFunction
function out = gainFunction(n,tuning,gainNoise)
%%
for thetai = 1:size(n,1)
    for speedi = 1:size(n,2)
        out(thetai,speedi,:) = log2(tuning.speed.pref)' ...
            * permute(n(thetai,speedi,:,:),[4,3,1,2]) + ...
            gainNoise*randn(size(n,3),1)';
        
    end
end

%% huang&Lisberger
function out = huangAndLisberger(n,tuning,dpool)
%%
for thetai = 1:size(n,1)
    for speedi = 1:size(n,2)
        out(thetai,speedi,:,1) = cosd(tuning.theta.pref(dpool))'*...
            permute(n(thetai,speedi,:,dpool),[4,3,1,2]);
        out(thetai,speedi,:,2) = sind(tuning.theta.pref(dpool))'*...
            permute(n(thetai,speedi,:,dpool),[4,3,1,2]);
    end
end
