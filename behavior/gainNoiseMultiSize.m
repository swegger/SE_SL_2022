function [wSDN, sigG, Gains, w_standard, wSDN_noSigG, B] = gainNoiseMultiSize(d,varargin)
%% gainNoiseMultiSize
%
%   gainNoiseMultiSize()
%
%
%%

%% Defaults
simulationTest_default.On = false;
saveTable_default.On = false;

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'d')
addParameter(Parser,'runs',NaN)
addParameter(Parser,'tWin',[50, 400])
addParameter(Parser,'analysisWin',[30, 80])
addParameter(Parser,'speeds',NaN)
addParameter(Parser,'directions',NaN)
addParameter(Parser,'sizes',NaN)
addParameter(Parser,'noSaccadeTime',150)
addParameter(Parser,'smoothWin',20)
addParameter(Parser,'gainNoiseFitTrainIndices',NaN)
addParameter(Parser,'gainNoiseFitTestIndices',NaN)
addParameter(Parser,'latencyFract',0.2)
addParameter(Parser,'startStop',[80 160])
addParameter(Parser,'eccThres',Inf)
addParameter(Parser,'simulationTest',simulationTest_default)
addParameter(Parser,'findLeeLatency',true)
addParameter(Parser,'outlierReject',true)
addParameter(Parser,'Ncv',1000)
addParameter(Parser,'alignToEye',false)
addParameter(Parser,'plotflg',false)
addParameter(Parser,'plotType','Default')
addParameter(Parser,'saveTable',saveTable_default)

parse(Parser,d,varargin{:})

d = Parser.Results.d;
runs = Parser.Results.runs;
tWin = Parser.Results.tWin;
analysisWin = Parser.Results.analysisWin;
speeds = Parser.Results.speeds;
directions = Parser.Results.directions;
sizes = Parser.Results.sizes;
noSaccadeTime = Parser.Results.noSaccadeTime;
smoothWin = Parser.Results.smoothWin;
gainNoiseFitTrainIndices = Parser.Results.gainNoiseFitTrainIndices;
gainNoiseFitTestIndices = Parser.Results.gainNoiseFitTestIndices;
latencyFract = Parser.Results.latencyFract;
startStop = Parser.Results.startStop;
eccThres = Parser.Results.eccThres;
simulationTest = Parser.Results.simulationTest;
findLeeLatency = Parser.Results.findLeeLatency;
outlierReject = Parser.Results.outlierReject;
Ncv = Parser.Results.Ncv;
alignToEye = Parser.Results.alignToEye;
plotflg = Parser.Results.plotflg;
plotType = Parser.Results.plotType;
saveTable = Parser.Results.saveTable;

if isnan(runs)
    runs = 1:length(d.meta);
end

OPTIONS = optimset('Display','off');

%% Unpack data
s = horzcat(d.target.s{runs});
dir = horzcat(d.target.d{runs});
sz = horzcat(d.target.sz{runs});
loc = horzcat(d.target.x{runs});

hvelocity = horzcat(d.eye.hvel{runs});
hposition = horzcat(d.eye.hpos{runs});
vposition = horzcat(d.eye.vpos{runs});
saccades = horzcat(d.eye.saccades{runs});
saccOn = [zeros(1,size(saccades,2)); diff(saccades,1,1)] > 0;
saccOff = [zeros(1,size(saccades,2)); diff(saccades,1,1)] < 0;
saccMarkers = false(size(saccOn));
for ti = 1:size(saccOn,2)
    OnInds = find(saccOn(:,ti));
    OffInds = find(saccOff(:,ti));
    for indi = 1:length(OnInds)
        if length(OffInds) >= indi
            saccMarkers(OnInds(indi)+...
                round((OffInds(indi)-OnInds(indi))/2),ti) = true;
        else
            saccMarkers(OnInds(indi),ti) = true;
        end
    end
end

trialNumbers = 1:size(s,2);

%% Determine speeds, directions, sizes
if isnan(speeds)
    speeds = unique(s);
end
if isnan(directions)
    directions = unique(dir);
end
if isnan(sizes)
    sizes = unique(sz);
end
colors = projectColorMaps('speeds','samples',1:length(speeds),...
    'sampleDepth',length(speeds));
szcolors = [0 1 0;...
                1 0 0;...
                0 0 0];

%% For each combination, generate mask and sort data
m_hvel = nan(length(tWin(1):tWin(2)),length(directions),length(sizes),length(speeds));
std_hvel = nan(length(tWin(1):tWin(2)),length(directions),length(sizes),length(speeds));

for di = 1:length(directions)
    for szi = 1:length(sizes)
        for si = 1:length(speeds)
            masks(:,di,szi,si) = ( ...
                s == speeds(si) & ...
                dir == directions(di) & ...
                sz == sizes(szi) & ...
                loc == 0 ...
                );
                
                eccMask = sqrt(hposition(1,:).^2+vposition(1,:).^2)' <= eccThres;
                sacMask = ~any(saccMarkers(tWin(1):tWin(2),:),1)';
                
                htemp = cosd(directions(di))*hvelocity(tWin(1)+1:tWin(2)+1,:) ...
                    >= latencyFract*speeds(si);
                for ti = 1:size(htemp,2)
                    if ~isempty(find(htemp(:,ti),1))
                        latTemp(ti) = find(htemp(:,ti),1);
                    else
                        latTemp(ti) = Inf;
                    end
                end
                
                trials{di,szi,si} = trialNumbers(~~(masks(:,di,szi,si).*sacMask.*eccMask));
                hvels{di,szi,si} = hvelocity(tWin(1):tWin(2),...
                    ~~(masks(:,di,szi,si).*sacMask.*eccMask));
                latency{di,szi,si} = latTemp(~~(masks(:,di,szi,si).*sacMask.*eccMask)) + tWin(1)-1;
                hpos{di,szi,si} = hposition(tWin(1):tWin(2),...
                    ~~(masks(:,di,szi,si).*sacMask.*eccMask));
                vpos{di,szi,si} = vposition(tWin(1):tWin(2),...
                    ~~(masks(:,di,szi,si).*sacMask.*eccMask));
                eccentricity{di,szi,si} = hpos{di,szi,si}.^2 + vpos{di,szi,si}.^2;
                
                m_hvel(:,di,szi,si) = mean(hvels{di,szi,si},2);
                std_hvel(:,di,szi,si) = std(hvels{di,szi,si},[],2);
                
                if simulationTest.On
                    for ti = 1:size(hvels{di,szi,si},2)
                        hvels{di,szi,si}(:,ti) = zeros(size(m_hvel(:,di,szi,si)));
                        lag{di,szi,si}(ti) = round(randn*simulationTest.stdlat(szi));
                        if lag{di,szi,si}(ti) > simulationTest.maxlat
                            lag{di,szi,si}(ti) = simulationTest.maxlat;
                        elseif lag{di,szi,si}(ti) < simulationTest.minlat
                            lag{di,szi,si}(ti) = simulationTest.minlat;
                        end
                        g = (simulationTest.stdg(szi)*randn+1);
                        if g < simulationTest.ming
                            g = simulationTest.ming;
                        end
                        hvels{di,szi,si}(-simulationTest.minlat+1:end-simulationTest.maxlat-1,ti) = ...
                            g*m_hvel((-simulationTest.minlat+1:end-simulationTest.maxlat-1)+lag{di,szi,si}(ti),di,szi,si);
                    end
                end
                
                bias(:,di,szi,si) = m_hvel(:,di,szi,si) - speeds(si);
                
                analysisV{di,szi,si} = mean(...
                    hvels{di,szi,si}(...
                    analysisWin(1)-tWin(1)+noSaccadeTime+1:analysisWin(2)-tWin(1)+noSaccadeTime+1,:),1);
                
                analysisV_pre{di,szi,si} = [];
                sacTime{di,szi,si} = [];
        end
        
    end    
end

%% Infer gain, latency and offset using Lee et. al.
if findLeeLatency
% Set up grid to search
lags = -40:2:50;
ggrid = linspace(0.2,1.9,41);
offgrid = linspace(-2,2,41);
[GGRID, OFFGRID] = meshgrid(ggrid,offgrid);
GGRID = GGRID(:);
OFFGRID = OFFGRID(:);
for di = 1:length(directions)
    for szi = 1:length(sizes)
        for si = 1:length(speeds)
            temp = hvels{di,szi,si};
            SSE = nan([1,size(temp,2),length(GGRID),length(lags)]);
            for lagi = 1:length(lags)   
                SSE(1,:,:,lagi) = Lee_algorithm(m_hvel(startStop(1):startStop(2),di,szi,si),...
                    temp((startStop(1):startStop(2))+lags(lagi),:),...
                    GGRID,OFFGRID);
            end
            for ti = 1:size(temp,2)
                [i1,i2] = find(squeeze(SSE(1,ti,:,:)) == min(min(SSE(1,ti,:,:))));
%                 [i1,i2,i3,i4] = ind2sub(size(SSE(:,ti,:,:)),minInd);
                if length(i1) > 1
                    disp([di szi si ti])
                    latency2{di,szi,si}(ti) = NaN;
                    gains{di,szi,si}(ti) = NaN;
                    offsets{di,szi,si}(ti) = NaN;
                else
                    latency2{di,szi,si}(ti) = lags(i2);
                    gains{di,szi,si}(ti) = GGRID(i1);
                    offsets{di,szi,si}(ti) = OFFGRID(i1);
                end
            end
            
        end
    end
end

for szi = 1:length(sizes)
    templat = [latency2{:,szi,:}];
    stdlat(szi) = nanstd(templat);
    mlat(szi) = nanmean(templat);
    tempg = [gains{:,szi,:}];
    stdg(szi) = nanstd(tempg);
    mg(szi) = nanmean(tempg);
    tempo = [offsets{:,szi,:}];
    stdo(szi) = nanstd(tempo);
    mo(szi) = nanmean(tempo);
end
end
%% Use inferred latency to align to motion onset
if findLeeLatency    
%     for di = 1:length(directions)
%         for szi = 1:length(sizes)
%             for si = 1:length(speeds)
%                 indx(di,szi,si) = find(abs(m_hvel(:,di,szi,si)) > 0.6,1);
%             end
%         end
%     end                
                
    for di = 1:length(directions)
        for szi = 1:length(sizes)
            for si = 1:length(speeds)
                indx(di,szi,si) = find(abs(m_hvel(:,di,szi,si)) > 0.6,1);
                hvels2{di,szi,si} = nan(size(hvels{di,szi,si}));
                for ti = 1:size(hvels{di,szi,si},2)
                    htemp = zeros(size(hvels{di,szi,si}(:,ti),1),1);
                    htemp(-min(lags)+1:end-max(lags)-1) = hvels{di,szi,si}((-min(lags)+1:end-max(lags)-1)+latency2{di,szi,si}(ti),ti);
                    if alignToEye
                        hvels2{di,szi,si}(1-(indx(di,szi,si)-max(indx(:))):end,ti) = htemp(1:end+(indx(di,szi,si)-max(indx(:))),:);
                    else
                        hvels2{di,szi,si}(:,ti) = htemp;
                    end
                end
                
                analysisV{di,szi,si} = mean(...
                    hvels2{di,szi,si}(...
                    analysisWin(1)-tWin(1)+noSaccadeTime+1:analysisWin(2)-tWin(1)+noSaccadeTime+1,:),1);
                
                
                if outlierReject
                    analysisV{di,szi,si} = analysisV{di,szi,si}(abs(analysisV{di,szi,si}-nanmean(analysisV{di,szi,si})) < 2.5*nanstd(analysisV{di,szi,si}));
                    
                    disp([di szi si sum(abs(analysisV{di,szi,si}-nanmean(analysisV{di,szi,si})) > 2.5*nanstd(analysisV{di,szi,si}))])
                end
            end
        end
    end
end

%% Fit gainSDN model across sizes
if isnan(Ncv)
    for di = 1:length(directions)
        for szi = 1:length(sizes)
            for si = 1:length(speeds)
                tempM(di,szi,si) = nanmean(analysisV{di,szi,si});
                tempV(di,szi,si) = nanvar(analysisV{di,szi,si});
            end
        end
        if isnan(gainNoiseFitTrainIndices)
            gainNoiseFitTrainIndices = 1:length(speeds);
        end
        % Fit gain noise model
        [wSDN(di,1),sigG(di,1)] = fit_gainSDN(speeds(gainNoiseFitTrainIndices),permute(tempM(di,:,gainNoiseFitTrainIndices),[1,3,2]),...
            permute(tempV(di,:,gainNoiseFitTrainIndices),[1,3,2]),0.1,0.5,OPTIONS);
        
        % Fit fixed w_s model
        wSDN_noSigG(di) = fit_gainSDN_constrainedsigG(speeds(gainNoiseFitTrainIndices),permute(tempM(di,:,gainNoiseFitTrainIndices),[1,3,2]),...
            permute(tempV(di,:,gainNoiseFitTrainIndices),[1,3,2]),0.1,0,OPTIONS);
        
        % Fit gain noise model with gains as free parameter
        [wSDN2(di),sigG2(di), Gfit(di,:), baselineVar(di)] = fit_gainSDNlinear(cosd(directions(di))*speeds(gainNoiseFitTrainIndices),permute(tempM(di,:,gainNoiseFitTrainIndices),[1,3,2]),...
            permute(tempV(di,:,gainNoiseFitTrainIndices),[1,3,2]),0.1,0.5,[0.1,0.2,0.4],0,OPTIONS);
        
        % Fit gain noise model with motor noise
        [wSDNwm(di),sigGwm(di),wm(di)] = fit_gainSDNwm(speeds(gainNoiseFitTrainIndices),permute(tempM(di,:,gainNoiseFitTrainIndices),[1,3,2]),...
            permute(tempV(di,:,gainNoiseFitTrainIndices),[1,3,2]),0.1,0.5,0,OPTIONS);
        
        if isnan(gainNoiseFitTestIndices)
            gainNoiseFitTestIndices = 1:length(speeds);
        end
        
        % Performance of gain noise model
        gainSDN_RMSE(di,1) = sqrt( mean( reshape( permute(tempV(di,:,gainNoiseFitTestIndices),[1,3,2]) - ...
            gainSDN(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            repmat(speeds(gainNoiseFitTestIndices),[1,1,size(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),3)]),...
            wSDN(di),sigG(di)), [numel(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2])),1]).^2 ) );
        
        rval = corrcoef(permute(tempV(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            gainSDN(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            repmat(speeds(gainNoiseFitTestIndices),[1,1,size(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),3)]),...
            wSDN(di),sigG(di)));
        gainNoiseR2(di) = rval(1,2).^2;
                
        % Performance of fixed w model
        gainSDN_RMSE_noSigG(di,1) = sqrt( mean( reshape( permute(tempV(di,:,gainNoiseFitTestIndices),[1,3,2]) - ...
            gainSDN(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            repmat(speeds(gainNoiseFitTestIndices),[1,1,size(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),3)]),...
            wSDN_noSigG(di),0), [numel(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2])),1]).^2 ) );
        
        rval = corrcoef(permute(tempV(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            gainSDN(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            repmat(speeds(gainNoiseFitTestIndices),[1,1,size(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),3)]),...
            wSDN_noSigG(di),0));
        fixedWR2(di) = rval(1,2).^2;
        
        % Performance of gain noise model with motor noise
        gainSDNwm_RMSE(di,1) = sqrt( mean( reshape( permute(tempV(di,:,gainNoiseFitTestIndices),[1,3,2]) - ...
            gainSDNwm(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            repmat(speeds(gainNoiseFitTestIndices),[1,1,size(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),3)]),...
            wSDNwm(di),sigGwm(di),wm(di)), [numel(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2])),1]).^2 ) );
        
        rval = corrcoef(permute(tempV(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            gainSDNwm(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),...
            repmat(speeds(gainNoiseFitTestIndices),[1,1,size(permute(tempM(di,:,gainNoiseFitTestIndices),[1,3,2]),3)]),...
            wSDNwm(di),sigGwm(di),wm(di)));
        gainNoiseWmR2(di) = rval(1,2).^2;
    end
    
else
    
    for cvi = 1:Ncv
        for di = 1:length(directions)
            for szi = 1:length(sizes)
                for si = 1:length(speeds)
                    fitInds = randsample(length(analysisV{di,szi,si}),ceil(length(analysisV{di,szi,si})/2),true);
                    tempVec = false(1,length(analysisV{di,szi,si}));
                    tempVec(fitInds) = true;
                    tempM(di,szi,si) = nanmean(analysisV{di,szi,si}(tempVec));
                    tempV(di,szi,si) = nanvar(analysisV{di,szi,si}(tempVec));
                    loM(di,szi,si) = nanmean(analysisV{di,szi,si}(~tempVec));
                    loV(di,szi,si) = nanvar(analysisV{di,szi,si}(~tempVec));
                end
                
                w_standard(di,szi,cvi) = fit_gainSDN_constrainedsigG(...
                    speeds,permute(tempM(di,szi,:),[1,3,2]),...
                    permute(tempV(di,szi,:),[1,3,2]),0.1,0,OPTIONS);
                
                standardSDN_Errs(di,szi,:,cvi) = loV(di,szi,:) - ...
                    w_standard(di,szi,cvi).^2.*loM(di,szi,:).^2;
            end
            [wSDN(di,cvi),sigG(di,cvi)] = fit_gainSDN(speeds,permute(tempM(di,:,:),[1,3,2]),...
                permute(tempV(di,:,:),[1,3,2]),0.1,0.5,OPTIONS);
            wSDN_noSigG(di,cvi) = fit_gainSDN_constrainedsigG(speeds,permute(tempM(di,:,:),[1,3,2]),...
                permute(tempV(di,:,:),[1,3,2]),0.1,0,OPTIONS);
            
            gainSDN_RMSE(di,cvi) = sqrt( mean( reshape( permute(loV(di,:,:),[1,3,2]) - ...
                gainSDN(permute(loM(di,:,:),[1,3,2]),...
                repmat(speeds,[1,1,size(permute(loM(di,:,:),[1,3,2]),3)]),...
                wSDN(di,cvi),sigG(di,cvi)), [numel(permute(loM(di,:,:),[1,3,2])),1]).^2 ) );
            
            
            gainSDN_RMSE_noSigG(di,cvi) = sqrt( mean( reshape( permute(loV(di,:,:),[1,3,2]) - ...
                gainSDN(permute(loM(di,:,:),[1,3,2]),...
                repmat(speeds,[1,1,size(permute(loM(di,:,:),[1,3,2]),3)]),...
                wSDN_noSigG(di,cvi),0), [numel(permute(loM(di,:,:),[1,3,2])),1]).^2 ) );
            
            
            standardSDN_RMSE(di,cvi) = sqrt( mean(...
                reshape(standardSDN_Errs(di,:,:,cvi).^2,[numel(standardSDN_Errs(di,:,:,cvi)),1]) ) );
        end
    end
    
    for di = 1:length(directions)
        mD = mean(gainSDN_RMSE_noSigG(di,:)-standardSDN_RMSE(di,:));
        seD = std(gainSDN_RMSE_noSigG(di,:)-standardSDN_RMSE(di,:));
        tstatStandard(di) = mD./seD;
        pvalStandard(di) = tcdf(-abs(tstatStandard(di)),Ncv/2-1);
        mD = mean(gainSDN_RMSE_noSigG(di,:)-gainSDN_RMSE(di,:));
        seD = std(gainSDN_RMSE_noSigG(di,:)-gainSDN_RMSE(di,:));
        tstatGainNoise(di) = mD./seD;
        pvalGainNoise(di) = tcdf(-abs(tstatGainNoise(di)),Ncv/2-1);
    end
    
    
    wSDNwm = nan(1,length(directions));
    sigGwm = nan(1,length(directions));
    wm = nan(1,length(directions));
    gainSDNwm_RMSE = nan(1,length(directions));
    gainNoiseWmR2 = nan(1,length(directions));
    fixedWR2 = nan(1,length(directions));
    standardSDN_R2 = nan(1,length(directions));
    gainNoiseR2 = nan(1,length(directions));
    Gfit = nan(length(directions),length(sizes));
    wSDN2 = nan(1,length(directions));
    sigG2 = nan(1,length(directions));
    baselineVar = nan(1,length(directions));
end

%% Fit standard noise model to the data
if isnan(Ncv)
    for di = 1:length(directions)
        for szi = 1:length(sizes)
            for si = 1:length(speeds)
                tempM(di,szi,si) = nanmean(analysisV{di,szi,si});
                tempV(di,szi,si) = nanvar(analysisV{di,szi,si});
            end
            w_standard(di,szi) = fit_gainSDN_constrainedsigG(...
                speeds(gainNoiseFitTrainIndices),permute(tempM(di,szi,gainNoiseFitTrainIndices),[1,3,2]),...
                permute(tempV(di,szi,gainNoiseFitTrainIndices),[1,3,2]),0.1,0,OPTIONS);
            
            standardSDN_Errs(di,szi,:) = tempV(di,szi,gainNoiseFitTestIndices) - ...
                w_standard(di,szi).^2.*tempM(di,szi,gainNoiseFitTestIndices).^2;
            
        end
        standardSDN_RMSE(di,1) = sqrt( mean(...
            reshape(standardSDN_Errs(di,:,:).^2,[numel(standardSDN_Errs(di,:,:)),1]) ) );
        
        rval = corrcoef(squeeze(tempV(di,:,gainNoiseFitTestIndices)),...
            squeeze(w_standard(di,:).^2.*tempM(di,:,gainNoiseFitTestIndices).^2));
        standardSDN_R2(di) = rval(1,2).^2;
    end
end

%% Perform regression for data w/in a size

for szi = 1:length(sizes)
    eyeVs{szi} = [];
    targVs{szi} = [];
    
    for di = 1:length(directions)
        eyeVs2 = [];
        targVs2 = [];
        for si = 1:length(speeds)
            targVs{szi} = [targVs{szi};...
                cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}(:)))];
            eyeVs{szi} = [eyeVs{szi}; analysisV{di,szi,si}(:)];
            
            targVs2 = [targVs2; ...
                cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}(:)))];
            eyeVs2 = [eyeVs2; analysisV{di,szi,si}(:)];
            
            condMEAN(di,si,szi) = mean(analysisV{di,szi,si}(:));
            condVAR(di,si,szi) = var(analysisV{di,szi,si}(:),1);
            condV(di,si,szi) = cosd(directions(di))*speeds(si);
            condN(di,si,szi) = numel(analysisV{di,szi,si});
            
        end
        
        [Gains(:,szi,di) GainsINT, Rgains, STATSgains] = regress(eyeVs2,...
            [targVs2 ones(size(targVs2))]);
    end
    
    [B{szi},BINT{szi},R{szi},~,STATS{szi}] = regress(eyeVs{szi},...
        [targVs{szi} ones(size(targVs{szi}))]);
    
    rmse(szi) = sqrt(nanmean(R{szi}.^2));
end

%% Generate table of parameters and results
monkeyInfo = {d.sname(1);' ';' ';' ';' ';' ';' '};
models = {'$w_s$ fixed';'$w_s$ flexible';'';'';'';'Gain noise';'$w_m$ included'};
sizeinfo = {'';'';'2 deg';'6 deg';'20 deg';'';''};
for di = 1:length(directions)
    if directions(di) == 0
        directionInfo = {'R';' ';' ';' ';' ';' ';' '};
    elseif directions(di) == 180
        directionInfo = {'L';' ';' ';' ';' ';' ';' '};
    else
        directionInfo = {num2str(directions(di));' ';' ';' ';' ';' ';' ';' '};
    end
    w_s = {wSDN_noSigG(di); ' '; ...
        w_standard(di,1)'; w_standard(di,2)'; w_standard(di,3)';...
        wSDN(di);...
        wSDNwm(di)};
    sigmaG = {0; ' ';...
        '-'; '-'; '-'; ...
        sigG(di); ...
        sigGwm(di)};
    w_ms = {0; '';...
        0; 0; 0;...
        0;...
        wm(di)};
    RMSEs = {gainSDN_RMSE_noSigG(di); standardSDN_RMSE(di);' ';' ';' ';...
        gainSDN_RMSE(di);...
        gainSDNwm_RMSE(di)};
    Rsquared = {fixedWR2(di); standardSDN_R2(di);' ';' ';' ';...
        gainNoiseR2(di);...
        gainNoiseWmR2(di)};
    
    resultsTable{di} = table(monkeyInfo,directionInfo,models,sizeinfo,w_s,sigmaG,w_ms,RMSEs,Rsquared);
    resultsTable{di}.Properties.VariableNames = {'Monkey','Direction','Model','Size','w_s','sigma','w_m','RMSE','R2'};
    if saveTable.On
        table2latex(vertcat(resultsTable{:}),[saveTable.directory '/' d.sname '_table' datestr(now(),'yyyymmdd')]);
    end
end

%% Plotting

if plotflg
    %% Single trial response
    figure('Name','Single trial eye response','Position',[197 98 1101 665])
    inds = 0;
    for di = 1:length(directions)
        for szi = 1:length(sizes)
            inds = inds + 1;
            subplotH(di,szi) = subplot(length(directions),length(sizes),inds);
        end
    end
    di = 1;
    if length(speeds) >= 4
        si = 4;
    else
        si = length(speeds);
    end
        for szi = 1:length(sizes)
            subplot(subplotH(1,szi))
            ctemp = colors(si,:);
            ctemp(ctemp > 1) = 1;
            sampInd = randsample(size(hvels{di,szi,si},2),1); %38;
            if ~isempty(hvels{di,szi,si})
                plot((tWin(1):tWin(2))-1,hvels{di,szi,si}(:,sampInd),'Color',ctemp)
            end
            hold on
            
            subplot(subplotH(2,szi))
            ctemp = colors(si,:);
            ctemp(ctemp > 1) = 1;
            if ~isempty(hvels{di,szi,si})
                plot((tWin(1):tWin(2))-1,cumsum(hvels{di,szi,si}(:,sampInd))/1000,'Color',ctemp)
            end
            hold on
        end
        
    
    for di = 1:length(directions)
        if di == 1
            ylimit = 30;
        else
            ylimit = 2;
        end
        for szi = 1:length(sizes)
            subplot(subplotH(di,szi))
            axis([tWin(1)-1,tWin(2),0 ylimit])
            
            if di == 1
                plotHorizontal(speeds(si));
            end
            
            
            if szi == 2 && di == 2
                xlabel('Time from motion onset (ms)')
            elseif szi == 1
                ylabel('Eye speed (deg/s)')
            end
            mymakeaxis(gca,...
                'xytitle',['Dir = ' num2str(directions(di)) ', Sz = ' num2str(sizes(szi))],...
                'xticks',[0 100 200])
        end
    end
    
    %% Eye traces
    figure('Name','Horizontal speeds','Position',[197 98 1101 665])
    inds = 0;
    for di = 1:length(directions)
        for szi = 1:length(sizes)
            inds = inds + 1;
            subplotH(di,szi) = subplot(length(directions),length(sizes),inds);
        end
    end
    
    switch plotType
        case {'default','Default',0}
            for di = 1:length(directions)
                for szi = 1:length(sizes)
                    subplot(subplotH(di,szi))
                    for si = 1:length(speeds)
                        ctemp = colors(si,:) + [0.2 0.2 0.2];
                        ctemp(ctemp > 1) = 1;
                        if ~isempty(hvels{di,szi,si})
                            plot((tWin(1):tWin(2))-1,hvels{di,szi,si}(:,randsample(size(hvels{di,szi,si},2),10)),'Color',ctemp)
                        end
                        hold on
                    end
                end
            end
            
            for di = 1:length(directions)
                for szi = 1:length(sizes)
                    subplot(subplotH(di,szi))
                    for si = 1:length(speeds)
                        plot((tWin(1):tWin(2))-1,m_hvel(:,di,szi,si),...
                            'LineWidth',2,'Color',colors(si,:))
                    end
                    axis([tWin(1),tWin(2),-30 30])
                    
                    plotVertical(analysisWin(1)+noSaccadeTime);
                    plotVertical(analysisWin(2)+noSaccadeTime);
                    
                    
                    if szi == 2 && di == 2
                        xlabel('Time from motion onset (ms)')
                    elseif szi == 1
                        ylabel('Eye speed (deg/s)')
                    end
                    mymakeaxis(gca,...
                        'xytitle',['Dir = ' num2str(directions(di)) ', Sz = ' num2str(sizes(szi))],...
                        'xticks',[0 100 200])
                end
            end
            
        case {'FigureComposer','fyp'}
            
            for di = 1:length(directions)
                for szi = 1:length(sizes)
                    subplot(subplotH(di,szi))
                    for si = 1:length(speeds)
                        plot((tWin(1):tWin(2)),m_hvel(:,di,szi,si),...
                            'LineWidth',2,'Color',colors(si,:))
                        hold on
                    end
                    axis([tWin(1)-1,tWin(2),-30 30])
                    
                    plotVertical(analysisWin(1)+noSaccadeTime);
                    plotVertical(analysisWin(2)+noSaccadeTime);
                    
                    
                    if szi == 2 && di == 2
                        xlabel('Time from motion onset (ms)')
                    elseif szi == 1
                        ylabel('Eye speed (deg/s)')
                    end
                end
            end
            
        otherwise
            error(['Plot type ' plotType ' not recognized!'])
    end
    
    %% Plot mean eye velocity in analysis window sorted by target velocity
    switch plotType
        case {'default','Default',0}
            figure('Name','mean speeds in window','Position',[197 98 1101 665])
            inds = 0;
            for szi = 1:length(sizes)
                inds = inds + 1;
                for di = 1:length(directions)
                    subplotH2(di,szi) = subplot(1,length(sizes),inds);
                end
            end
            for szi = 1:length(sizes)
                for di = 1:length(directions)
                    subplot(subplotH2(di,szi))
                    for si = 1:length(speeds)
                        if ~isempty(analysisV{di,szi,si})
                            plot(cosd(directions(di))*speeds(si),randsample(analysisV{di,szi,si},15),'o',...
                                'Color',colors(si,:),'MarkerFaceColor',colors(si,:))
                        end
                        hold on
                    end
                end
                vs = -30:30;
                plot(vs,B{szi}(1)*vs + B{szi}(2),'k')
            end
            
            for di = 1:length(directions)
                for szi = 1:length(sizes)
                    subplot(subplotH2(di,szi))
                    axis([-30 30 -30 30])
                    plotUnity;
                    axis square
                    if szi == 2 && di == 2
                        xlabel('Target speed (deg/s)')
                    elseif szi == 1
                        ylabel('Eye speed (deg/s)')
                    end
                    mymakeaxis(gca,...
                        'xytitle',['Sz = ' num2str(sizes(szi))])
                end
            end
    
        case {'FigureComposer','fyp'}
            figure('Name','mean speeds in window','Position',[197 98 500 430])
            
            for szi = 1:length(sizes)
                szind = length(sizes)-szi+1;
                for di = 1:length(directions)
                    for si = 1:length(speeds)
                        if ~isempty(analysisV{di,szind,si})
                            if size(analysisV{di,szind,si}) < 15
                                plot(cosd(directions(di))*speeds(si),analysisV{di,szind,si},'o',...
                                    'Color',szcolors(szind,:),'MarkerFaceColor',szcolors(szind,:))
                            else
                                plot(cosd(directions(di))*speeds(si),randsample(analysisV{di,szind,si},15),'o',...
                                    'Color',szcolors(szind,:),'MarkerFaceColor',szcolors(szind,:))
                            end
                        end
                        hold on
                    end
                end
                vs = -30:30;
                plot(vs,B{szind}(1)*vs + B{szind}(2),'Color',szcolors(szind,:))
            end
            plotUnity;
            axis square
            xlabel('Target speed (deg/s)')
            ylabel('Eye speed (deg/s)')
            
        otherwise
            error(['Plot type ' plotType ' not recognized!'])
    end
    
    figure('Name','Slopes +/- 95% CI')
    % Determine significance of slope differences
    z = nan(length(sizes));
    pDiff = nan(length(sizes));
    for ii = 1:length(sizes)
        for jj = 1:length(sizes)
            if jj > ii
                se1(ii,jj) = (BINT{ii}(1,2) - B{ii}(1))/1.96;
                se2(ii,jj) = (BINT{jj}(1,2) - B{jj}(1))/1.96;
                z(ii,jj) = ( B{ii}(1)-B{jj}(1) )/sqrt( se1(ii,jj).^2+se2(ii,jj).^2 );
                pDiff(ii,jj) = 2*normcdf(-abs(z(ii,jj)),0,1);
            end
        end
    end           
    
    for szi = 1:length(sizes)
        errorbar(szi,B{szi}(1),BINT{szi}(1,1)-B{szi}(1),BINT{szi}(1,2)-B{szi}(1),...
            'k.','MarkerFaceColor',[0 0 0])
        xticklabels{szi} = num2str(sizes(szi));
        hold on
    end
    
    allB = horzcat(B{:});
    if pDiff(1,2) < 0.01
        plot(1.5,max(allB(1,:))*1.1,'k*')
        text(1.5,max(allB(1,:))*1.2,['p = ' num2str(pDiff(1,2))])
    end
    if pDiff(2,3) < 0.01
        plot(2.5,max(allB(1,:))*1.1,'k*')
        text(2.5,max(allB(1,:))*1.2,['p = ' num2str(pDiff(2,3))])
    end
    
    
    xlabel('Target size')
    ylabel('Slope (+/- 95% CI)')
    axis square
    mymakeaxis(gca,'xticks',[1,2,3],'xticklabels',xticklabels)
    
    %% Variance as a function of eye speed
    inds = 0;
       
        figure('Name','Horizontal speeds','Position',[197 98 414 665])
        for di = 1:2
            inds = inds + 1;
            subplotH_var(di,1) = subplot(length(directions),1,inds);
        end
        
        switch plotType
            case {'default','Default',0}
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    for szi = 1:length(sizes)
                        ys = [];
                        xs = [];
                        for si = 1:length(speeds)
                            tempM(di,szi,si) = nanmean(analysisV{di,szi,si});
                            tempV(di,szi,si) = nanvar(analysisV{di,szi,si});
                            plot(nanmean(analysisV{di,szi,si}),...
                                nanvar(analysisV{di,szi,si}),'s','MarkerSize',10,...
                                'Color',colors(si,:),'MarkerFaceColor',szcolors(szi,:))
                            hold on
                            
                            ys = [ys analysisV{di,szi,si}];
                            xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
                        end
                        
                        ys = ys(~isnan(ys));
                        xs = xs(~isnan(ys));
                        betas(:,szi) = regress(ys(:),[xs(:)]);
                        x = linspace(0,cosd(directions(di))*max(speeds),100);
                        plot(betas(1,szi)*x,...
                            gainSDN(betas(1,szi)*x,x,wSDN(di),sigG(di)),...
                            'Color',szcolors(szi,:));
                        hold on
                    end
                    ax(di,:) = axis;
                    hold on
                end
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    if cosd(directions(di)) < 0
                        axis([cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 0 max(ax(:,4))])
                        xticks = [-10 0];
                    else
                        axis([0 cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 max(ax(:,4))])
                        xticks = [0 10];
                    end
                    xlabel('Eye speed (deg/s)')
                    ylabel('VAR (deg/s)^2')
                    ylim([0 3.5])
                    axis square
                    mymakeaxis(gca,'xytitle',['Dir = ' num2str(directions(di))],...
                        'xticks',xticks)%,'yticks',[0 2])
                end
                
            case {'FigureComposer','fyp'}
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    for szi = 1:length(sizes)
                        ys = [];
                        xs = [];
                        for si = 1:length(speeds)
                            tempM(di,szi,si) = nanmean(analysisV{di,szi,si});
                            tempV(di,szi,si) = nanvar(analysisV{di,szi,si});
                            ys = [ys analysisV{di,szi,si}];
                            xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
                        end
                        plot(permute(tempM(di,szi,:),[3,1,2]),...
                            permute(tempV(di,szi,:),[3,1,2]),'o',...
                            'Color',szcolors(szi,:),'MarkerFaceColor',szcolors(szi,:))
                        hold on
                        
                        ys = [ys analysisV{di,szi,si}];
                        xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
                        
                        ys = ys(~isnan(ys));
                        xs = xs(~isnan(ys));
                        betas(:,szi) = regress(ys(:),[xs(:)]);
                        x = linspace(0,cosd(directions(di))*max(speeds),100);
                        plot(betas(1,szi)*x,...
                            gainSDN(betas(1,szi)*x,x,mean(wSDN(di,:)),mean(sigG(di,:))),...
                            'Color',szcolors(szi,:));
                        hold on
                    end
                    ax(di,:) = axis;
                    hold on
                end
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    if cosd(directions(di)) < 0
                        axis([cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 0 max(ax(:,4))])
                        xticks = [-10 0];
                    else
                        axis([0 cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 max(ax(:,4))])
                        xticks = [0 10];
                    end
                    xlabel('Eye speed (deg/s)')
                    ylabel('VAR (deg/s)^2')
                    ylim([0 3.5])
                    axis square
                end
                
            otherwise
                error(['Plot type ' plotType ' not recognized!'])                
        end

    %% Flexible SDN model fit
    inds = 0;
        figure('Name','Flexible SDN','Position',[197 98 414 665])
        for di = 1:2
            inds = inds + 1;
            subplotH_var(di,1) = subplot(length(directions),1,inds);
        end
        
        switch plotType
            case {'default','Default',0}
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    for szi = 1:length(sizes)
                        ys = [];
                        xs = [];
                        for si = 1:length(speeds)
                            tempM(di,szi,si) = nanmean(analysisV{di,szi,si});
                            tempV(di,szi,si) = nanvar(analysisV{di,szi,si});
                            plot(nanmean(analysisV{di,szi,si}),...
                                nanvar(analysisV{di,szi,si}),'s','MarkerSize',10,...
                                'Color',colors(si,:),'MarkerFaceColor',szcolors(szi,:))
                            hold on
                            
                            ys = [ys analysisV{di,szi,si}];
                            xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
                        end
                        
                        ys = ys(~isnan(ys));
                        xs = xs(~isnan(ys));
                        betas(:,szi) = regress(ys(:),[xs(:)]);
                        x = linspace(0,cosd(directions(di))*max(speeds),100);
                        plot(betas(1,szi)*x,...
                            w_standard(di,szi).^2*(betas(1,szi)*x).^2,...
                            'Color',szcolors(szi,:));
                        hold on
                    end
                    ax(di,:) = axis;
                    hold on
                end
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    if cosd(directions(di)) < 0
                        axis([cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 0 max(ax(:,4))])
                        xticks = [-10 0];
                    else
                        axis([0 cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 max(ax(:,4))])
                        xticks = [0 10];
                    end
                    xlabel('Eye speed (deg/s)')
                    ylabel('VAR (deg/s)^2')
                    ylim([0 3.5])
                    axis square
                    mymakeaxis(gca,'xytitle',['Dir = ' num2str(directions(di))],...
                        'xticks',xticks)%,'yticks',[0 2])
                end
                
            case {'FigureComposer','fyp'}
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    for szi = 1:length(sizes)
                        ys = [];
                        xs = [];
                        for si = 1:length(speeds)
                            tempM(di,szi,si) = nanmean(analysisV{di,szi,si});
                            tempV(di,szi,si) = nanvar(analysisV{di,szi,si});
                            ys = [ys analysisV{di,szi,si}];
                            xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
                        end
                        plot(permute(tempM(di,szi,:),[3,1,2]),...
                            permute(tempV(di,szi,:),[3,1,2]),'o',...
                            'Color',szcolors(szi,:),'MarkerFaceColor',szcolors(szi,:))
                        hold on
                        
                        ys = [ys analysisV{di,szi,si}];
                        xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
                        
                        ys = ys(~isnan(ys));
                        xs = xs(~isnan(ys));
                        betas(:,szi) = regress(ys(:),[xs(:)]);
                        x = linspace(0,cosd(directions(di))*max(speeds),100);
                        plot(betas(1,szi)*x,...
                            mean(w_standard(di,szi,:)).^2*(betas(1,szi)*x).^2,...
                            'Color',szcolors(szi,:));
                        hold on
                    end
                    ax(di,:) = axis;
                    hold on
                end
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    if cosd(directions(di)) < 0
                        axis([cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 0 max(ax(:,4))])
                        xticks = [-10 0];
                    else
                        axis([0 cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 max(ax(:,4))])
                        xticks = [0 10];
                    end
                    xlabel('Eye speed (deg/s)')
                    ylabel('VAR (deg/s)^2')
                    ylim([0 3.5])
                    axis square
                end
                
            otherwise
                error(['Plot type ' plotType ' not recognized!'])                
        end
        
    %% Model comparison
    figure('Name','Model comparision','Position',[1000 1024 867 305])
    barProperties.FaceColor = [0 0 0];
    barProperties.EdgeColor = 'none';
    barProperties.ShowBaseLine = 'off';
    barProperties.BarWidth = 0.4;
    
    for di = 1:length(directions)
        subplot(1,2,di)
        if isnan(Ncv)
            hbar = mybargraph([1,2,3,4],[mean(standardSDN_RMSE(di,:)) mean(gainSDN_RMSE(di,:))...
                mean(gainSDNwm_RMSE(di,:)) mean(gainSDN_RMSE_noSigG(di,:))],...
                'barProperties',barProperties);
        else
            hbar = mybargraph([1,2,3],[mean(standardSDN_RMSE(di,:)) mean(gainSDN_RMSE(di,:)) mean(gainSDN_RMSE_noSigG(di,:))],...
                'barProperties',barProperties);
        end
        hold on
        if ~isnan(Ncv)
            plot(1,standardSDN_RMSE(di,:),'ro')
            plot(2,gainSDN_RMSE(di,:),'ro')
            plot(3,gainSDN_RMSE_noSigG(di,:),'ro')
        end
        xlabel([num2str(directions(di)) ' deg'])
        ylabel('RMSE')
        axTemp(di,:) = axis;
    end
    
    for di = 1:length(directions)
        subplot(1,2,di)
        ylim([0,max(axTemp(:,4))]);
        sigStr = num2str(sigG(di));
        if isnan(Ncv)
            barlabels = {['w flex'],['\sigma_G = ' sigStr(1:5)],...
                ['w_m'],['\sigma_G = 0']};
            
            mymakeaxis(gca,'xticks',[1,2,3,4],'xticklabels',barlabels)
        else
            barlabels = {['w flex'],['\sigma_G = ' sigStr(1:5)],['\sigma_G = 0']};
            
            mymakeaxis(gca,'xticks',[1,2],'xticklabels',barlabels)
        end
    end
    

    %% Latency analysis
    figure('Name','Latency Distirbutions','Position',[195 315 995 368])
    
    edges = linspace(0,250,25);
    for szi = 1:length(sizes)
        subplot(3,2,2*szi-1)
        Nlats = histcounts([latency{:,szi,:}],edges);
        
        barProperties.FaceColor = szcolors(szi,:);
        barProperties.EdgeColor = 'none';
        barProperties.FaceAlpha = 1;
        barProperties.ShowBaseLine = 'off';
        barProperties.BarWidth = 0.9;
        
        hold on
        
        mybargraph(edges(1:end-1)+(edges(2)-edges(1))/2,Nlats/sum(Nlats),...
            'barProperties',barProperties);
        plotVertical(analysisWin(1));
        plotVertical(analysisWin(2));
        if szi == 2
           ylabel(['Cummulative probability > ' num2str(latencyFract) '*speed'])
        end
    end
    xlabel('Time from motion onset (ms)')
    
    for szi = 1:length(sizes)
        subplot(3,2,2*szi)
        hold on
        for di = 1:length(directions)
            for si = 1:length(speeds)
                inds = find(latency{di,szi,si} > analysisWin(2));
                if ~isempty(inds)
                    plot((tWin(1):tWin(2))-1,hvels{di,szi,si}(:,inds)/speeds(si),...
                        'Color',colors(si,:))
                end
            end
        end
        axis([tWin(1)-1 tWin(2) -1.2 1.2])
        plotVertical(analysisWin(1));
        plotVertical(analysisWin(2));
        plotHorizontal(latencyFract);
        plotHorizontal(-latencyFract);
        if szi == 2
            ylabel('(Eye speed) / (target speed)')
        end
    end
    xlabel('Time from motion onset (ms)')
    
    %%
    if findLeeLatency
    figure('Name','Lee et. al. analysis results')
    for szi = 1:length(sizes)
        subplot(3,3,3*szi-2)
        histogram([latency2{:,szi,:}],lags);
        if szi == 2
            ylabel('N trials')
        elseif szi == 3
            xlabel('Latency (ms)')
        end
        
        subplot(3,3,3*szi-1)
        histogram([gains{:,szi,:}],ggrid);
        if szi == 3
            xlabel('Gain')
        end
        
        subplot(3,3,3*szi)
        histogram([offsets{:,szi,:}],offgrid);
        if szi == 3
            xlabel('Offset (deg/s)')
        end
    end
    end
    
    %% Gfit vs B
    figure('Name','Slope of data vs that inferred from variance and target speed','Position',[440 24 348 774])
    inds = 0;
    for di = 1:length(directions)
        inds = inds + 1;
        subplotH_var(di,1) = subplot(length(directions)+1,1,inds);
    end
    inds = inds + 1;
    subplotH_var(length(directions)+1,1) = subplot(length(directions)+1,1,inds);
    for di = 1:length(directions)
        subplot(subplotH_var(di,1))
        for szi = 1:length(sizes)
            ys = [];
            xs = [];
            for si = 1:length(speeds)
                tempM(di,szi,si) = nanmean(analysisV{di,szi,si});
                tempV(di,szi,si) = nanvar(analysisV{di,szi,si});
                ys = [ys analysisV{di,szi,si}];
                xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
            end
            plot(permute(tempM(di,szi,:),[3,1,2]),...
                permute(tempV(di,szi,:),[3,1,2]),'o',...
                'Color',szcolors(szi,:),'MarkerFaceColor',szcolors(szi,:))
            hold on
            
            ys = [ys analysisV{di,szi,si}];
            xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
            
            ys = ys(~isnan(ys));
            xs = xs(~isnan(ys));
            betas(:,szi) = regress(ys(:),[xs(:)]);
            x = linspace(0,cosd(directions(di))*max(speeds),100);
            plot(Gfit(di,szi)*x,...
                gainSDNlinear(Gfit(di,szi),Gfit(di,szi)*x,wSDN2(di),sigG2(di),baselineVar(di)),...
                'Color',szcolors(szi,:));
            hold on
        end
        ax(di,:) = axis;
        hold on
    end
    subplot(subplotH_var(end,1))
    btemp = [B{:}];
    for di = 1:length(directions)
        plot(btemp(1,:),Gfit(di,:),'o')
        hold on
    end
    plotUnity;
    axis square
    xlabel('Speed sensitivity')
    ylabel('Inferred from variance')
    
    %% Results of gain noise model with motor noise
    inds = 0;
       
        figure('Name','Horizontal speeds','Position',[197 98 414 665])
        for di = 1:2
            inds = inds + 1;
            subplotH_var(di,1) = subplot(length(directions),1,inds);
        end
        
        switch plotType
            case {'default','Default',0}
                
                
            case {'FigureComposer','fyp'}
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    for szi = 1:length(sizes)
                        ys = [];
                        xs = [];
                        for si = 1:length(speeds)
                            tempM(di,szi,si) = nanmean(analysisV{di,szi,si});
                            tempV(di,szi,si) = nanvar(analysisV{di,szi,si});
                            ys = [ys analysisV{di,szi,si}];
                            xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
                        end
                        plot(permute(tempM(di,szi,:),[3,1,2]),...
                            permute(tempV(di,szi,:),[3,1,2]),'o',...
                            'Color',szcolors(szi,:),'MarkerFaceColor',szcolors(szi,:))
                        hold on
                        
                        ys = [ys analysisV{di,szi,si}];
                        xs = [xs cosd(directions(di))*speeds(si)*ones(size(analysisV{di,szi,si}))];
                        
                        ys = ys(~isnan(ys));
                        xs = xs(~isnan(ys));
                        betas(:,szi) = regress(ys(:),[xs(:)]);
                        x = linspace(0,cosd(directions(di))*max(speeds),100);
                        plot(betas(1,szi)*x,...
                            gainSDNwm(betas(1,szi)*x,x,wSDNwm(di),sigGwm(di),wm(di)),...
                            'Color',szcolors(szi,:));
                        hold on
                    end
                    ax(di,:) = axis;
                    hold on
                end
                for di = 1:2
                    subplot(subplotH_var(di,1))
                    if cosd(directions(di)) < 0
                        axis([cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 0 max(ax(:,4))])
                        xticks = [-10 0];
                    else
                        axis([0 cosd(directions(di))*max(max(abs(ax(:,1:2)))) 0 max(ax(:,4))])
                        xticks = [0 10];
                    end
                    xlabel('Eye speed (deg/s)')
                    ylabel('VAR (deg/s)^2')
                    ylim([0 3.5])
                    axis square
                end
                
            otherwise
                error(['Plot type ' plotType ' not recognized!'])                
        end
end

%% Fuctions

%% Lee_algorithm
function SSE = Lee_algorithm(m,h,GGRID,OFFGRID)
    %%
    m = repmat(m,[1,size(h,2),length(GGRID)]);
    h = repmat(h,[1,1,size(GGRID,1)]);
    G = repmat(permute(GGRID,[2, 3, 1]),[size(h,1), size(h,2), 1]);
    O = repmat(permute(OFFGRID,[2, 3, 1]),[size(h,1), size(h,2), 1]);
    
    E = m - (G.*h+O);
    SSE = sum( E.^2 , 1);

%% fit_gainSDN
function [w, sigG] = fit_gainSDN(speeds,VeM,VeVAR,w0,sigG0,OPTIONS)
    %%
    minimizer = @(p)( sum( sum( (VeVAR - gainSDN(VeM,repmat(speeds(:)',[size(VeM,1),1,size(VeM,3)]),p(1),p(2))).^2 ) ) );
    p = fmincon(minimizer,[w0 sigG0],[],[],[],[],[0,0],[Inf,Inf],[],OPTIONS);
    w = p(1);
    sigG = p(2);
    
%% fit_gainSDNwm
function [w, sigG, wm] = fit_gainSDNwm(speeds,VeM,VeVAR,w0,sigG0,wm0,OPTIONS)
    %%
    minimizer = @(p)( sum( sum( (VeVAR - gainSDNwm(VeM,repmat(speeds(:)',[size(VeM,1),1,size(VeM,3)]),p(1),p(2),p(3))).^2 ) ) );
    p = fmincon(minimizer,[w0 sigG0 wm0],[],[],[],[],[0,0,0],[Inf,Inf,Inf],[],OPTIONS);
    w = p(1);
    sigG = p(2);
    wm = p(3);
    
%% fit_gainSDN
function [w, sigG, G, baselineVar] = fit_gainSDNlinear(speeds,VeM,VeVAR,w0,sigG0,G0,baselineVar0,OPTIONS)
    %%
%     minimizer = @(p)( sum( sum( (VeVAR - gainSDNlinear(p(3:5),repmat(speeds(:)',[size(VeM,1),1,size(VeM,3)]),p(1),p(2),0)).^2 ) ) );
    minimizer = @(p)( sum( sum( (VeVAR - gainSDNlinear(p(3:5),VeM,p(1),p(2),0)).^2 ) ) );
    p = fmincon(minimizer,[w0 sigG0 G0],[],[],[],[],[0,0,0,0,0],[Inf,Inf,1,1,1],[],OPTIONS);
    w = p(1);
    sigG = p(2);
    G = p(3:5);
    baselineVar = 0;
    
%% fit_gainSDNlinearComplex
function [w, sigG, G, baselineVar] = fit_gainSDNlinearComplex(speeds,VeM,VeVAR,w0,sigG0,G0,baselineVar0,OPTIONS)
    %%
    minimizer = @(p)( sum( sum( (VeM - repmat(permute(p(3:5),[1,3,2]),[size(VeM,1),size(VeM,2),1]).*repmat(speeds(:)',[size(VeM,1),1,size(VeM,3)])).^2 ) ) + ...
        sum( sum( (VeVAR - gainSDNlinear(p(3:5),repmat(speeds(:)',[size(VeM,1),1,size(VeM,3)]),p(1),p(2),p(6))).^2 ) ) );
    p = fmincon(minimizer,[w0 sigG0 G0 baselineVar0],[],[],[],[],[0,0,0,0,0,0],[Inf,Inf,1.5,1.5,1.5,Inf],[],OPTIONS);
    w = p(1);
    sigG = p(2);
    G = p(3:5);
    baselineVar = p(6);
    
%% gainSDN
function v = gainSDN(ve,vs,w,sigG)
    %%
    v = w.^2.*ve.^2 + (sigG.^2 + sigG.^2.*w.^2).*vs.^2;
    
%% gainSDNlinear
function v = gainSDNlinear(G,ve,w,sigG,baselineVar)
    %%
%     v = baselineVar +  w.^2.*repmat(permute(G(:),[3,2,1]),[size(vs,1),size(vs,2),1]).^2.*vs.^2 + (sigG.^2 + sigG.^2.*w.^2).*vs.^2;
    v = baselineVar +  w.^2.*ve.^2 + (sigG.^2 + sigG.^2.*w.^2).*ve.^2./repmat(permute(G(:),[3,2,1]),[size(ve,1),size(ve,2),1]).^2;

%% gainSDNwm
function v = gainSDNwm(ve,vs,w,sigG,wm)
    %%
    v = (w.^2 + wm.^2).*ve.^2 + (sigG.^2 + sigG.^2.*w.^2).*vs.^2;
    
%% fit_gainSDN_constrainedsigG
function w = fit_gainSDN_constrainedsigG(speeds,VeM,VeVAR,w0,sigG,OPTIONS)
    %%
    minimizer = @(p)( sum( sum( (VeVAR - gainSDN(VeM,repmat(speeds(:)',[size(VeM,1),1,size(VeM,3)]),p(1),sigG)).^2 ) ) );
    p = fmincon(minimizer,[w0,0.1],[],[],[],[],[0,0],[Inf,Inf],[],OPTIONS);
    w = p(1);  
    
    
    