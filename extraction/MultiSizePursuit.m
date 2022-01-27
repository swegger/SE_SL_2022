function d = MultiSizePursuit(irun,d,date,varargin)
%% MultiSizePursuit
%
%   d = MultiSizePursuit(d,sname,date)
%
%   Analysis of trials collected using multiSize trial sets.
%
%%

%% Defaults
ConditionFile_default = 'None';

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'irun')
addRequired(Parser,'d')
addRequired(Parser,'date')
addParameter(Parser,'ConditionFile',ConditionFile_default)
addParameter(Parser,'plotflg',false)

parse(Parser,irun,d,date,varargin{:})

irun = Parser.Results.irun;
d = Parser.Results.d;
date = Parser.Results.date;
ConditionFile = Parser.Results.ConditionFile;
plotflg = Parser.Results.plotflg;

%% Meta data
% irun = find(strcmp(d.meta(:).datafile(end-7:end),date));
d.meta(irun).ConditionFile = ConditionFile;


%% Load condition file
switch ConditionFile
    case 'None'
        CondSwitch = false;
        load([d.projectDir '/analysis/extraction/multiSize_ConditionTable'])
        szs = unique(vertcat(condition(:).sz));
    otherwise
        load([d.projectDir '/analysis/extraction/' ConditionFile])
        szs = unique(vertcat(condition(:).sz));
end


%% Move to data directory
cd(d.meta(irun).datafile)
files = dir;
for i = 1:length(files)
    temp(i) = strcmp(files(i).name(1),'.');
end
nonsenseFiles = sum(temp);

%%
for i = 1:length(files)-nonsenseFiles
    file = readcxdata(files(i+nonsenseFiles).name);
    d.meta(irun).trialInfo(i) = file.trialInfo;
    
    if CondSwitch
        if str2double(date) < 20190521
            underInd = strfind(file.trialname,'_');
            conditionInd = str2double(file.trialname(underInd+1:end));
        else
            [startIndex,endIndex] = regexp(file.trialname,'t\d{3}');
            conditionInd = str2double(file.trialname(startIndex+1:endIndex));
        end
        
        % Assign target data
        d.target.x{irun}(i) = condition(conditionInd).x;
        d.target.y{irun}(i) = condition(conditionInd).y;
        d.target.s{irun}(i) = condition(conditionInd).s;
        d.target.d{irun}(i) = condition(conditionInd).d;
        d.target.sz{irun}(i) = condition(conditionInd).sz;
        
        targInd = find(szs == condition(conditionInd).sz)+1;
        
    else
        filename = file.trialname;
        % Assign target data
        [startIndex,endIndex] = regexp(filename,'x\d{3}');
        d.target.x{irun}(i) = str2double(filename(startIndex+1:endIndex))-100;
        
        [startIndex,endIndex] = regexp(filename,'y\d{3}');
        d.target.y{irun}(i) = str2double(filename(startIndex+1:endIndex))-100;
        
        [startIndex,endIndex] = regexp(filename,'s\d{3}');
        d.target.s{irun}(i) = str2double(filename(startIndex+1:endIndex));
        
        [startIndex,endIndex] = regexp(filename,'d\d{3}');
        d.target.d{irun}(i) = str2double(filename(startIndex+1:endIndex));
        
        [startIndex,endIndex] = regexp(filename,'z\d{3}');
        d.target.sz{irun}(i) = str2double(filename(startIndex+1:endIndex));
        
        targInd = find(szs == str2double(filename(startIndex+1:endIndex)))+1;
    end
    d.target.hpos{irun}(:,i) = file.targets.hpos(targInd,:);
    d.target.vpos{irun}(:,i) = file.targets.vpos(targInd,:);
    d.target.hvel{irun}(:,i) = file.targets.hvel(targInd,:);
    d.target.vvel{irun}(:,i) = file.targets.vvel(targInd,:);
    d.target.patvelH{irun}(:,i) = file.targets.patvelH(targInd,:);
    d.target.patvelV{irun}(:,i) = file.targets.patvelV(targInd,:);
    graylevel = dec2hex(file.tgtdefns(targInd).params.iRGBMean(1));
    d.target.graylevel{irun}(:,i) = sscanf(graylevel.','%2x');
    
    % Detect saccades
    start=file.trialInfo.segStart(2);
    sacs = saccadeDetect(file.data(3,:)*.09189,file.data(4,:)*.09189,...
        'accelerationThreshold',1.1,'windowSize',40);
    temph = file.data(3,:);
    temph(sacs) = NaN;
    mh = 0;%nanmean(temph(start-150:start-100));
    tempv = file.data(4,:);
    tempv(sacs) =  NaN;
    mv = 0;%nanmean(tempv(start-150:start-100));
    
    % Assign eye data
    d.eye.saccades{irun}(:,i) = sacs;
    d.eye.hpos{irun}(:,i) = (file.data(1,:) - mean(file.data(1,start-150:start)))*.025;
    d.eye.vpos{irun}(:,i) = (file.data(2,:) - mean(file.data(2,start-150:start)))*.025;
    d.eye.hvel{irun}(:,i) = (file.data(3,:) - mh)*.09189;
    d.eye.vvel{irun}(:,i) = (file.data(4,:) - mv)*.09189;
    
    d.time.startTime{irun}(i) = file.trialInfo.segStart(3)-file.trialInfo.segStart(2);
    d.time.endTime{irun}(i) = file.trialInfo.segStart(end);
end

%% Plotting
if plotflg
    inds = find(d.target.x{irun} == 0 & ...
        d.target.y{irun} == 0 & ...
        d.target.d{irun} == 0 & ...
        d.target.sz{irun} == 2);
    figure('Name','Target and eye position velocity')
    for i = 1:length(inds)
        subplot(2,2,1)
        t = 0:size(d.target.hpos{irun}(:,inds(i)),1)-1;
        t = t-d.time.startTime{irun}(inds(i));
        plot(t,d.target.hpos{irun}(:,inds(i))+d.target.x{irun}(inds(i)),'k');
        hold on
        plot(t,d.eye.hpos{irun}(:,inds(i)),'b');
        plotVertical(150);
        
        subplot(2,2,2)
        plot(t,d.target.vpos{irun}(:,inds(i))+d.target.y{irun}(inds(i)),'k');
        hold on
        plot(t,d.eye.vpos{irun}(:,inds(i)),'b');
        plotVertical(150);

        subplot(2,2,3)
        plot(t,d.target.hvel{irun}(:,inds(i)),'k');
        hold on
        hvel = d.eye.hvel{irun}(:,inds(i));
        hvel(d.eye.saccades{irun}(:,inds(i))) = NaN;
        plot(t,hvel,'b');
        plotVertical(150);
        
        subplot(2,2,4)
        vvel = d.eye.vvel{irun}(:,inds(i));
        vvel(d.eye.saccades{irun}(:,inds(i))) = NaN;        
        plot(t,d.target.vvel{irun}(:,inds(i)),'k');
        hold on
        plot(t,vvel,'b')
        plotVertical(150);
    end
end