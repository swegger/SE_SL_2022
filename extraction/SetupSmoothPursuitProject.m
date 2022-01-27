function d = SetupSmoothPursuitProject(sname,projectName,directory,varargin)
%% SetupSmoothPursuitProject
%
%   d = SetupSmoothPursuitProject(sname,projectName,directory)
%
%   Sets up a data structure for smooth pursuit. Returns data structure
%   with the following fields:
%       sname - subject name
%       projectName - name of the project
%       projectFnc - function designed to extract Maesto data
%       projectDir - project directory
%       meta - meta data about each data collection run
%   and saves the contents of the structure to 
%       directory/data/sname/sname_projectName.mat
%
%   If directory/data/sname/sname_projectName.mat already exists, the data
%   will be loaded into the variable d.
%
%   If unextracted data exists under the directory/data/sname, the new data
%   will be extracted from the corresponding Maesto files and saved to
%       directory/data/sname/sname_projectName.mat
%
%
%%

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'sname')
addRequired(Parser,'projectName')
addRequired(Parser,'directory')

parse(Parser,sname,projectName,directory,varargin{:})

sname = Parser.Results.sname;
projectName = Parser.Results.projectName;
directory = Parser.Results.directory;

%% Check if project file for sname exists; if not create 
cd([directory '/data']);
test = exist(sname,'file');
if ~test
    mkdir(sname);
end

%% Check if project data file exists for sname; if not create
cd([directory '/data/' sname])
datafile = [directory '/data/' sname '/' sname '_' projectName '.mat'];
test = exist(datafile,'file');
if ~test
    d.sname = sname;
    d.projectName = projectName;
    d.projectFnc = eval(['@' projectName]);
    d.projectDir = directory;
    save(datafile,'-struct','d')
    
else
    d = load(datafile);
end

%% Determine if any data files exist that have not been processed
files = dir([d.projectDir '/data/' d.sname '/']);
datafileInds = find(vertcat(files.isdir));

if length(datafileInds) > 2
    for i = 3:length(datafileInds)
        irun = i-2;
        fileName = files(datafileInds(i)).name;
        
        if isfield(d,'meta')
            if length(d.meta) < irun
                dataDir = [d.projectDir '/data/' sname '/' fileName];
                d.meta(irun).datafile = dataDir;
                
                d = d.projectFnc(irun,d,fileName);
            elseif ~strcmp(d.meta(irun).datafile(end-7:end),fileName)
                dataDir = [d.projectDir '/data/' sname '/' fileName];
                d.meta(irun).datafile = dataDir;
                
                d = d.projectFnc(irun,d,fileName);
            end
        else
            dataDir = [d.projectDir '/data/' sname '/' fileName];
            d.meta(irun).datafile = dataDir;
            d = d.projectFnc(irun,d,fileName);
        end
    end
else
    disp(['No data files for subject ' sname ', project ' projectName])
end

%% Save output
save(datafile,'-struct','d','-v7.3')


