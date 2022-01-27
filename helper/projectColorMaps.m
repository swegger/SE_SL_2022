function colors = projectColorMaps(type,varargin)
%% projectColorMaps
%
%   colors = projectColorMaps(type)
%
%   Returns a matrix of values corresponding to the RGB values (columns)
%   for 64 linearly sample positions in the color map.
%
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'type')
addParameter(Parser,'samples',NaN)
addParameter(Parser,'sampleN',64)
addParameter(Parser,'sampleDepth',64)

parse(Parser,type,varargin{:})

type = Parser.Results.type;
samples = Parser.Results.samples;
sampleN = Parser.Results.sampleN;
sampleDepth = Parser.Results.sampleDepth;

if isnan(samples)
    samples = round(linspace(1,64,sampleN));
end

%% Find colors
switch type
    case 'speeds'
        cmap = myrgb(64,[119,154,212]/255,[193,119,176]/255);
        ctemp = cmap(round(linspace(1,64,sampleDepth)),:);
        colors = ctemp(samples,:);
        
        
    case 'popN'
        cmap = myrgb(64,[100,100,100]/255,[193,0,0]/255);
        ctemp = cmap(round(linspace(1,64,sampleDepth)),:);
        colors = ctemp(samples,:);
        
    otherwise
        error(['Sample type ' type ' not supported!'])

end