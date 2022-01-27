function sacInds = saccadeDetect(h,v,varargin)
%% saccadeDetect
%
%   sacTimes = saccadeDetect(h,v)
%
%   Detects the indices of a saccadic eye movement by crossing of default
%   threshold of acceleration and interpolating across time.
%
%%

%% Defaults

%% Parse inputs
Parser = inputParser;

addRequired(Parser,'h')
addRequired(Parser,'v')
addParameter(Parser,'accelerationThreshold',2)
addParameter(Parser,'windowSize',10)

parse(Parser,h,v,varargin{:})

h = Parser.Results.h;               % array of horizontal eye velocities (row = trial, column = time)
v = Parser.Results.v;               % array of vertical eye velocities (row = trial, column = time)
accelerationThreshold = Parser.Results.accelerationThreshold;
windowSize = Parser.Results.windowSize;

%% Find saccdes
ha = [zeros(size(h,1),1),diff(h,1,2)];
va = [zeros(size(v,1),1),diff(v,1,2)];

inds = abs(ha) > accelerationThreshold | abs(va) > accelerationThreshold;
sinds = smooth(inds,windowSize);
sacInds = sinds > 0;

%% Find times of saccades

