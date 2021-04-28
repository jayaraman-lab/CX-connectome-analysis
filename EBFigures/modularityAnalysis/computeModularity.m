function [Q,Qshuffle,p,A] = computeModularity(varargin)

% COMPUTEMODULARITY computes the EB-wedge-specific modularity of inputs from 
% ring neurons to EPGs, or vice versa
%
% USAGE:
%       computeModularity('erType','ER4d','prePost','pre','weightType','weightRelative','simMeas','corr','nShuffles',1000);
%
% OPTIONAL INPUTS:
%
%  erType (default: 'ER2_a')
%            specifies ring neuron type to consider
% 
%  weightType (default: 'weightRelative')
%            specifies whether to use relative weights or synpatic counts for
%            computing similarity of inputs
% 
%  simMeas (default: 'corr')
%            specifies whether to use correlation or cosine similarity as a
%            measure for computing similarity
% 
%  nShuffles (default = 1000)
%            number of shuffles to perform of adjacency matrix for 
%            computing statistical significance. 
%
% OUTPUTS:
%
%         Q: modularity of measured (unshuffled) adjacency matrix
%   
%  Qshuffle: vector of modularity values for shuffled adjacency matrices
%
%         p: significance of true modularity (measured as fraction of shuffles 
%            whose modularity exceeds the measured value)
%
%         A: measured adjacency matrix

%----------------------------- parse inputs ------------------------------%

defaultERtype = 'ER2_a';
validERtypes  = {'ER1_a','ER1_b','ER2_a','ER2_b','ER2_c','ER2_d','ER4d','ER4m'};
checkERtype   = @(x) any(validatestring(x,validERtypes));

defaultWeightType = 'weightRelative';
validWeightTypes  = {'weightRelative','ROIweight'};
checkWeightType   = @(x) any(validatestring(x,validWeightTypes));

defaultSimType = 'corr';
validSimTypes  = {'corr','cos'};
checkSimType   = @(x) any(validatestring(x,validSimTypes));

defaultPrePostType = 'pre';
validPrePostTypes  = {'pre','post'};
checkPrePostType   = @(x) any(validatestring(x,validPrePostTypes));

defaultNshuffles = 1000;

p = inputParser();
p.addOptional('nShuffles',defaultNshuffles,@isnumeric)
p.addOptional('erType',defaultERtype,checkERtype)
p.addOptional('prePost',defaultPrePostType,checkPrePostType)
p.addOptional('weightType',defaultWeightType,checkWeightType)
p.addOptional('simMeas',defaultSimType,checkSimType)

p.parse(varargin{:});
erType     = p.Results.erType;
prePost    = p.Results.prePost;
weightType = p.Results.weightType;
simMeas    = p.Results.simMeas;
nShuffles  = p.Results.nShuffles;
%-------------------------------------------------------------------------%

%define simlarity measure
if strcmp(simMeas,'corr')
    similarity = @corrcoef;
elseif strcmp(simMeas,'cos')
    similarity = @cossim;
end

%select CSV file to load
if strcmp(prePost,'pre')
    csvfile = 'er2epg.csv';
elseif strcmp(prePost,'post')
    csvfile = 'epg2erEB.csv';
end

%sort EPGs
sortVec = {'EPG(PB08)_R1','EPG(PB08)_L8','EPG(PB08)_R2','EPG(PB08)_L7',...
    'EPG(PB08)_R3','EPG(PB08)_L6','EPG(PB08)_R4','EPG(PB08)_L5',...
    'EPG(PB08)_R5','EPG(PB08)_L4','EPG(PB08)_R6','EPG(PB08)_L3',...
    'EPG(PB08)_R7','EPG(PB08)_L2','EPG(PB08)_R8','EPG(PB08)_L1'};

%find columns to load
opts = detectImportOptions(csvfile);
varNames = opts.SelectedVariableNames;
colNames = {'from','to','name_from','name_to','databaseType_from','databaseType_to'};
cols = [];
for i=1:numel(colNames)
    cols(i) = find(strcmp(colNames{i}, varNames));
end
opts.SelectedVariableNames = cols;
nameMat = readmatrix(csvfile,opts,'outputType','string');

opts.SelectedVariableNames = find(strcmp(weightType, varNames));
weightMat = readmatrix(csvfile,opts,'outputType','double');


if strcmp(prePost,'pre')==1
    erCol  = find(strcmp('databaseType_from',colNames));
    epgCol = find(strcmp('name_to',colNames));
    erID   = find(strcmp('from',colNames));
    erSide = find(strcmp('name_from',colNames));
    epgID  = find(strcmp('to',colNames));
else
    erCol  = find(strcmp('databaseType_to',colNames));
    epgCol = find(strcmp('name_from',colNames));
    erID   = find(strcmp('to',colNames));
    erSide = find(strcmp('name_to',colNames));
    epgID  = find(strcmp('from',colNames));
end
indER = find(nameMat(:,erCol)==erType);

%select ring neurons
nameMat   = nameMat(  indER,:);
weightMat = weightMat(indER);

%identify unique ring neurons/EPGs
[erU, ierU]  = unique(nameMat(:,erID ));
[epgU,iepgU] = unique(nameMat(:,epgID));
[~,ierS]     = sort(nameMat(ierU,erSide));

%build connectivity matrix
connMat = zeros(numel(epgU),numel(erU));
for i=1:numel(epgU)
    for j=1:numel(erU)
        ind = find(nameMat(:,epgID)==epgU(i) & nameMat(:,erID)==erU(ierS(j)));
        if numel(ind)>0
            connMat(i,j) = weightMat(ind);
        elseif numel(ind)>1
            error('multiple entries')
        end
    end
end


%sort based on wedges
isort = [];
Anull = zeros(numel(epgU));
m = 1;
nEPGperWedge = [];
for i=1:numel(sortVec)
    ii = find(strcmp(nameMat(iepgU,epgCol),sortVec{i}));
    isort = [isort;ii];
    inds = m:(m+numel(ii)-1);
    Anull(inds,inds) = 1;
    m = m + numel(ii);
    nEPGperWedge(i) = numel(inds);
end

%compute adjacency matrix
connMat = connMat(isort,:);
A = similarity(connMat');
A = normalizeSim(A,simMeas);
Q = modularity(A,Anull);
Qnull = modularity(Anull,Anull);

%compute shuffles
Qshuffle = zeros(1,nShuffles);
for i=1:nShuffles
    %shufle er inputs to each epg:
    cShuffle = Shuffle(connMat,1);
    Ashuffle = similarity(cShuffle');
    Ashuffle = (Ashuffle+1)./2;
    Qshuffle(i) = modularity(Ashuffle,Anull);
end
p = numel(find(Qshuffle>Q))./nShuffles;

end

%auxiliary functions
function s = cossim(connMat)
num = connMat'*connMat;
k   = sqrt(sum(connMat.^2,2));
den = k*k';
s = num./den;
end

function A = normalizeSim(A,simMeas)
if strcmp(simMeas,'corr')
    A = (A+1)./2;
end
end

function Q = modularity(A,Anull)
twom = nansum(A(:));
k    = nansum(A,2);
Qmat = (A - (k*k')./twom).*Anull;
Q    = nansum(nansum(Qmat))./twom;
end


