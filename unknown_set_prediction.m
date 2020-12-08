clear all;
%% Load data
load('b.mat');
unknown_set = readtable("Unknown_set_Bacillus.xlsx");
unknown_set = table2struct(unknown_set);

%% features
tAI = importdata('unknown_tAI.txt');
tAI = tAI.data(:,2);
AGGFrequency = (extractfield(unknown_set,'AGGFrequency'))';
mRNAfold = mRNAfoldAVG(unknown_set)';
CCGATATfreq = promoter_freq(unknown_set,'CCGATAT');

% AA frequencies
ORFs = (extractfield(unknown_set,'ORF'))';
AAseq = cellfun(@nt2aa, ORFs, 'UniformOutput', false);
AAseq_length = cellfun(@length, AAseq);
allAAcount = cellfun(@aacount, AAseq);

Qfreq = (extractfield(allAAcount,'Q'))'./AAseq_length;
Hfreq = (extractfield(allAAcount,'H'))'./AAseq_length;
Ifreq = (extractfield(allAAcount,'I'))'./AAseq_length;
Lfreq = (extractfield(allAAcount,'L'))'./AAseq_length;
Mfreq = (extractfield(allAAcount,'M'))'./AAseq_length;
Sfreq = (extractfield(allAAcount,'S'))'./AAseq_length;
Wfreq = (extractfield(allAAcount,'W'))'./AAseq_length;
Yfreq = (extractfield(allAAcount,'Y'))'./AAseq_length;
Nfreq = (extractfield(allAAcount,'N'))'./AAseq_length;

%% train and validation data for linear regression

X = [ones(700,1) tAI Qfreq Hfreq Ifreq Lfreq Mfreq Sfreq Wfreq Yfreq Nfreq AGGFrequency mRNAfold CCGATATfreq];

predictedPA = X*b;

%% sorting genes by predicted PA
% lowest = 1, highest = 700
[~, sortedPAind] = sort(predictedPA);
GeneIndex = (extractfield(unknown_set,'GeneIndex'))';
GeneRank = table([1:700]', GeneIndex(sortedPAind), 'VariableNames', {'Rank', 'GeneIndex'});
