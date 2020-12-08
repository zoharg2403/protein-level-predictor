clear all;
%% data preprocess
known_set = readtable("Known_set_Bacillus.xlsx");
known_set = table2struct(known_set);

%% features
% given features
CAI = (extractfield(known_set,'CAI'))';
SDScore = (extractfield(known_set,'SDScore'))';
AGGFrequency = (extractfield(known_set,'AGGFrequency'))';
ORFlength = (extractfield(known_set,'ORFLength'))';

% additional features
tAI = importdata('known_tAI.txt');
tAI = tAI.data(:,2);
mRNAfold = mRNAfoldAVG(known_set)';
GCcontent = GCcontent(known_set);

% AA frequencies
ORFs = (extractfield(known_set,'ORF'))';
AAseq = cellfun(@nt2aa, ORFs, 'UniformOutput', false);
AAseq_length = cellfun(@length, AAseq);
allAAcount = cellfun(@aacount, AAseq);

Afreq = (extractfield(allAAcount,'A'))'./AAseq_length;
Rfreq = (extractfield(allAAcount,'R'))'./AAseq_length;
Nfreq = (extractfield(allAAcount,'N'))'./AAseq_length;
Dfreq = (extractfield(allAAcount,'D'))'./AAseq_length;
Cfreq = (extractfield(allAAcount,'C'))'./AAseq_length;
Qfreq = (extractfield(allAAcount,'Q'))'./AAseq_length;
Efreq = (extractfield(allAAcount,'E'))'./AAseq_length;
Gfreq = (extractfield(allAAcount,'G'))'./AAseq_length;
Hfreq = (extractfield(allAAcount,'H'))'./AAseq_length;
Ifreq = (extractfield(allAAcount,'I'))'./AAseq_length;
Lfreq = (extractfield(allAAcount,'L'))'./AAseq_length;
Kfreq = (extractfield(allAAcount,'K'))'./AAseq_length;
Mfreq = (extractfield(allAAcount,'M'))'./AAseq_length;
Ffreq = (extractfield(allAAcount,'F'))'./AAseq_length;
Pfreq = (extractfield(allAAcount,'P'))'./AAseq_length;
Sfreq = (extractfield(allAAcount,'S'))'./AAseq_length;
Tfreq = (extractfield(allAAcount,'T'))'./AAseq_length;
Wfreq = (extractfield(allAAcount,'W'))'./AAseq_length;
Yfreq = (extractfield(allAAcount,'Y'))'./AAseq_length;
Vfreq = (extractfield(allAAcount,'V'))'./AAseq_length;

% promotor frequencies
CTAAAfreq = promoter_freq(known_set,'CTAAA');
CCGATATfreq = promoter_freq(known_set,'CCGATAT');
TATAATfreq = promoter_freq(known_set,'TATAAT');

%% train and validation data for linear regression
max_b = [];
max_rho = 0;
max_pval = 0;
for i = 1:100

    [trainInd,valInd,~] = dividerand(2775,0.9,0.1,0);

    X = [ones(2775,1) tAI Qfreq Hfreq Ifreq Lfreq Mfreq Sfreq Wfreq Yfreq Nfreq AGGFrequency mRNAfold CCGATATfreq];
    trainX = X(trainInd,:);
    valX = X(valInd,:);

    Y = (extractfield(known_set,'PA'))';
    Y = cellfun(@str2num,Y);
    trainY = Y(trainInd,:);
    valY = Y(valInd,:);

%% train
    b = regress(trainY, trainX);

    predictedPA = valX*b;

%% corolation
    [rho,pval] = corr(predictedPA,valY,'Type','Spearman');
    if rho>max_rho
        max_rho = rho;
        max_pval = pval;
        max_b = b;
    end
end