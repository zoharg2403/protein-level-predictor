clear all; 
%% Load Known data set
known_set = readtable("Known_set_Bacillus.xlsx");
known_set = table2struct(known_set);

%% AA frequency
% for R frequency
R_frequencys = getAAfreq(known_set,'R');

% test
known_R_frequencys = (extractfield(known_set,'argeninFrequnecy'))';
X = sum(known_R_frequencys == R_frequencys);

if X == length(known_R_frequencys)
    disp('getAAfreq correct!!');
else
    disp ('getAAfreq fail :(');
end

%% codon frequency
% for AGG
AGG_frequencys = getCodonFreq(known_set,'AGG');

% test
known_AGG_frequencys = (extractfield(known_set,'AGGFrequency'))';
X = sum(known_AGG_frequencys == AGG_frequencys);

if X == length(known_AGG_frequencys)
    disp('getCodonFreq correct!!');
else
    disp ('getCodonFreq fail :(');
end

%% coding seq length
ORFlengths = ORFlength(known_set);

% test
known_ORF_length = (extractfield(known_set,'ORFLength'))';
X = sum(known_ORF_length == ORFlengths);

if X == length(known_ORF_length)
    disp('ORFlength correct!!');
else
    disp ('ORFlength fail :(');
end

