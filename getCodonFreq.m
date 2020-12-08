function [CodonFrequencys] = getCodonFreq(known_set,codon)
% this function returnes the frequency of chosen codon - "codon"
% get ORFs from known data set and their length in codons
ORF = (extractfield(known_set,'ORF'))';
totalCodonNum = cellfun(@(x) length(x)/3,ORF); % ORFs length in codons
% count all codons in every ORF
allCodonCount = cellfun(@codoncount, ORF);
% number of chosen codon in every ORF
CodonCount = (extractfield(allCodonCount,codon))';
% calculating codon frequency for every ORF
CodonFrequencys = CodonCount./totalCodonNum;
end

