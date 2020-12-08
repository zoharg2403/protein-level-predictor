function [AA_frequency] = getAAfreq(known_set,AA)
% this function returnes the frequency of chosen animo acid - 'AA'
% get ORFs from known data set
ORF = (extractfield(known_set,'ORF'))';
% convert NT sequence to AA sequence
AAseq = cellfun(@nt2aa, ORF, 'UniformOutput', false);
AAseq_length = cellfun(@length, AAseq); %% number of AA in AAseq
% counting number of every AA in every AAseq
allAAcount = cellfun(@aacount, AAseq);
% number of chosen AA in every AAseq
AAcount = (extractfield(allAAcount,AA))';  
% calculating AA_frequency for every AAseq
AA_frequency = AAcount./AAseq_length;
end

