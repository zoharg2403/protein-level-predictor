function [seq_length] = ORFlength(known_set)
% get all ORFs from known data set
ORF = (extractfield(known_set,'ORF'))';
% get length for every ORF
seq_length = cellfun(@length, ORF);
end

