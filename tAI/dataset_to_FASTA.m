clear all;

unknown_set = readtable("Unknown_set_Bacillus.xlsx");
unknown_set = table2struct(unknown_set);

H = (extractfield(unknown_set,'GeneIndex'))';
init = repmat({'>'},length(H),1);
Header = cellfun(@(x,y) horzcat(x,y),init , H, 'UniformOutput', false);
seq = (extractfield(unknown_set,'ORF'))';
fastawrite('my_fasta.txt', Header, seq);