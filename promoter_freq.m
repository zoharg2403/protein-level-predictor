function [promoter_freq] = promoter_freq(Data,promoter)
% get all ORFs from known data set
ORFs = (extractfield(Data, 'ORF'))';
% counting number of promoter appearances every ORF
promo_ind = cellfun(@(x) strfind(x,promoter), ORFs, 'UniformOutput', false);
promo_count = cellfun(@length,promo_ind);
% normalizing promoter count by ORF length
ORFlength = cellfun(@length , ORFs);
promoter_freq = promo_count./ORFlength;
end

