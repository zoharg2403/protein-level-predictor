function [GCcontent] = GCcontent(data)
% get ORFs from data and ORFs length
ORF = (extractfield(data,'ORF'))';
orfLength = cellfun(@length, ORF);
% count all bases
BaseCount = cellfun(@basecount, ORF);
% C and G number of appeariences
Ccount = (extractfield(BaseCount,'C'))';
Gcount = (extractfield(BaseCount,'G'))';
% GCcontent = (#G+#C)/totalNTnumber
GCcontent = (Ccount + Gcount)./orfLength;
end

