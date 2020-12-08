function [avgFE,FE] = mRNAfoldAVG(known_set)
FE = [];
windows = [];
for i = 1:100
    try
        windows = extractfield(known_set,'window'+string(i));
    end
    try
        windows = extractfield(known_set,'wndow'+string(i));
    end
    FE(i,:) = windows';
end
avgFE = mean(FE);
end

