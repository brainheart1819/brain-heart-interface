function [ SIndex ] = SWavedetection (RIndex,signal)
xDiff = differFilter(signal);

SIndex = zeros(1, length(xDiff));
SIndexNo = 0;

for i = 1:length(RIndex)
    j = RIndex(i);
    if (i == length(RIndex))
        rightLimit = length(xDiff);
    else
        rightLimit = RIndex(i+1);
    end
    % Find the minPeak of xDiff right RIndex
    while ((j < rightLimit) && xDiff(j) > xDiff(j+1))
        j =j+1;
    end
    if (j < rightLimit)
        % Find the zero point
        while ((j < rightLimit) && abs(xDiff(j)) > abs(xDiff(j+1)))
            j = j+1;
        end
        if (j < rightLimit)
            SIndexNo = SIndexNo+1;
            SIndex(SIndexNo) = j;
        end        
    end
end

SIndex = SIndex(1:SIndexNo);

end