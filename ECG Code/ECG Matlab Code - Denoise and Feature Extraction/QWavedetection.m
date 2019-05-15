function [QIndex] = QWavedetection(RIndex,signal)
%xDiff = differFilter(signal);
xDiff = differFilter(signal);

%RIndex = rlocaiton;
QIndex = zeros(1, length(xDiff));
QIndexNo = 0;

for i = 1:length(RIndex)
    j = RIndex(i);
    if (i == 1)
        leftLimit = 0;
    else
        leftLimit = RIndex(i-1);
    end
    % Find the minPeak of xDiff left RIndex
    while ((j > leftLimit) && xDiff(j) < xDiff(j-1))
        j =j-1;
    end
    if (j > leftLimit)
        % Find the zero point
        while ((j > leftLimit) && abs(xDiff(j)) > abs(xDiff(j-1)))
            j = j-1;
        end
        if (j > leftLimit)
            QIndexNo = QIndexNo+1;
            QIndex(QIndexNo) = j;
        end        
    end
end
QIndex = QIndex(1:QIndexNo);