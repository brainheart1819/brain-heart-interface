function [TIndex] = TWavedetection(SIndex, QIndex,signal)
xDiff = differFilter(signal);

TIndex = zeros(1, length(xDiff));
TIndexNo = 0;

for i = 1:length(SIndex)-1
    leftLimit = SIndex(i)+1;
    if QIndex(i) > SIndex(i)
        rightLimit = QIndex(i)+1;
    else
        rightLimit = QIndex(i+1)+1;
    end
    % Find the first maxPeak of xDiff between SIndex and the next QIndex
    [tpeak,locs] = findpeaks(xDiff(leftLimit:rightLimit));
    TIndexNo = TIndexNo + 1;
    TIndex(TIndexNo) = leftLimit+locs(1)+4;
end
TIndex = TIndex(1:TIndexNo);

end