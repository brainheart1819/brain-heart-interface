function [ RIndex ] = RPeekDetect( x, Ts )
%RPeekDetect Summary of this function goes here
%   Detailed explanation goes here

xDiff = differFilter(x);

% Variables used
RIndexNo = 0; % Storing R-Peek Index No(Count)
RIndex = zeros(1, length(xDiff)); % Storing R-Peek Index
minThresholdAry = zeros(1, 6);
rrIntervalAry = zeros(1, 5);

% Get the prefixed sampled threshold
minTemp = zeros(1, 10);
for i = 1:10
    beginIndex = (i-1)*3/Ts + 1; % For 3 Seconds divided into 10 groups
    endIndex = i*3/Ts;
    array = xDiff(beginIndex:endIndex);
    minTemp(i) = min(array);    
end
minThreshold = 0.5 * (sum(minTemp) - max(minTemp) - min(minTemp)) / (length(minTemp)-2);

% Get the first 6 R-Peeks
j = 1;
while ((RIndexNo < 6) && (j < length(xDiff)))
    if (xDiff(j) < minThreshold)
        jTemp = j;
        while ((jTemp>0) && (abs(xDiff(jTemp))>abs(xDiff(jTemp-1))))
            jTemp = jTemp-1;
        end
        if (jTemp > 0)
            RPeekIndex = jTemp;
            jTemp = j;
            while ((jTemp<length(xDiff)) && (abs(xDiff(jTemp))...
                    < abs(xDiff(jTemp+1))))
                jTemp = jTemp+1;
            end
            if(jTemp < length(xDiff))
                RIndexNo = RIndexNo+1;
                RIndex(RIndexNo) = RPeekIndex;
                minThresholdAry(RIndexNo) = xDiff(jTemp);
                if (RIndexNo > 1)
                    rrIntervalAry(RIndexNo-1) = RIndex(RIndexNo)-RIndex(RIndexNo-1);
                end
            end
            j = jTemp + int32(0.15/Ts);
        else
            % We know that the max heart rate is smaller than 250Beats/min
            % Which the smallest RR interval is 0.24s
            % And the max QRS length is 0.1s
            % Here we set add j with 0.15/Ts, which maks the algorithm
            % faster
            while ((jTemp<length(xDiff)) && (abs(xDiff(jTemp))...
                    < abs(xDiff(jTemp+1))))
                jTemp = jTemp+1;
            end
            j = jTemp + int32(0.15/Ts);
        end
    end
    j = j+1;
end

% Calculate the RR-Interval and minThreshold
minThreshold = 0.5 * ((sum(minThresholdAry)-max(minThresholdAry)-min(minThresholdAry))/(length(minThresholdAry)-2));
rrInterval = sum(rrIntervalAry)/length(rrIntervalAry);

m = RIndex(RIndexNo) + int32(0.15/Ts); % Start from 7th point
repeatTime = 0;
doSkip = 0; % 0-Normal; 1-Too long;
while (m < length(xDiff))
    if (xDiff(m) < minThreshold)
        mTemp = m;
        while ((abs(xDiff(mTemp))>abs(xDiff(mTemp-1))))
            mTemp = mTemp-1;
        end
        RPeekIndex = mTemp;
        if ((RPeekIndex-RIndex(RIndexNo)) > 1.66*rrInterval && repeatTime < 6)
            repeatTime = repeatTime+1;
            minThreshold = minThreshold * 0.8^repeatTime;
            m = RIndex(RIndexNo) + int32(0.15/Ts);
            doSkip = 1;
        elseif ((RPeekIndex-RIndex(RIndexNo)) < 0.6*rrInterval && doSkip == 0)
            repeatTime = 0;
            mTemp = m;
            while (mTemp<length(xDiff) && xDiff(mTemp) < minThreshold)
                mTemp = mTemp+1;
            end
            m = mTemp + int32(0.15/Ts);
        else
            diSkip = 0;
            repeatTime = 0;
            RIndexNo = RIndexNo+1;
            RIndex(RIndexNo) = RPeekIndex;
            minThresholdAry(1:(length(minThresholdAry)-1))...
                = minThresholdAry(2:(length(minThresholdAry)));
            mTemp = m;
            while ((mTemp<length(xDiff)) && (abs(xDiff(mTemp))...
                    < abs(xDiff(mTemp+1))))
                mTemp = mTemp+1;
            end
            if (mTemp < length(xDiff))
                minThresholdAry(length(minThresholdAry)) = xDiff(mTemp);
            else
                minThresholdAry(length(minThresholdAry)) = 0;
            end
            rrIntervalAry(1:(length(rrIntervalAry)-1))...
                = rrIntervalAry(2:length(rrIntervalAry));
            rrIntervalAry(length(rrIntervalAry)) = RIndex(RIndexNo)...
                - RIndex(RIndexNo-1);
            minThreshold = 0.5 * ((sum(minThresholdAry)-max(minThresholdAry)-min(minThresholdAry))/(length(minThresholdAry)-2));
            rrInterval = sum(rrIntervalAry)/length(rrIntervalAry);
            m = RIndex(RIndexNo) + int32(0.15/Ts);
        end
    end
    m = m+1;
end

RIndex = RIndex(1:RIndexNo);

end