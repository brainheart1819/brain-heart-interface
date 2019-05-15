function [ y ] = differFilter( x )
%differFilter Summary of this function goes here
%   Detailed explanation goes here

% The 2-index delay makes larger mistake
% But can less the radom difference
% numd = [1 0 -1];
% dend = [1 0 0];

% Hear we take the 1-index instead
numd = [1 -1];
dend = [1 0];

y = filter(numd, dend, x);

end

