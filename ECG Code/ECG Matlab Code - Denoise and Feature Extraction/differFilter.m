function [ y ] = differFilter( x )

numd = [1 -1];
dend = [1 0];

y = filter(numd, dend, x);

end

