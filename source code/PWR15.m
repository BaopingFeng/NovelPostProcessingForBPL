function [free_sets,frozen_sets] = PWR15(N,K)
load('reliability.mat');   %R15 channel reliability
if N<1024
    val = val(1:N,1);
end
[~,index] = sort(val);
frozen_sets = sort(index(1:N-K));
free_sets = sort(index(K+1:N));
end
 