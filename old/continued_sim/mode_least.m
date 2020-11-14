function [M] = mode_least(x)
uv = unique(x);
n = histc(x,uv);
[~,i] = min(n);
M = uv(i);
end