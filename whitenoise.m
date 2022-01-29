function [y, t] = whitenoise(L, var, fs)

y = var*randn(L,1);
t = (1:L).*(1/fs);

end