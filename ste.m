function std_error = ste(X,dim)
% Takes standard error of sample mean

if nargin < 2
    dim = 1;
end

if min(size(X)) == 1
    std_error = nanstd(X)./sqrt(length(X));
else
    
    std_error = nanstd(X,0,dim)./sqrt(size(X,dim));
end