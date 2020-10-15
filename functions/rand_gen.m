function [ y  ]=rand_gen(F,x,n)

% the integral of PDF is the cumulative distribution function
%x = interp(x,4);
%F = interp(F,4);

% remove negative elements due to the spline interpolation
%F(F < 0) = 0;
cdf = cumsum(F);

%cdf = interp(cdf,4);

% remove non-unique elements
[cdf, mask] = unique(cdf);
x = x(mask);

%%
% create an array of 2500 random numbers
randomValues = rand(1, n);

% inverse interpolation to achieve P(x) -> x projection of the random values
y = interp1(cdf, x, randomValues);

end