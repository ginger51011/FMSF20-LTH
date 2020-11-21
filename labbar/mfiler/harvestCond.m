function [X,P] = harvestCond(n,p,mu,y)
% HARVESTCOND Plot conditional density for binomial sum of Poisson
%
% harvestCond(n,p,mu,y)
%
% Plots the conditional density for the number of original seeds that grew
% (X) given the number of harvested seeds (y). Assuming a model where
%   X ~ Bin(n,p)
% and, assuming X=k,
%   Y = sum_i^k Z_i ~ Po(k*mu)
% since number of harvested seeds for each original seed that grew is
%   Z_i ~ Po(mu)
%
% Alternative use
%   [X,P] = harvestCond(n,p,mu,y)
% which returns the computed conditional probabilities with corresponding X
% vector. Example:
%   [X,P] = harvestCond(10,0.75,10,50);
%   bar(X, P)    

% Johan Lindström

%%betingad p_X|Y
x = 0:n;
%what should be computed
% pX_Y = binopdf(x,n,p).*poisspdf(y,x*mu);
%However, we use a numerically stable log approach
pX_Y = log(binopdf(x,n,p)) - x*mu;
if y~=0 %for y!=0 we also need to add (x*mu)^y.
  pX_Y = pX_Y + y*log(x*mu);
  %this term is y*log(x*mu) = 0*log(0) = nan for the case of y=0 and x=0
  %which gives numerical issues in remaining computations.
end
%convert to standard scale
pX_Y = exp(pX_Y-max(pX_Y));
%and normalize
pX_Y = pX_Y / sum(pX_Y);
% Plot density
if nargout==0
    bar(x, pX_Y)
    xlabel('x')
    title( sprintf('p_{X|Y}(k | Y=%u)', y) )
    axis([0-0.5 n+0.5 0 Inf])
else
    X = x;
    P = pX_Y;
end
