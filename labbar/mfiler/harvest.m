function [Y,P] = harvest(n,p,mu)
% HARVEST Plot marginal density for binomial sum of Poisson
%
% harvest(n,p,mu)
%
% Plots the marginal density for the number of harvested seeds (Y).
% Assuming a model where
%   X ~ Bin(n,p)
% and, assuming X=k,
%   Y = sum_i^k Z_i ~ Po(k*mu)
% since number of harvested seeds for each original seed that grew is
%   Z_i ~ Po(mu)
%
% Alternative use
%   [Y, P] = harvest(n,p,mu)
% which returns the computed marignal probabilities with corresponding Y
% vector. Example:
%   [Y, P] = harvest2D(7, 0.4, 7);
%   bar(Y, P)    

% Johan Lindström (modified from originial file by Anna Lindgren)

% Determin Ymax, based on E(Y) and V(Y)
tot_E=n*p*mu;
tot_V=tot_E*(1+mu*(1-p));
Ymax=tot_E+4*sqrt(tot_V);
y = (0:Ymax)';

% Compute probabilities using Law of Total Prob.
pY=zeros(size(y));
for k=0:n
    pY=pY+poisspdf(y,k*mu)*binopdf(k,n,p);
end

% Plot density
if nargout==0
    bar(y,pY)
    titel=['Sannolikhetsfunktionen för antal nya frön, n='...
        num2str(n) ', p=' num2str(p) ', \mu=' num2str(mu)];
    xlabel('antal nya frön (k)')
    ylabel('p(k)')
    title(titel)
    axis([0 Ymax 0 Inf])
else
    Y = y;
    P = pY;
end