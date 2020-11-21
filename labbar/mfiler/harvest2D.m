function [X, Y, pXY] = harvest2D(n,p,mu)
% HARVEST2D Plot joint density for binomial sum of Poisson
%
% harvest2D(n,p,mu)
%
% Plots the joint density for the number of harvested seeds (Y) and the
% number of original seeds that grew (X). Assuming a model where
%   X ~ Bin(n,p)
% and, assuming X=k,
%   Y = sum_i^k Z_i ~ Po(k*mu)
% since number of harvested seeds for each original seed that grew is
%   Z_i ~ Po(mu)
%
% Alternative use
%   [X, Y, pXY] = harvest2D(n,p,mu)
% which returns the computed joint probabilities with corresponding X and Y
% vectors. Example:
%   [X, Y, pXY] = harvest2D(7, 0.4, 7);    
%   stem3(X, Y, pXY)    

% Johan Lindström

% Determin Ymax, based on E(Y) and V(Y)
tot_E=n*p*mu;
tot_V=tot_E*(1+mu*(1-p));
Ymax=tot_E+4*sqrt(tot_V);
x = 0:n;
y = 0:Ymax;

% Compute probabilities using Law of Total Prob.
pXY_in = bsxfun(@(x,y) binopdf(x,n,p).*poisspdf(y,x*mu), x, y');

X_in = repmat(x,[length(y) 1]);
Y_in = repmat(y',[1 length(x)]);

%plots
if nargout==0
    subplot(121)
    stem3(x,y,pXY_in, 'marker', 'none')
    xlabel('x')
    ylabel('y')
    zlabel('p_{XY}(x,y)')
    axis([0 n 0 Ymax 0 Inf])
    
    subplot(122)
    scatter(X_in(:), Y_in(:), 25, pXY_in(:), 'filled')
    colorbar
    xlabel('x')
    ylabel('y')
    title('p_{XY}(x,y)')
    axis tight
else
    X = X_in;
    Y = Y_in;
    pXY = pXY_in;
end