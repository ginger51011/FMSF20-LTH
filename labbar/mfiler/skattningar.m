function skattningar(mu, sigma, n1, n2, alternativ)
%skattningar Illustrerar mu och sigma^2-skattning samt konfidensintervall
%
% skattningar(mu, sigma, n1, n2, alternativ)
%
% Ritar histogram for mu- och sigma^2-skattning samt illustrerar
% konfidensintervall for mu. Skattningarna baseras pa n1 respektive n2
% observationer fran en normalfordelning, N(mu, sigma).
%
% Inparametrar:
%   mu, sigma - Parametrar i normalfordelningen.
%   n1, n2 - Antal observationer att simulera i de tva stickproven. Gor det
%            mojligt att jamfora hur skattningarna beter sig med olika
%            antal observationer.
%  alternativ - Textstrang som talar om vad som ska illustreras:
%       'muskatt': Histogram for skattningar av mu
%    'sigmaskatt': Histogram for skattningar av sigma^2
%       'konfint': Illustrerar konfidensintervall for mu.
%          'alla': Samtligga tre plottar
%  
% Exempel:
%   skattningar(35,2,5,25,'muskatt')

% Copyright 2017 Johan Lindström, Olof Zetterqvist, Lena Zetterqvist; Lunds Universitet.

%Check alternativ
if nargin<5, alternativ=[]; end
if ~isempty(alternativ)
  alternativ = validatestring(alternativ, ...
    {'alla', 'muskatt', 'sigmaskatt', 'konfint'}, ...
    'skattningar', 'alternativ', 5);
end
if strcmp(alternativ,'alla'), alternativ=[]; end

%Antal simuleringar som gors
n = 1000;
%simulera tva sample
x = normrnd(mu,sigma,n1,n);
y = normrnd(mu,sigma,n2,n);

if isempty(alternativ) || strcmp(alternativ,'muskatt')
  figure;
  subplot(211)
  hist(mean(x), 20)
  
  %intervallens bredd, for att satta axlar.
  width = 3.5*sigma/sqrt(min(n1,n2));
  axis([mu-width mu+width 0 200])
  hold on;
  plot([mu mu],[1 200],'r')
  hold off;
  %title(sprintf('n obs'))
  title(['$\mu$-skattningens f\"ordelning, n1 = ', sprintf('%2.0f',n1)],...
    'interpreter', 'latex');
  subplot(212)
  hist(mean(y),20)
  axis([mu-width mu+width 0 200])
  %title(sprintf('n obs'))
  title(['$\mu$-skattningens f\"ordelning, n2 = ', sprintf('%2.0f',n2)],...
    'interpreter', 'latex');
  hold on;
  plot([mu mu],[1 200],'r')
  hold off;
end

% sigma-skattningar
if isempty(alternativ) || strcmp(alternativ,'sigmaskatt')
  figure
  subplot(211)
  hist(var(x),20)
  width = sigma^2 * max(chi2inv(0.9995, n1-1)/(n1-1), ...
    chi2inv(0.9995, n2-1)/(n2-1));
  axis([0 width 0 200])
  title(['$\sigma^2$-skattningens f\"ordelning, n1 = ', ...
    sprintf('%2.0f',n1)], 'interpreter', 'latex');
  hold on;
  plot([sigma sigma].^2,[1 200],'r')
  hold off;

  subplot(212)
  hist(var(y),20)
  axis([0 width 0 200])
  hold on;
  plot([sigma sigma].^2,[1 200],'r')
  hold off;
  title(['$\sigma^2$-skattningens f\"ordelning, n2 = ', ...
    sprintf('%2.0f',n2)], 'interpreter', 'latex');
end

%SIMULERING AV KONFIDENSINTERVALL

if isempty(alternativ) || strcmp(alternativ,'konfint')
  % berakna de 1000 konfidens intervallen
  CI_x = bsxfun(@plus, mean(x), norminv(0.975)*sigma/sqrt(n1)*[-1;1]);
  CI_y = bsxfun(@plus, mean(y), norminv(0.975)*sigma/sqrt(n2)*[-1;1]);
  %testa vilka intervall som tacker mu
  I_x = CI_x(1,:)<mu & mu<CI_x(2,:);
  I_y = CI_y(1,:)<mu & mu<CI_y(2,:);
  %computer intervall widths for plotting
  width = 1.2*max([max(abs(CI_x(:,1:100)-mu)) max(abs(CI_y(:,1:100)-mu))]);

  %use only the first 100 for plotting
  IndGood = find(I_x(1:100));
  IndBad = find(~I_x(1:100));
  
  figure
  subplot(121)
  line(CI_x(:,IndGood), [IndGood;IndGood], 'color', 'b')
  line(CI_x(:,IndBad), [IndBad;IndBad], 'color', 'r')
  line([mu mu],[0 101], 'color', 'k')
  axis([mu-width mu+width 0 101])
  % markera det sanna vardet
  title({['100 intervall f\"or $\mu$, n1 =', sprintf(' %2.0f',n1)], ...
    'konfidensgrad = 0.95'}, 'interpreter', 'latex');
  xlabel({'Andel av 1000 intervall som', ...
    ['missar sant $\mu$: ', sprintf('%1.3f',mean(~I_x))]}, ...
    'interpreter', 'latex');

  %use only the first 100 for plotting
  IndGood = find(I_y(1:100));
  IndBad = find(~I_y(1:100));
  
  subplot(122)
  line(CI_y(:,IndGood), [IndGood;IndGood], 'color', 'b')
  line(CI_y(:,IndBad), [IndBad;IndBad], 'color', 'r')
  line([mu mu],[0 101], 'color', 'k')
  axis([mu-width mu+width 0 101])
  % markera det sanna vardet
  title({['100 intervall f\"or $\mu$, n2 =', sprintf(' %2.0f',n2)], ...
    'konfidensgrad = 0.95'}, 'interpreter', 'latex');
  xlabel({'Andel av 1000 intervall som', ...
    ['missar sant $\mu$: ', sprintf('%1.3f',mean(~I_y))]}, ...
    'interpreter', 'latex');
end
