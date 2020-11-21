function styrkefkn(sigma, n, mu_0, alfa, riktning, mu_sant)
%styrkefkn Illustrerar styrekfunktion for hypotestest
%
% styrkefkn(sigma, n, mu_0, alfa, riktning)
% styrkefkn(sigma, n, mu_0, alfa, riktning, mu_sant)
%
% Illustrerar kritiskt omrade och styrekfunktion for hypotestest av mu under
% antagande om n observationer fran en normalfordelning, N(mu_0, sigma).
%
% Inparametrar:
%   mu_0, sigma - Parametrar for normalfordelningen under H0.
%   n - Antal observationer.
%   alfa - Signifikansniva for testet
%   riktning - Vilken mothypotes att anvanda. H0: mu=mu_0
%      '<' - H1: mu < mu_0
%      '>' - H1: mu > mu_0
%     '!=' - H1: mu != mu_0
%   mu_sant - Varda pa mu att berakna styrkan for, d.v.s.
%             P(forkasta H_0 | mu=mu_sant)
%  
% Exempel:
%   styrkefkn(1, 10, 5, .05, '<', 4.52)

% Copyright 2017 Johan Lindström, Olof Zetterqvist, Lena Zetterqvist; Lunds Universitet.

%Check mu_sant
if nargin<6, mu_sant=[]; end
riktning = validatestring(riktning, {'<', '>', '!='}, ...
  'hypotes', 'riktning', 5);
%replace != with \neq for plotting
if strcmp(riktning,'!='), riktning='\neq'; end

%compute standard deviation of estimate
s = sigma/sqrt(n);

%setup figures
clf; if isempty(mu_sant), subplot(111); else subplot(212); end

%Beräkna och plota styrke funktionen.
if strcmp(riktning,'<')
  k = norminv(alfa, mu_0, s);
  x = linspace(k-3*s, mu_0+s, 1000);
  plot(x, normcdf(k,x,s), 'LineWidth', 2);
elseif strcmp(riktning,'>')
  k = norminv(1-alfa, mu_0, s);
  x = linspace(mu_0-s, k+3*s, 1000);
  plot(x, 1-normcdf(k,x,s), 'LineWidth', 2);
else %if strcmp(riktning,'\neq')
  k1 = norminv(1-alfa/2, mu_0, s);
  k2 = norminv(alfa/2, mu_0, s);
  x = linspace(k2-3.5*s, k1+3.5*s, 1000);
  plot(x, 1-normcdf(k1,x,s)+normcdf(k2,x,s), 'LineWidth', 2);
end
hold off;
grid on;
axis([x(1) x(end) -0.1 1.1]);
xlabel('\mu');
if isempty(mu_sant)
  title('h($\mu$) = P(f\"orkasta H$_0$)', 'interpreter', 'latex')
else
  if strcmp(riktning,'<')
    f = normcdf(k, mu_sant, s);
  elseif strcmp(riktning,'>')
    f = 1-normcdf(k, mu_sant, s);
  else %if strcmp(riktning,'\neq')
    f = 1 - normcdf(k1, mu_sant, s) + normcdf(k2, mu_sant, s);
  end
  title(sprintf('h($\\mu$) = P(f\\"orkasta H$_0$); h(%2.1f) = %2.2f', ...
    mu_sant, f), 'interpreter', 'latex')
  line([mu_sant mu_sant x(1)], [0 f f], 'LineWidth', 2, 'color', [1 0 0]);
end
textning(sigma, n, mu_0, riktning, alfa);

if ~isempty(mu_sant)
  %compute maximum of density function
  f_0 = normpdf(0,0,s);

  %compute density functions
  y = normpdf(x, mu_0, s );
  y_H1 = normpdf(x, mu_sant, s);

  subplot(211)
  plot(x, y_H1, 'b', 'LineWidth', 2)
  hold on;
  plot(x, y, 'k', 'LineWidth', 2)
  axis([min(x) max(x) [0 1.2]*f_0])

  if strcmp(riktning,'\neq')
    x1 = linspace(k1, max(x));
    x2 = linspace(min(x), k2);
    set(area(x1, normpdf(x1,mu_0,s)), 'FaceColor', [1 0 0])
    set(area(x2, normpdf(x2,mu_0,s)), 'FaceColor', [1 0 0])
    x1_H1 = linspace(k2, k1);
  else
    if strcmp(riktning,'>')
      x1 = linspace(k, max(x));
      x1_H1 = linspace(min(x), k);
    else %if strcmp(riktning,'<')
      x1 = linspace(min(x), k);
      x1_H1 = linspace(k, max(x));
    end
    %arean forkasta H0
    set(area(x1, normpdf(x1,mu_0,s)), 'FaceColor', [1 0 0])
  end
  %arean forkasta inte H0
  set(area(x1_H1, normpdf(x1_H1,mu_sant,s)), 'FaceColor', [0 0 1])
  plot([mu_0 mu_0], [0 f_0], 'k--', [mu_sant mu_sant], [0 f_0], 'k--')
  hold off;
  title(sprintf('Sannolikheter f\\"or fel av typ 1 och typ 2; Styrka d{\\aa} $\\mu$ = %4.1f',mu_sant), 'interpreter', 'latex')
end

function textning(sigma, n, mu_0, riktning, alfa)
y = max(get(gca, 'ylim'));
if strcmp(riktning,'>')
  x = mu_0+0.5*sigma/sqrt(n);
elseif strcmp(riktning,'<')
  x = mu_0-sigma/sqrt(n);
else %if strcmp(riktning,'\neq')
  x = mu_0;
end
text(x, 0.9*y, sprintf('H$_0$: $\\mu$ = %4.1f', mu_0), 'interpreter', 'latex')
text(x, 0.8*y, sprintf('H$_1$: $\\mu$ $%s$ %4.1f', riktning, mu_0), 'interpreter', 'latex')
text(x, 0.7*y, sprintf('n = %4.0f',n), 'interpreter', 'latex')
text(x, 0.6*y, sprintf('$\\sigma$ = %4.1f', sigma), 'interpreter', 'latex')
text(x, 0.5*y, sprintf('$\\alpha$ = %4.3f', alfa), 'interpreter', 'latex')
