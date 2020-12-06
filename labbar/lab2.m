% Funktioner som löser vissa av mozquizto-quizzen finns i botten
% av filen

% Förvarning: Denna labben är helt kaos, men jag tänker inte sitta
% och fixa en labb jag fick godkännt på

%% Stänger ner alla fönster och rensar terminalen
clc
close all

%% 1.2

% Mozquizto 1: X tillhör Bin(7, 0.75)
disp(['1.2.2 E(X)=', num2str(7 * 0.75), ' V(X)=', num2str(7*0.75*(1-0.75))])
% Mozquizto 2: Y = Po(µ_1+µ_2+...)

%% 2.1
% Ändrat för mitt Mozquizto
n = 5;
p = 0.6
U = rand(1, n);
X = sum(U<=p);

figure(1)
stem(U)
refline(0, p)

X = binornd(n, p);

N = 100;
X = binornd(n, p, N, 1);
antal = hist(X, 0:n);
disp(['Antal där X=3: ', num2str(antal(4))])
sum(X==3)
figure(2)
bar(0:n, antal)
xlabel('Antal frön som gror')
ylabel('Antal tillfällen')

bar(0:n,[antal/N; binopdf(0:n,n,p)]) 
xlabel('Antal frön som gror')
ylabel('Andel påsar')
legend('Simulering', 'Exakt')

% För att få n: Kolla vad högsta värdet är i grafen
% För p: Kolla vad mittenvärdet är och kör mittenvärdet / n

%% 2.2
% Ändra talen för att matcha uppgiften du får av Mozquizto

n = 85;
p = 0.8;
np = n * p;
npq = np * (1 - p);
x = linspace(np-4*sqrt(npq),np+4*sqrt(npq)); % µ +- sigma

figure(3)
stairs(0:n,binocdf(0:n,n,p))
hold on
plot(x,normcdf(x,np,sqrt(npq)))
hold off

% Mozquizto 7 - Ändra siffror och läs av i grafen
%disp('Mozquizto 7 - Sannolikheten att högst z frön gror')

%exact_func = @(x) binocdf(x, n, p);
%exact_prob = integral(exact_func, 0, 75)

%norm_func = @(x) normcdf(x, np, sqrt(npq));
%norm_prob = integral(norm_func, 0, 75)

%% 2.3

figure(4)
subplot(211)
bar(0:2, binopdf(0:2,2,.75))
title('Antal frön som gror')
ylabel('\itp(x)')

mu = 10;
x = 0:4*mu;

figure(4)
subplot(234)
bar(x, fromplanted(x, 0, mu))
title('Skörd om 0 frön gror')
ylabel('\itp(y|x=0)')

subplot(235)
bar(x, fromplanted(x, 1, mu))
title('Skörd om 1 frö gror')
ylabel('\itp(y|x=1)')

subplot(236)
bar(x, fromplanted(x, 2, mu))
title('Skörd om 2 frön gror')
ylabel('\itp(y|x=2)')

figure(5)
bar(x, probgrow(x, 2, mu, 0.75))
xlabel('Antal nya frön')

y = 0:100;
figure(6)
bar(y, probgrow(y, 7, 10, 0.75))
xlabel('Antal nya frön')

harvest(7, 0.75, 10)

%% 3

%% test
mu = 3;
x = poissrnd(mu, 2, 1);
disp(['Medel är ', num2str(mean(x))])

%% Drrrr
n = 5;
py = harvest(5, 0.6, 3);
E = 0;
for k=0:n
    E = E + k * py(k);
end
disp(num2str(E))

%% Drrr2
X = poissrnd(1.7, 34, 1);
xmedel = mean(X)

%% 3 igen

mu = 1.7;
n = 34;
M = 10000;
x = poissrnd(mu, n, M);
xmedel = mean(x);
subplot(2, 1, 1)
hist(x(1,:), 0:15)
subplot(2,1,2)
hist(xmedel, 0:0.01:15);

%% Fråga 9
[Y, P] = harvest(5, 0.6, 3);
sum(Y .* P)

%% Fråga 10
% För Poisson gäller µ + µ + µ... Dvs. ta så många µ
% att sum(µ) > 15

% För binomialfördelning gäller Bin(n1+n2, p), dvs.
% vi vill att (n1+n2+n2... * p * (1-p)) >= 10
% => vi löser x*5 * 0.2 * 0.8 >= 10 => x=13
%% Funktioner

% Världens mest onödiga funktion? Kanske
function y = fromplanted(x, k, mu)
    y = poisspdf(x, k * mu);
end

function answer = question9(n, p, mu)
    [Y, P] = harvest(5, 0.6, 3);
    answer = sum(Y .* P);
end

function pY = probgrow(x, n, mu, p)
    pY = recursive_probgrow(x, n, mu, n, p);

    function pY = recursive_probgrow(x, k, mu, n, p)
        if k <= 0
            pY = poisspdf(x, n * mu) * binopdf(0, n, p);
        else
            pY = recursive_probgrow(x, k-1, mu, p, n) + poisspdf(x, n * mu) * binopdf(n, n, p);
        end
    end
end


