% Funktioner som löser vissa av mozquizto-quizzen finns i botten
% av filen

%% Stänger ner alla fönster och rensar terminalen
clc
close all

%% 2.1
load kroppsTemp.mat;
% whos T; % Helt onödigt

figure(1)
hist(T(:,1));
figure()
plot(T(:,1),'.');

figure()
subplot(211)
hist(T(T(:,2)==1, 1))
subplot(212)
hist(T(T(:,2)==2, 1))
figure()
boxplot(T(:,1), T(:,2)) % Boxplot

% Sorterad data
figure(6)
plot(sort(T(:,1)), 1:size(T,1), '.r')

ratio = mean(T(:,1)<=37);
disp(['Sannolikheten att ha kroppstemp <= 37 är: ', num2str(ratio)]);

% Dessa skriver bara ut boolska arrays. Ointressant. Bättre att
% ta medel av dem
% T(:,1)
% T(:,1)<=37
% T(:,1)<=36.5
% T(:,1)<=38
disp(['Medelkroppstemperaturen är: ', num2str(mean(T(:,1)))])
andel_under(37, T)
andel_under(36.5, T)
andel_under(38, T)
andel_under(36.8, T)
andel_under(37.5, T)
andel_under(36.3, T)

figure()
stairs(sort(T(:,1)),(1:size(T(:,1),1))/size(T(:,1),1))
grid on

%% 2.2

% µ
mu = mean(T(:,1));
% Standardavvikelse
sigma = std(T(:,1));
disp(['Standardavvikelsen för kroppstemperatur är ', num2str(sigma), ' grader'])

figure(1)
hist(T(:,1), 'Normalization', 'pdf') % Normalisera till area == 1
x = linspace(35.5, 38.5, 1e2); % Skapar x-vektor
hold on
plot(x, normpdf(x, mu, sigma))
hold off

figure(6)
stairs(sort(T(:,1)),(1:size(T(:,1),1))/size(T(:,1),1), '-r')
grid on
hold on
plot(x, normcdf(x,mu,sigma))
hold off

%% 2.3

data = normrnd(mu, sigma, 1, 2000); % Simulerar 2000 nya datapunkter med föregående fördelning
figure()
hist(data);
hold on
plot(data, normpdf(data, mu, sigma))
hold off

figure()
stairs(sort(data), (1:length(data))/length(data), '.-g')
hold on
plot(data, normcdf(data, mu, sigma))
hold off
grid on

under_37 = normspec([-Inf 37], mu, sigma);
disp(['Andelen som är under 37 grader är för simuleringen: ', num2str(under_37)])

% Under 36.6 som jag fick i min labb
normspec([-Inf 36.6], mu, sigma)

%% 2.4

x_005 = norminv(1 - 0.05, mu, sigma);
disp(['Kvantilen för 0.05 är: ', num2str(x_005)])

% Kvantilen jag fick för min labb, ändra för ditt värde:
x_004 = norminv(1 - 0.04, mu, sigma);
disp(['Kvantilen för 0.04 är: ', num2str(x_004)])

%% 2.5

x2 = linspace(0, 10, 1000);

% Ritar ut normalfördelning för lite olika µ och sigma
figure()
plot(x2, normpdf(x2, 2, 0.5))
hold on
plot(x2, normpdf(x2,7,0.5))
plot(x2, normpdf(x2, 5, 2))
plot(x2, normpdf(x2, 5, 0.2))
hold off
xlabel('\it x')
title('Täthetsfunktioner, f(x)')
% Lägger till labels
label = 'µ=%g , σ=%g';
legend(sprintf(label, 2, 0.5), sprintf(label, 7, 0.5), sprintf(label, 5, 2), sprintf(label, 5, 0.2))

% Fördelningsfunktioner
figure()
plot(x2, normcdf(x2, 2, 0.5))
hold on
plot(x2, normcdf(x2,7,0.5))
plot(x2, normcdf(x2, 5, 2))
plot(x2, normcdf(x2, 5, 0.2))
hold off
xlabel('x')
title('Fördelningsfunktioner, F(x)')
% Lägger till labels
legend(sprintf(label, 2, 0.5), sprintf(label, 7, 0.5), sprintf(label, 5, 2), sprintf(label, 5, 0.2))

%% 3

u = rand(1000, 1);

% Ritar histogram
figure()
subplot(211)
hist(u)
title('Histogram över 1000 slumpmässiga tal')

% u = 1 - e^(-lambda*x) => x = -1 * ln(1-u) / lambda
% Vi använder detta för att transformera u -> x3
x3 = inverse2exp(u, 3);
subplot(212)
hist(x3)
title({'Histogram över 1000 slumpmässiga tal efter', 'transform till exponentialfördelning, \lambda = 3'})

% För min labb
y = inverse2exp(0.566, 0.87)

%% Funktioner
function andel_under(t, T)
    disp(['Andel under ', num2str(t), ' grader: ', num2str(mean(T(:,1)<=t))])
end

function y = inverse2exp(x, lambda)
    y = -1 .* log(1 - x) / lambda;
    return
end



