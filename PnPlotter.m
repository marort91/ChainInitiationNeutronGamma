%Plotting Script for Pn,g Functions

clc, clear, clf

PnData = load('ProbN.txt');
PgData = load('ProbG.txt');

chains = 100000;

lf = 0.25;
lc = 1 - lf;
ti = 0; tf = 20; 
N = 50;

%DO NOT CHANGE ANYTHING ABOVE THIS LINE

%close all, clf

t = linspace(ti,tf,N+1);

%n = 1:1:11;
n = 0:1:7;

for i = 1:length(n)

figure(1)
plot(t,PnData(n(i)+1,:),'o');

strPn = sprintf('P_n%i( \\tau ) ',n(i));
strTitle = sprintf(' Neutrons \\lambda _f = %.2f, \\lambda _c = %.2f, \\tau = 1, %i Chains',lf,lc,chains);

set(gca,'DefaultTextFontSize',18)

xlabel('Mean Neutron Lifetimes \tau','FontSize',16)
ylabel(strPn,'FontSize',16)
title(strTitle,'FontSize',16);
legend('Monte Carlo');

export_fig(sprintf('Neutrons_P%i_sigf%.2f_sigc%.2f_%ichains',n(i),lf,lc,chains), '-pdf', '-m2.5', '-nocrop');

end

for i = 1:length(n)
    
figure(2)
plot(t,PgData(n(i)+1,:),'o');
strPg = sprintf('P_g%i( \\tau ) ',n(i));
strTitle = sprintf(' Gammas \\lambda _f = %.2f, \\lambda _c = %.2f, \\tau = 1, %i Chains',lf,lc,chains);
set(gca,'DefaultTextFontSize',18)

xlabel('Mean Neutron Lifetimes \tau','FontSize',16)
ylabel(strPg,'FontSize',16)
title(strTitle,'FontSize',16);
legend('Monte Carlo');

export_fig(sprintf('Gammas_P%i_sigf%.2f_sigc%.2f_%ichains',n(i),lf,lc,chains), '-pdf', '-m2.5', '-nocrop');

end

ntrn = load('ntrnfission.txt');
gamma = load('gammafission.txt');

nubar = sum(ntrn)/length(ntrn);
mubar = sum(gamma)/length(gamma);
