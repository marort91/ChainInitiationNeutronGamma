%Plotting Script for Pn,g Functions

clc, clear, clf

PnData = load('ProbN.txt');
PgData = load('ProbG.txt');

chains = 1e3;

lf = 0; lc = 1 - lf;
ti = 0; tf = 20; N = 51;

n = 1;

t = linspace(ti,tf,N);

figure(1)
plot(t,PnData(n,:),'o');

strPn = sprintf('P_n%i( \\tau ) ',n-1);
strTitle = sprintf(' \\lambda _f = %.2f, \\lambda _c = %.2f, \\tau = 1, %i Chains',lf,lc,chains);

set(gca,'DefaultTextFontSize',18)

xlabel('Mean Neutron Lifetimes \tau','FontSize',16)
ylabel(strPn,'FontSize',16)
title(strTitle,'FontSize',16);
legend('Monte Carlo');

figure(2)
plot(t,PgData(n,:),'o');
strPg = sprintf('P_g%i( \\tau ) ',n-1);
set(gca,'DefaultTextFontSize',18)

xlabel('Mean Neutron Lifetimes \tau','FontSize',16)
ylabel(strPg,'FontSize',16)
title(strTitle,'FontSize',16);
legend('Monte Carlo');
