%Data Analysis

close all

lf = 0.0;
lc = 1 - lf;

ProbsN = load('ProbNStats.txt');
ProbsG = load('ProbGStats.txt');

NtrnFiss = load('ntrnfission.txt');
GammaFiss = load('gammafission.txt');

y = load('ntrn.mult');
v = 0:1:7;
nubaract = sum(y.*v);
%y = load('gamma.mult');
%mubaract = sum(y.*v)

nubar = sum(NtrnFiss)/length(NtrnFiss);
mubar = sum(GammaFiss)/length(GammaFiss);

ntrndev = (1/(sqrt(5e4)))*sum((y-nubar).^2);
lens = length(ProbsN);
chains = 1000;
N = 5;
ti = 0; tf = 20;
t = linspace(ti,tf,N+1);
Pn = 1;

idxn = 1+Pn:11:lens;

PnMatrix = ProbsN(idxn,:);

PnSums = sum(PnMatrix)/length(idxn);

dims = size(PnMatrix);

for i = 1:dims(1)
    
    diffmat(i,:) = PnMatrix(i,:) - PnSums;
    
end

diffmatsqc = diffmat.^2;
diffmatsqcsum = sum(diffmatsqc);
stddev = (1/sqrt(chains)).*diffmatsqcsum;

%plot(t,PnSums,'o')
%errorb(t,PnSums,stddev)
%errorbar(t,PnSums,stddev,'.')

set(gca,'DefaultTextFontSize',18)

figure(1)

errorbar(t,PnSums,stddev,'o-')
axis([0 20 0 1])

strPn = sprintf('P_n%i( \\tau ) ',Pn);
strTitle = sprintf(' \\lambda _f = %.2f, \\lambda _c = %.2f, \\tau = 1, %i Chains',lf,lc,chains);
xlabel('Mean Neutron Lifetimes \tau','FontSize',16)
ylabel(strPn,'FontSize',16)
title(strTitle,'FontSize',16);
legend('Monte Carlo');

Pg = 0:1:25;

for i = 1:length(Pg)

idxg = 1+Pg(i):26:length(ProbsG);

PgMatrix = ProbsG(idxg,:);

PgSums(i,:) = sum(PgMatrix)/length(idxg);

dims = size(PgMatrix);

for j = 1:dims(1)
    
    diffmatg(j,:) = PgMatrix(j,:) - PgSums(i,:);
    
end

diffmatsqcg = diffmatg.^2;
diffmatsqcsumg = sum(diffmatsqcg);
stddevg = (1/sqrt(chains)).*diffmatsqcsumg;

hold on

figure(2)
%errorbar(t,PgSums,stddevg,'o')
plot(t,PgSums,'o');
axis([0 20 0 1])

strPn = sprintf('P_g%i( \\tau ) ',Pg);
strTitle = sprintf(' \\lambda _f = %.2f, \\lambda _c = %.2f, \\tau = 1, %i Chains',lf,lc,chains);
xlabel('Mean Neutron Lifetimes \tau','FontSize',16)
ylabel(strPn,'FontSize',16)
title(strTitle,'FontSize',16);
legend('Monte Carlo');

end
