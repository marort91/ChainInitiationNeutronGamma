%Plotting Script for Pn,g Functions

clc, clear, clf

%PnData = load('ProbN.txt');
%PgData = load('ProbG.txt');

chains = 1000;

lf = 0.0;
lc = 1 - lf;
ti = 0; tf = 20; 
N = 5;

%DO NOT CHANGE ANYTHING ABOVE THIS LINE

PNDF = load('ProbNStats.txt');
PGDF = load('ProbGStats.txt');

%close all, clf

t = linspace(ti,tf,N+1);

n = 0:1:3;
%n = 0;

for i = 1:length(n)
    
idx = (n(i)+1):31:length(PNDF);
Pn_idx = PNDF(idx,:);
Pn = sum(Pn_idx)./length(Pn_idx);

for j = 1:length(Pn_idx)
stdevn(j,:) = (Pn_idx(j,:)-Pn).^2;
end
stdevn = (1/sqrt(chains)).*sum(stdevn);

figure(1)
%plot(t,PnData(n(i)+1,:),'o');
errorbar(t,Pn,stdevn,'.');
axis([0 max(t) 0 max(Pn)]);
%axis('TIGHT')

strPn = sprintf('P_n%i( \\tau ) ',n(i));
strTitle = sprintf(' Neutrons \\lambda _f = %.2f, \\lambda _c = %.2f, \\tau = 1, %i Chains',lf,lc,chains);

set(gca,'DefaultTextFontSize',18)

xlabel('Mean Neutron Lifetimes \tau','FontSize',16)
ylabel(strPn,'FontSize',16)
title(strTitle,'FontSize',16);
legend('Monte Carlo');

export_fig(sprintf('Neutrons_P%i_sigf%.2f_sigc%.2f_%ichains',n(i),lf,lc,chains), '-pdf', '-m2.5', '-nocrop');

end

g = 0:1:3;

for i = 1:length(g)
    
gidx = g(i)+1:100:length(PGDF);
Pg_idx=PGDF(gidx,:);
Pg = sum(Pg_idx)./length(Pg_idx);

for j = 1:length(Pg_idx)
stdevg(j,:) = (Pg_idx(j,:)-Pg).^2;
end
stdevg = (1/sqrt(chains)).*sum(stdevg);
    
figure(2)
%plot(t,PgData(n(i)+1,:),'o');
%plot(t,Pg,'o');
errorbar(t,Pg,stdevg,'o');
axis([0 20 0 max(Pg)]);
strPg = sprintf('P_g%i( \\tau ) ',g(i));
strTitle = sprintf(' Gammas \\lambda _f = %.2f, \\lambda _c = %.2f, \\tau = 1, %i Chains',lf,lc,chains);
set(gca,'DefaultTextFontSize',18)

xlabel('Mean Neutron Lifetimes \tau','FontSize',16)
ylabel(strPg,'FontSize',16)
title(strTitle,'FontSize',16);
legend('Monte Carlo');

export_fig(sprintf('Gammas_P%i_sigf%.2f_sigc%.2f_%ichains',g(i),lf,lc,chains), '-pdf', '-m2.5', '-nocrop');

end

ntrn = load('ntrnfission.txt');
gamma = load('gammafission.txt');

nubar = sum(ntrn)/length(ntrn);
mubar = sum(gamma)/length(gamma);
