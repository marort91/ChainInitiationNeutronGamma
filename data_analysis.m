clc, clear, clf
%clf

data = load('fortran.txt');
%load 'for1e6.mat';

dim = size(data);

chains = 2e4;
lf = 0.3;
lc = 1 - lf;

N = 0:1:7;
%N = 1;
%N = 1;

for k = 1:length(N)
    
clf
%time = 0:0.2525:25+0.2525;
time = linspace(0,20,1001);

for i = 1:dim(2)
    
    Pn(i) = length(data(data(:,i)==N(k))); %/dim(1);
    
end

n = 0:1:10;

for j = 1:length(n)
    
    Prob(j) = length(data(data(:,2)==n(j)))/dim(1);
    
end

%Pn(29:end) = 0;

%figure(1)
hold on

%area = trapz(time,Pn./dim(1));

y = Pn./dim(1); %/area;

strPn = sprintf('P_%i( \\tau ) ',N(k));
strTitle = sprintf(' \\lambda _f = %.2f, \\lambda _c = %.2f, \\tau = 1, %i Chains',lf,lc,chains);

plot(time,y,'o-')
set(gca,'DefaultTextFontSize',18)

if N(k) ~= 0
    %str = sprintf('Area = %f',area);
    %annotation('textbox',[0.4,0.6,0.1,0.1],'String', str);
    
end

xlabel('Mean Neutron Lifetimes \tau','FontSize',16) %,'fontweight', 'bold'); 
ylabel(strPn,'FontSize',16) %,'fontweight', 'bold');
%title(' \lambda _f = 0,  \lambda _c = 1,  \tau = 1, 1e4 Chains','FontSize',16) %, 'fontweight', 'bold')
title(strTitle,'FontSize',16);
%annotation('textbox',[0.2,0.4,0.1,0.1],...
%           'String', str);
legend('Monte Carlo');

%strsave = 

export_fig(sprintf('P%i_sigf%.2f_sigc%.2f_%ichains',N(k),lf,lc,chains), '-pdf', '-m2.5', '-nocrop');

end
