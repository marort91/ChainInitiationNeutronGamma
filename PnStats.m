%Data Analysis

clc, clear, clf

ProbsN = load('ProbNStats.txt');
ProbsG = load('ProbGStats.txt');

NtrnFiss = load('ntrnfission.txt');
GammaFiss = load('gammafission.txt');

y = load('ntrn.mult');
v = 0:1:7;
nubaract = sum(y.*v)
%y = load('gamma.mult');
%mubaract = sum(y.*v)

nubar = sum(NtrnFiss)/length(NtrnFiss);
mubar = sum(GammaFiss)/length(GammaFiss);

ntrndev = (1/(sqrt(5e4)))*sum((y-nubar).^2)

lens = length(ProbsN);
chains = 5e4;

ti = 0; tf = 20; N = 51;
t = linspace(ti,tf,N);

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
errorbar(t,PnSums,stddev,stddev,'o-')
axis([0 20 0 1])