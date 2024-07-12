%%
delta = 0.001:0.001:0.01;
d = 4;

N8 = (1/eps)*(exp(1)/(exp(1)-1))*(log(1./delta)+(d-1)); %Online Sampling Lorenzen 2017 (8)

%Lorenzen (12)
N12 = 5*(log(4./delta)+d*log(40/eps))/eps;

%Lorenzen (13)
p = 0.8; % risk parameter
eps = 1-p;
N_0 = 4.1*(log(21.64./delta)+4.39*d*log2(8*exp(1)/eps))/eps;

p = 0.9; % risk parameter
eps = 1-p;
N_1 = 4.1*(log(21.64./delta)+4.39*d*log2(8*exp(1)/eps))/eps;

p = 0.95; % risk parameter
eps = 1-p;
N_2 = 4.1*(log(21.64./delta)+4.39*d*log2(8*exp(1)/eps))/eps;

p = 0.99; % risk parameter
eps = 1-p;
N_3 = 4.1*(log(21.64./delta)+4.39*d*log2(8*exp(1)/eps))/eps;

p = 0.999; % risk parameter
eps = 1-p;
N_4 = 4.1*(log(21.64./delta)+4.39*d*log2(8*exp(1)/eps))/eps;


figure
semilogy(delta,N_0)
hold on
semilogy(delta,N_1)
semilogy(delta,N_2)
semilogy(delta,N_3)
semilogy(delta,N_4)
legend("p is 0.8","p is 0.9","p is 0.95", "p is 0.99", "p is 0.999")
ylim([100,max(N_4)*10])
xlabel("$\delta$", "Interpreter","Latex")
ylabel("Number of samples")

%%
%% Parameters

%finding minimal Ns such that cond_l <= cond_u, and then find minimal r:
params.p=0.95;
params.eps_u = 0.12;%1 - 0.98*params.p;
params.eps_l = 0.08;%1 - 1.02*params.p;
params.beta = 1e-2;

eps_u = params.eps_u;
eps_l = params.eps_l;
beta = params.beta;

Ns = 100:100:200000;

cond_u = eps_u*Ns - sqrt(2*eps_u*Ns*log(1/beta));         %>=r
cond_l = eps_l*Ns - 1 + sqrt(3*eps_l*Ns*log(2/beta));     %<=r

zci = @(v) find(diff(sign(v)));
diffs_b = @(beta) eps_u*Ns - sqrt(2*eps_u*Ns*log(1/beta)) -(eps_l*Ns - 1 + sqrt(3*eps_l*Ns*log(2/beta)));

betas = 0.001:0.001:0.02;
Nsplot = zeros(size(betas));
for k = 1:length(betas)
    Nsplot(k)=Ns(zci(diffs_b(betas(k))));
end
%figure
plot(betas,Nsplot)
hold on
ylabel('required samples')
xlabel('beta = 1-confidence')

%figure
%plot(Ns,cond_u-cond_l)

