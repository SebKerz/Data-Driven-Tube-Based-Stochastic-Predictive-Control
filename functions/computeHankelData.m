% This script is used to compute the Hankel matrices and necessary data
% Disturbance and state data are affected by noise
function data = computeHankelData(sys, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of states:
n = sys.n;

%Number of inputs:
m = sys.m;

%Number of disturbances:
md = sys.md;

%Desired prediction horizon / order of persistency of excitation:
L = params.L;
disp(['Prediction horizon set to ',num2str(L)])

if isfield(params,'N') && params.N > 0
        N = params.N;  
        disp(['Length of data trajectory user-specified as ',num2str(N),'. Minimum is ',num2str((m+md+1)*(L+1) + n - 1),'.'])  
else
        N = (m+md+1)*(L+1) + n - 1;
        N = N + 1;
        disp(['Length of data trajectory set to ',num2str(N),', 1 above minimum.'])
end

rng('shuffle'); % random seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      COLLECT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Random input sequence:
u = sys.input.u_max*(2*rand(m,N)-1);

%Random disturbance sequence:
d = random(sys.disturbance.d_pdf,[md,N]);

%Generate states with dynamics (sys.A, sys.B, sys.Bd)
x0 = zeros(n,1); %Initial state
x = zeros(n,N); 
x(:,1)=x0;
for k = 2:N
    x(:,k)=sys.A*x(:,k-1)+sys.B*u(:,k-1)+sys.Bd*d(:,k-1);
end

%Add noise if specified
if params.noisy_data.flag
    eps_x = random(sys.noise.epsx_pdf,[n,N]); %eps_x = epsd_max; %eps_x = 0;
    eps_d = random(sys.noise.epsd_pdf,[md,N]); %eps_d = epsd_max; %eps_d = 0;
    x = x + eps_x;    
    d = d + eps_d;  %no measurement noise on disturbance data

    %Alternative: Faulwasser disturbance data estimation
%     S = [x(:,1:end-1);u(:,1:end-1)];
%     dest = x(:,2:end)*(eye(N-1) - pinv(S)*S);
%     %last entry can be chosen freely:
%     d = [dest,rand(n,1)];
%     %the estimate gives dimension n instead of md:
%     md = n;
end

%Hankel matrices:
Hd = hankelmatrix(d,L+1);
Hx = hankelmatrix(x,L+1);
Hu = hankelmatrix(u,L+1);

%Check persistency of excitation
if rank(Hd) == md*(L+1) && rank(Hu) == m*(L+1)
    disp('Generated data with pers. exciting inputs.')
else
    error('Input or disturbance sequence not persistently exciting')
end

%Put everything into return struct
data.Hx = Hx;
data.Hu = Hu;
data.Hd = Hd;
data.L = params.L;
data.U = u;
data.X = x;
data.D = d;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Reconstruction errors:
if(norm(e_Kref - e_K) > 1e-9)    
    disp(['Reconstruction error - e_K: ', num2str(norm(e_Kref - e_K))])
end

if(norm(x_Kref - x_K) > 1e-9)
    disp(['Reconstruction error - x_K: ', num2str(norm(x_Kref - x_K))])
end

if(norm(z_Kref - z_K) > 1e-9)
    disp(['Reconstruction error - z_K: ', num2str(norm(z_Kref - z_K))])
end

%Plot prediction vs. reference:
figure;
subplot(3,1,1)
plot(x_K')
hold on
plot(x_Kref')
legend('hankel_1','hankel_2','hankel_3','hankel_4','ref_1','ref_2','ref_3','ref_4')
title('state')

subplot(3,1,3)
plot(e_K')
hold on
plot(e_Kref')
legend('hankel_1','hankel_2','hankel_3','hankel_4','ref_1','ref_2','ref_3','ref_4')
title('state error')

subplot(3,1,2)
plot(z_K')
hold on
plot(z_Kref')
legend('hankel_1','hankel_2','hankel_3','hankel_4','ref_1','ref_2','ref_3','ref_4')
title('nominal state')

%Plot difference of prediction and reference:
figure;
subplot(3,1,1)
plot((x_K - x_Kref)')
legend('e_1','e_2','e_3','e_4')
title(['deviation state, max = ', num2str(max(abs(x_K - x_Kref),[],'all'))]);

subplot(3,1,3)
plot((e_K - e_Kref)')
legend('e_1','e_2','e_3','e_4')
title(['deviation state error, max = ', num2str(max(abs(e_K - e_Kref),[],'all'))]);

subplot(3,1,2)
plot((z_K - z_Kref)')
legend('e_1','e_2','e_3','e_4')
title(['deviation nominal state, max = ', num2str(max(abs(z_K - z_Kref),[],'all'))]);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    REPEATED HANKEL PREDICTION + PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function 'DataPreprocessing.m' can also be used to retrieve a Hankel
% prediction with arbitrary input and disturbance signals:
%{

u = 10*sin((1:N)/10);
d = 1*sin((1:N)*2) + eps_d;
x0 = zeros(n,1);

%Stepsize for recursive data-preprocessing: 1 <= step <= L
step = L;

%Nominal and error state considering state-feedback:
[Hz_K, z_K, He_K, e_K] = DataPreprocessing(Hu, Hd, Hx, u, d, x0+eps_x(:,1), K, step);
x_K = z_K + e_K;

%%Reference signals:
%Stabilized state error dynamics:
e_Kref = dynamics(zeros(n,1), zeros(m,N), d-eps_d, K);

%Stabilized nominal state dynamics:
z_Kref = dynamics(x0, u, 0*d, K);

%Stabilized dynamics:
x_Kref = dynamics(x0, u, d-eps_d, K);

%Reconstruction errors:
if(norm(e_Kref - e_K) > 1e-9)    
    disp(['Reconstruction error - e_K: ', num2str(norm(e_Kref - e_K))])
end

if(norm(x_Kref - x_K) > 1e-9)
    disp(['Reconstruction error - x_K: ', num2str(norm(x_Kref - x_K))])
end

if(norm(z_Kref - z_K) > 1e-9)
    disp(['Reconstruction error - z_K: ', num2str(norm(z_Kref - z_K))])
end

%Plot prediction vs. reference:
figure;
subplot(3,1,1)
plot(x_K')
hold on
plot(x_Kref')
legend('hankel_1','hankel_2','hankel_3','hankel_4','ref_1','ref_2','ref_3','ref_4')
title('state')

subplot(3,1,3)
plot(e_K')
hold on
plot(e_Kref')
legend('hankel_1','hankel_2','hankel_3','hankel_4','ref_1','ref_2','ref_3','ref_4')
title('state error')

subplot(3,1,2)
plot(z_K')
hold on
plot(z_Kref')
legend('hankel_1','hankel_2','hankel_3','hankel_4','ref_1','ref_2','ref_3','ref_4')
title('nominal state')

%Plot difference of prediction and reference:
figure;
subplot(3,1,1)
plot((x_K - x_Kref)')
legend('e_1','e_2','e_3','e_4')
title(['deviation state, max = ', num2str(max(abs(x_K - x_Kref),[],'all'))]);

subplot(3,1,3)
plot((e_K - e_Kref)')
legend('e_1','e_2','e_3','e_4')
title(['deviation state error, max = ', num2str(max(abs(e_K - e_Kref),[],'all'))]);

subplot(3,1,2)
plot((z_K - z_Kref)')
legend('e_1','e_2','e_3','e_4')
title(['deviation nominal state, max = ', num2str(max(abs(z_K - z_Kref),[],'all'))]);
%}
