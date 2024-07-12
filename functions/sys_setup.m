%Define system parameters
function sys = sys_setup(params,sys)

% Real LQR (to check if data driven gives the same)
[K,P] = dlqr(sys.A,sys.B,params.Q,params.R);
sys.controller.realK = -K;
sys.controller.realP=P;

[sys.n, sys.m] = size(sys.B); %number of states and inputs
sys.md = size(sys.Bd,2); %number of disturbances

% %State constraints
% G_X = [eye(sys.n);
%       -eye(sys.n)];
% x12m = params.x_max;
% x34m = params.v_max;
% g_X =  [x12m;x12m;x34m;x34m;x12m;x12m;x34m;x34m]; %State constraints
% sys.X = Polyhedron('A', G_X, 'b', g_X); 
% sys.G_X = G_X;
% sys.g_X = g_X;

% %Input constraints
% u_max = params.u_max;
% G_U = [eye(sys.m);
%       -eye(sys.m)];
% g_U = u_max*ones(2*sys.m,1);
% sys.input.U = Polyhedron('A', G_U, 'b', g_U);
% sys.input.G_U = G_U;
% sys.input.g_U = g_U;
% sys.input.u_max = u_max;


%% Disturbance and noise
%Bound of disturbance: -d_max <= d <= d_max
d_max = params.d_max;
G_d = [eye(sys.md);
      -eye(sys.md)];
g_d = d_max*ones(2*sys.md,1);
sys.disturbance.d_max = d_max;
sys.disturbance.G_d = G_d;
sys.disturbance.g_d = g_d;

%parameters of the pdf (truncated normal distribution):
sigma = params.d_sigma;  %standard deviation
mu = 0;           %mean
d_pdf = makedist('Normal', mu, sigma);
sys.disturbance.d_pdf = truncate(d_pdf, -d_max+mu, d_max+mu);
sys.disturbance.mu = mu;
sys.disturbance.sigma = sigma;

%Bound of measurement noise: infnorm(eps) <= eps_max
eps_max = params.eps_max;
sys.noise.G_eps = [eye(sys.n);
        -eye(sys.n)];
sys.noise.g_eps = eps_max*ones(2*sys.n,1);
sys.noise.eps_max = eps_max;
%Distribution of measurement noise:
sys.noise.eps_pdf = truncate(makedist('normal','mu',0,'sigma',params.eps_sigma),-eps_max,eps_max);

if params.noisy_data.flag
%Noise inside disturbance data:
    sys.noise.epsd_max = params.noisy_data.epsd_max;
    sys.noise.epsd_pdf = truncate(makedist('normal','mu',0,'sigma',sys.noise.epsd_max/3),-sys.noise.epsd_max,sys.noise.epsd_max);
    
    %Noise inside state data:
    sys.noise.epsx_max = params.noisy_data.epsx_max;
    sys.noise.epsx_pdf = truncate(makedist('normal','mu',0,'sigma',sys.noise.epsx_max/3),-sys.noise.epsx_max,sys.noise.epsx_max);
end