%% Parameters
params.L = 20; %prediction horizon
params.N = 100;%50; %length of data trajectory (remove or set to 0 for automatic choice)
params.p = 0.8; %risk parameter

params.use_first_step_constraint = 1;

params.Q = diag([1,10]);
params.R = 1;

params.d_max = 0.2; %Bound of disturbance: infnorm(d) <= d_max, used to truncate distribution
params.d_sigma = params.d_max; % Variance of distribution (prior truncation)
params.eps_max =  0.001; %Bound of measurement noise: infnorm(eps) <= eps_max, used to truncate distribution
params.eps_sigma = params.eps_max; % Variance of distribution (prior truncation)

%Sampling-based probabilstic tightening, confidence:
params.eps_u = (1-params.p)*1.05;
params.eps_l = (1-params.p)*0.95;
params.beta = 1e-4;


params.noisy_data.flag = 0; %Specify if measurement noise should be added to data
if params.noisy_data.flag
    params.noisy_data.epsd_max = 0.1; %Bound on noise inside disturbance data (only relevant if flag = 1)
    params.noisy_data.epsx_max = 0.1; %Bound on noise inside state data (only relevant if flag = 1)
end

% Plant model
[sys.A, sys.B, sys.Bd, sys.T] = doubleMassOscillator();
params.u_max = 1; %Input constraints
params.x_max = 2*pi; %Max angles (box constraints for x1 and x2)
params.v_max = pi/2; %Max velocities (box constraints for x3 and x4)

sys.A = [1, 0.0075; -0.143, 0.996];
sys.B = [4.798; 0.115];
sys.Bd = eye(2);
sys.T = 1;

sys.n=2;
sys.m=1;

%State constraints
sys.G_X = [eye(sys.n);-eye(sys.n)];
sys.g_X = repmat([2;3],[2 1]);
sys.X = Polyhedron('A', sys.G_X, 'b', sys.g_X);

%Input constraints
params.u_max = 1;
sys.input.G_U = [eye(sys.m);
    -eye(sys.m)];
sys.input.g_U = params.u_max*ones(2*sys.m,1);
sys.input.U = Polyhedron('A', sys.input.G_U, 'b', sys.input.g_U);
sys.input.u_max = params.u_max;

%% Automated Code starts here
sys = sys_setup(params,sys); % specify distributions and constraints from parameters
data = computeHankelData(sys,params); %Get Hankel matrices

%Compute K and P
[controller.K,controller.P] = get_KP_fromData(data,params.Q,params.R);
disp(['Norm deviation for K: ',num2str(norm(controller.K-sys.controller.realK)),...
    ', for P: ',num2str(norm(controller.P-sys.controller.realP))])

%Compute constraints
% M_matrices = compute_M_matrices(sys, data, params, controller); %This still includes use of A,B!
% stoch_constraints = computeStochasticConstraints(sys, params, data, controller, M_matrices);
% constraints = computeConstraints(sys, params, data, M_matrices, controller, stoch_constraints);
ce_matrices = compute_ce_matrices(sys, data, controller); %This still includes use of A,B!
stoch_constraints = computeStochasticConstraints_ce(sys, params, data, controller, ce_matrices);
constraints = computeConstraints_ce(sys, params, data, ce_matrices, controller, stoch_constraints);

save("MPC_init")
%if sys.n==2
%%
figure(77)

hold on
if isfield(constraints,'polyfs')
    constraints.polyfs.plot('Color','b','alpha',0.6)
    %   constraints.polyinit.plot('Color','g','alpha',0.1)
end
%constraints.Xf.plot('Color', 'black','alpha',0.5)

Z1 = Polyhedron('A',sys.G_X,'b',sys.g_X-stoch_constraints.qx(:,2));
Z1.plot('Color','red','alpha',0.8)

%Z2 = Polyhedron('A',sys.G_X,'b',sys.g_X-stoch_constraints.qx(:,3));
%Z2.plot('Color','yellow','alpha',0.5)
sys.X.plot('Color','black','alpha',0.5)

legend('$\hat{\mathcal{Z}}_I$','$\hat{\mathcal{Z}}_1$','$\mathcal{X}$',...
    'Interpreter','latex','fontsize',13,'location','northwest')
set(gca,'fontsize',11)
xlabel('$x_1$','Interpreter','latex','FontSize',13)
ylabel('$x_2$','Interpreter','latex','FontSize',13)
%
%     figure(78)
%     sys.X.plot()
%     hold on
%     constraints.CL.plot('Color','yellow','alpha',1)
%     constraints.CLR.plot('Color','g','alpha',0.96)
%     if isfield(constraints,'polyfs')
%         constraints.polyfs.plot('Color','b','alpha',1)
%         %constraints.polyinit.plot('Color','g','alpha',0.1)
%     end
%     %constraints.Xf.plot('Color', 'yellow','alpha',1)
%
%     %Z1 = Polyhedron('A',sys.G_X,'b',sys.g_X-stoch_constraints.qx(:,2));
%     %Z1.plot('Color','yellow','alpha',0.5)
%
%     %Z2 = Polyhedron('A',sys.G_X,'b',sys.g_X-stoch_constraints.qx(:,3));
%     %Z2.plot('Color','yellow','alpha',0.5)
%     constraints.Xf.plot('Color', 'black','alpha',0.8)
%
%     legend('X','CL','CLinf','First Step Constraint','Terminal Set')
%     xlabel('x1')
%     ylabel('x2')
%     title("FS analysis")
% end
