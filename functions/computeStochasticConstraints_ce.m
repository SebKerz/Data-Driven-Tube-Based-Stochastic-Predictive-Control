function stochastic_constraints = computeStochasticConstraints_ce(sys, params, data, controller, ce_matrices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      PARAMETERS OF SAMPLING TECHNIQUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = sys.m;
n = sys.n;
md = sys.md;

L = params.L;

%p = params.p;
eps_u = params.eps_u;
eps_l = params.eps_l;
beta = params.beta;

%finding minimal Ns such that cond_l <= cond_u, and then find minimal r:

Nmax = 1e7;
for Ns=1:Nmax
    cond_u = eps_u*Ns - sqrt(2*eps_u*Ns*log(1/beta));         %>=r
    cond_l = eps_l*Ns - 1 + sqrt(3*eps_l*Ns*log(2/beta));     %<=r
    if(floor(cond_l) < floor(cond_u))
        for j=1:Ns
            r=j;
            if(r <= cond_u && r>= cond_l)
                disp(['Risk parameter p=',num2str(params.p),' and confidence 1-Beta=',num2str(1-params.beta),...
                    ' means we use ',num2str(Ns),' samples to solve chance constraint optimization.'])
                disp([num2str(r),' sample-based constraints can be violated.']);
                break;
            elseif (r==Ns && (r > cond_u || r < cond_l))
                disp('r could not be found with this set of parameters!');
            end
        end
        break;
    end
    if Ns==Nmax
        disp('Ns could not be found with this set of parameters.');
    end
end

% Get blkdiag expansion of stabilizing state feedback
K = controller.K;
Kcell = repmat({K}, 1, L);
BigK = blkdiag(Kcell{:});


%% Compute state and input constraint tightening parameters
Hu = data.Hu(1:m*L,:);
Hd = data.Hd(1:md*L,:);
Hx0 = data.Hx(1:n,:);
Hx = data.Hx(n+1:end,:);
G_X = sys.G_X;
G_U = sys.input.G_U;
dim_a = size(Hx,2);

%Get Ns samples of disturbance sequences
ds = random(sys.disturbance.d_pdf,[sys.md*L,Ns]);

E = zeros(size(sys.G_X,1),L,Ns);
E_K = zeros(size(sys.input.G_U,1),L,Ns);

u = zeros(m*L,1); %set input to zero
e0 = zeros(n,1); %initial state error

Aeq = [Hu - BigK*[Hx0;Hx(1:end-n,:)]; Hd; Hx0];
Aeq_pinv = pinv(Aeq);

for j=1:Ns %For each sampled disturbance sequence
    d = ds(:,j);

    %compute resulting state error sequence
    beq = [u; d; e0];
    if params.noisy_data.flag
        alpha = quadprog(100*eye(dim_a), zeros(dim_a,1), [], [], Aeq, beq, [], [], [], optimset('Display','off'));
    else
        alpha = Aeq_pinv*beq;
    end
    e = reshape(Hx*alpha,n,L);

    %And multiply by GX or GU for state or input tightening
    E(:,:,j) = G_X*e;
    E_K(:,:,j) = G_U*K*e;
end

%State constraint tightening parameter eta = g_X - qx
%Input constraint tightening parameter mu = g_U - qu
quant = 1-r/Ns; %ratio of samples for which constraints need to hold
qx = [zeros(size(G_X,1),1), quantile(E,quant,3)];
qu = [zeros(size(G_U,1),1), quantile(E_K,quant,3)];
stochastic_constraints.qx = qx;
stochastic_constraints.qu = qu;

%% Comput Terminal Constraint:
%Terminal state constraints:
A_xf = [G_X;
    G_X*ce_matrices.Acl_ce;
    G_U*K];
b_xf = [sys.g_X;
    sys.g_X - qx(:,2);
    sys.input.g_U];
Xf = Polyhedron('A', A_xf, 'b', b_xf);

%positive invariance:
system = LTISystem('A', ce_matrices.Acl_ce);
Xf = system.invariantSet('X', Xf);
G_Xf = Xf.A;
g_Xf = Xf.b;

if(~sys.X.contains(Xf) || ~sys.input.U.contains(K*Xf))
    disp('Xf infeasible! Change feedback')
else
    %Probabilstic tightening for terminal state set
    Ef = zeros(size(G_Xf,1),1,Ns);
    u = zeros(m*L,1);
    e0 = zeros(n,1);

    %find optimal alpha:
    Aeq = [Hu - BigK*[Hx0;Hx(1:end-n,:)]; Hd; Hx0];
    Aeq_pinv = pinv(Aeq);
    ds = random(sys.disturbance.d_pdf,[md*L,Ns]);
    for j=1:Ns
        %draw random disturbance:
        d = ds(:,j);
        %compute resulting state error
        beq = [u; d; e0];
        if params.noisy_data.flag
            alpha = quadprog(100*eye(dim_a), zeros(dim_a,1), [], [], Aeq, beq, [], [], [], optimset('Display','off'));
        else
            alpha = Aeq_pinv*beq;
        end
        e = reshape(Hx*alpha,n,L);

        %store sample:
        Ef(:,:,j) = G_Xf*e(:,end);
    end
    qxf = quantile(Ef,quant,3);
    stochastic_constraints.qxf = qxf;
    stochastic_constraints.G_Xf = G_Xf;
    stochastic_constraints.g_Xf = g_Xf;
end

