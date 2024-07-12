function constraints = computeConstraints_ce(sys, params, data, ce_matrices, controller, stoch_constraints)
A_ce = ce_matrices.A_ce;
Acl_ce = ce_matrices.Acl_ce;
B_ce = ce_matrices.B_ce;
Bd_ce = ce_matrices.Bd_ce;
L = params.L;
n = sys.n;
m = sys.m;
md = sys.md;
Hu = data.Hu;
Hd = data.Hd;
Hx = data.Hx;
Hudx = [Hu;Hd;Hx(1:n,:)];
Pi = eye(size(Hudx,2)) - pinv(Hudx)*Hudx;
K = controller.K;
Kcell = repmat({K}, 1, L+1);
BigK = blkdiag(Kcell{:});
Hudx_cl = [Hu- BigK*Hx;Hd;Hx(1:n,:)];
Pi_cl = eye(size(Hudx,2)) - pinv(Hudx_cl)*Hudx_cl;

rho = zeros(L+1,1);
for i = 1:L+1
    %rho(i) can be used instead of A^i as an overapproximation
    Ai = A_ce^i;
    rho(i) = norm(Ai, inf);
end
%rho_inf = max(rho);

%Set for alpha_eps:
G_eps = sys.noise.G_eps;
g_eps = sys.noise.g_eps;
Eps = Polyhedron('A', G_eps, 'b', g_eps);
% Alpha_eps = Polyhedron('A',G_eps*Hx(1:n,:),'b',g_eps,...
%     'Ae', [Hu;Hd;Pi],'be',zeros((m+md)*(L+1) + size(Pi,1),1));


%Set for alpha_d:
G_d = sys.disturbance.G_d;
g_d = sys.disturbance.g_d;
D = Polyhedron('A', G_d, 'b', g_d);
% G_dcell = cell(L+1,1);
% Big_g_d = [];
% for i=1:L+1
%     G_dcell{i} = G_d;
%     Big_g_d = [Big_g_d; g_d];
% end
% Big_G_d = blkdiag(G_dcell{:});
% Alpha_d = Polyhedron('A',Big_G_d*Hd,'b',Big_g_d,...
%     'Ae', [Hu;Hx(1:n,:);Pi_cl],'be',zeros((md)*(L+1) + n + size(Pi,1),1));
% 

%%  TIGHTENED CONSTRAINT SETS DURING PREDICTION - Data-driven
U = cell(L+1,1);
X = cell(L+1,1);
G_U = sys.input.G_U;
g_U = sys.input.g_U;
G_X = sys.G_X;
g_X = sys.g_X;
U{1} = Polyhedron('A', G_U, 'b', g_U);
X{1} = Polyhedron('A', G_X, 'b', g_X);% - Eps;
for i=2:L+1
    Eps_pred = A_ce^(i-1)*Eps;
    U{i} = Polyhedron('A', G_U, 'b', g_U - stoch_constraints.qu(:,i));
    X{i} = Polyhedron('A', G_X, 'b', g_X - stoch_constraints.qx(:,i)) - Eps_pred;
    X{i} = X{i}.minHRep;
end
Xf = Polyhedron('A', stoch_constraints.G_Xf, 'b', stoch_constraints.g_Xf - stoch_constraints.qxf)- Eps_pred;
Xf = intersect(Xf,X{end});
Xf = Xf.minHRep;
if Xf.isEmptySet()
    disp("Terminal set is empty, this won't work!")
else
    disp("Terminal set okay!")
end

%%  COMPUTE First Step constraint
% First get reachable set - BACKWARDS RECURSION

disp(['Reachable Set - Step ', num2str(L)])

%inverse of A:
Ainv = A_ce\eye(size(A_ce));

%Reachable set:
CL = Ainv*(intersect(Xf,X{end}) + (-B_ce)*U{L});
CL = intersect(CL,X{L});

dec = 3;
%Simplification:
[G_CL, g_CL] = noredund(round(CL.A,dec),round(CL.b,dec));
CL = Polyhedron('A',G_CL,'b',g_CL);

for i = L-1:-1:1
    disp(['Reachable Set - Step ', num2str(i)])
    %Reachable set:
    CL = Ainv*(CL + (-B_ce)*U{i});
    if i>1
        CL = intersect(CL, X{i});
    end

    %Simplification:
    %CL=CL.minHRep();
    [G_CL, g_CL] = noredund(round(CL.A,dec),round(CL.b,dec));
    CL = Polyhedron('A',G_CL,'b',g_CL);
end
disp(['Model-based reachable set volume = ',num2str(CL.volume)])

% %Alternative: Reachable Set - PROJECTION (data-driven)
% %Stack constraint sets:
% G_Xcell = cell(L+1,1);
% G_Ucell = cell(L,1);
% Big_g_X = [];
% Big_g_U = [];
% 
% for i=1:L
%     G_Xcell{i} = X{i}.A;
%     Big_g_X = [Big_g_X; X{i}.b];
%     G_Ucell{i} = U{i}.A;
%     Big_g_U = [Big_g_U; U{i}.b];
% end
% G_Xcell{L+1} = Xf.A;
% Big_g_X = [Big_g_X; Xf.b];
% G_Ucell{L+1} = U{L+1}.A;
% Big_g_U = [Big_g_U; U{L+1}.b];
% 
% Big_G_X = blkdiag(G_Xcell{:});
% Big_G_U = blkdiag(G_Ucell{:});
% 
% %Constraint set for alpha:
% dec = 2;
% G_alpha = [Big_G_X*Hx;
%     Big_G_U*Hu];
% g_alpha = [Big_g_X;
%     Big_g_U];
% 
% 
% %[G_alpha, g_alpha] = noredund(round(G_alpha,dec),round(g_alpha,dec));
% %Noredund uses convhulln, this is hopeless for alpha dimension > 8
% 
% if params.noisy_data.flag
%     Alpha_CL = Polyhedron('A',G_alpha,'b',g_alpha,...
%         'Ae', [Hd;Pi],'be',zeros(md*(L+1) + size(Pi,1),1));
% else
%     Alpha_CL = Polyhedron('A',G_alpha,'b',g_alpha,...
%         'Ae', Hd,'be',zeros(md*(L+1),1));
% end
% 
% %Reachable set:
% %CL = Hx(1:n,:)*Alpha_CL; %Will be impossible to compute
% 
% %Simplification:
% %[G_CL, g_CL] = noredund(round(CL.A,dec),round(CL.b,dec));
% %CL = Polyhedron('A',G_CL,'b',g_CL);
% %disp(['DD Projection reachable set volume = ',num2str(CL.volume)])
% 
% %Compute with small steps
% Lsmall = 3;
% Nsmall = (sys.m+sys.md+1)*(Lsmall+1) + sys.n;
% Hxsmall = hankelmatrix(data.X(:,1:Nsmall),Lsmall);
% Husmall = hankelmatrix(data.U(:,1:Nsmall),Lsmall);
% Hdsmall = hankelmatrix(data.U(:,1:Nsmall),Lsmall);
% 
% %and so on


% COMPUTE RPI SUBSET OF REACHABLE SET
imax = 20; %maximal number of iteration
Dist = Bd_ce*D + A_ce*Eps + Eps; %Extended disturbance Polyhedron

[G_Dist, g_Dist] = noredund(round(Dist.A,2),round(Dist.b,2));
Dist = Polyhedron('A',G_Dist,'b',g_Dist).minVRep; %Simplify

%usys = ULTISystem('A',sys.A,'B',sys.B,'E', eye(sys.n)); %Autonomous disturbed system

%InvSet = usys.invariantSet('X',CL,'U',U{1},'D',Dist);

%RPI subset:
s = mptopt;
s.rel_tol = 1e-6;

constraints.CL = CL;

CLR = CL;
for i = 1:imax
    disp(['RPI Reachable Set - Iteration ', num2str(i-1), ', Volume: ', num2str(CLR.volume)])
    CLnew = Ainv*(CLR - Dist + (-B_ce)*U{1}); %One step backwards in the dynamics
    CLnew = intersect(CLnew,CL); %Intersection with Reachable set CL

    %Simplification:
    [G_CLnew, g_CLnew] = noredund(round(CLnew.A,dec),round(CLnew.b,dec));
    CLnew = Polyhedron('A',G_CLnew,'b',g_CLnew);

    if CLnew.isEmptySet
        disp("Set empty, stopping.")
        break
    end
    %Setdiff1 = CLnew\CLR;
    %Setdiff2 = CLR\CLnew;
    if CLR==CLnew %sum([sum(Setdiff1.volume),sum(Setdiff2.volume)])<1e-4 
        disp(['RPI Reachable Set - Iteration ', num2str(i), ', Volume: ', num2str(CLR.volume)])
        disp('RPI computation converged!')
        s.rel_tol = 1e-6;
        break;
    else
        CLR = CLnew;
    end

    if(i >= imax)
        disp(['RPI Reachable Set - Iteration ', num2str(i), ', Volume: ', num2str(CLR.volume)])
        disp('RPI computation did not converge! Increase maximal iteration number..')
    elseif(Dist.contains(CLR))
        disp('RPI computation infeasible! Smaller disturbance bounds needed..');
        break;
    end
end

constraints.CLR = CLR;

%% Combine all constraints:
%First step constraint:
if params.use_first_step_constraint
    CLRfs = CLR - Dist;
    [G_CLRfs, g_CLRfs] = noredund(round(CLRfs.A,dec),round(CLRfs.b,dec));
end

G_X = cell(L+1,1);
g_X = cell(L+1,1);
G_U = cell(L,1);
g_U = cell(L,1);

for i=1:L+1
    if i == 2 && params.use_first_step_constraint %first step constraint
        G_X{i} = [X{i}.minHRep.A;
            G_CLRfs;
            ];
        g_X{i} = [X{i}.minHRep.b;
            g_CLRfs;
            ];
        G_U{i} = U{i}.minHRep.A;
        g_U{i} = U{i}.minHRep.b;
    elseif i == L+1 %terminal constraint
        G_X{i} = Xf.minHRep.A;
        g_X{i} = Xf.minHRep.b;
    else
        G_X{i} = X{i}.minHRep.A;
        g_X{i} = X{i}.minHRep.b;
        G_U{i} = U{i}.minHRep.A;
        g_U{i} = U{i}.minHRep.b;
    end
end

constraints.G_X = blkdiag(G_X{:});
constraints.g_X = sum(blkdiag(g_X{:}),2);
constraints.G_U = blkdiag(G_U{:});
constraints.g_U = sum(blkdiag(g_U{:}),2);
constraints.Xf=Xf;
if params.use_first_step_constraint
    constraints.polyfs = Polyhedron(G_CLRfs,g_CLRfs);
    %BU = (-sys.B)*sys.input.U;
    %minusInput = constraints.polyfs+(BU);
    constraints.polyinit = inv(A_ce)*(constraints.polyfs+((-B_ce)*sys.input.U));
    if constraints.polyfs.isEmptySet()
        disp("First step constraint set empty, this won't work!")
    else
        if constraints.polyinit.isEmptySet()
            disp("Set of feasible initial states, from which first step constraint set is reachable, is empty. This won't work!")
        else
            disp("First step constraint okay and reachable!")
        end
    end
end

