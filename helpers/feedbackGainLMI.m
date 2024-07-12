function [K, P] = feedbackGainLMI(u, z, Q, R)
%This function computes a stabilizing state feedback gain according to De
%Persis et al. "Formulas for Data-Driven Control: Stabilization,
%Optimality, and Robustness" (2020), Section IV-A.
% INPUTS:
% - u, z:   input and nominal state data
%
% OUTPUTS:
% - K:      resulting stabilizing state feedback
% - P:      P-matrix of the discrete algebraic ricatti equation
%
% Author: J. Teutsch

% if (any(eig(P)<=0) || ~issymmetric(P))
%     error("The matrix P has to be positive definite and symmetric!")
% end

%Data:
U = u(:,1:end-1);
X = z(:,1:end-1);
Xp = z(:,2:end);
tol = 1e-6;

[n, N] = size(z);
m = size(u,1);

%for de-noising:
S = [X;U];
Pi = eye(size(S,2)) - pinv(S)*S;

% %%OPTION 1: data-driven LQR with certainty equivalence:
% cvx_begin sdp quiet
%     
%     variable W(N-1,n) 
%     variable V(m,m) symmetric
% 
%     minimize trace(Q*X*W) + trace(R*V)
%     subject to
%         0.5*[V, U*W; W'*U', X*W] + 0.5*[V, U*W; W'*U', X*W]'>= 0;
%         0.5*[X*W - eye(n), Xp*W; W'*Xp', X*W] + 0.5*[X*W - eye(n), Xp*W; W'*Xp', X*W]' >= 0;
%         Pi*W == 0;
%         diag(V) >= tol;
%         diag(X*W) >= tol;
%         
% cvx_end

%%OPTION 2: data-driven LQR with regularization:
%regularization parameters:
rho = 0;
lambda = 100;
cvx_begin sdp quiet
    
    variable W(N-1,n) 
    variable V(m,m) symmetric
    variable V2(N-1,N-1) symmetric

    minimize trace(Q*X*W) + trace(R*V) + rho*trace(V2) + lambda*norm(Pi*W)
    subject to
        0.5*[V, U*W; W'*U', X*W] + 0.5*[V, U*W; W'*U', X*W]' >= 0;
        0.5*[X*W - eye(n), Xp*W; W'*Xp', X*W] +  0.5*[X*W - eye(n), Xp*W; W'*Xp', X*W]' >= 0;
        0.5*[V2, W; W', X*W] + 0.5*[V2, W; W', X*W]' >= 0;
        diag(V) >= tol;
        diag(V2) >= tol;
        diag(X*W) >= tol;
        
cvx_end

%State feedback gain:
K = U*W/(X*W);
Qp = Q + K'*R*K;
Pw = Qp*X*W;

%find P for given K such that Ak'*P*Ak - P + Q + K'*R*K = 0

% cvx_begin sdp
%     
%     variable Y(N-1,n)
%     
%     minimize norm(Y,1)
%     subject to
%         X*Y == eye(n);
%         U*Y == K;
% 
% cvx_end
Y = pinv([U;X])*[K;eye(n)];

% cvx_begin sdp quiet
%     
%     variable P(n,n) symmetric
%     minimize trace(P)
%     subject to
%         [P - Qp, Y'*Xp'*P; P*Xp*Y, P] >= 0;
%         P >= Qp;
%         
% cvx_end
% 
P = dlyap((Xp*Y)',Qp);

if(any(isnan(K) | isnan(P)))
    error("No solution for K or P found. Try again!")
end

end


