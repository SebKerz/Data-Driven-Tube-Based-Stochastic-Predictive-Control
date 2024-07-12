function [A, b, Aeq, beq, lb, ub] = linconstraints(sys, data_hankel, constraints, K, x0)
% This function defines the linear constraint parameters for the 
% optimization variable [x;u;alpha].
%
% A*[x;u;alpha] <= b, Aeq*[x;u;alpha] = beq, lb <= [x;u;alpha] <= ub
%

%Hankel data:
n = sys.n;
m = sys.m;
md = sys.md;
Hu = data_hankel.Hu;
Hd = data_hankel.Hd;
Hx = data_hankel.Hx;
L = data_hankel.L;

Kcell = repmat({K}, 1, L+1);
BigK = blkdiag(Kcell{:});

%Constraint data:
G_X = constraints.G_X;
g_X = constraints.g_X;
G_U = constraints.G_U;
g_U = constraints.g_U;

dim_alpha = size(Hx,2);

%% Constraints

A = blkdiag(G_X, G_U, zeros(m), zeros(dim_alpha));
A(size(G_X,1)+(1:size(G_U,1)), 1:n*(L+1)) = G_U*BigK(1:m*L,:);

b = [g_X;
     g_U;
     zeros(m,1);
     zeros(dim_alpha,1)];

%We dont need to constrain the initial step, because it is fixed anyway
%(given by the measurement). Remove the associated rows from A,b
A(1:size(sys.G_X,1),:)=[];
b(1:size(sys.G_X,1))=[];

Aeq = [eye(n), zeros(n, n*L + m*(L+1) + dim_alpha);
       eye(n*(L+1)), zeros(n*(L+1),m*(L+1)), -Hx;
       zeros(m*(L+1),n*(L+1)), eye(m*(L+1)), -(Hu - BigK*Hx);
       zeros(md*(L+1),n*(L+1)+m*(L+1)), -Hd];
beq = [x0;
       zeros(n*(L+1) + m*(L+1) + md*(L+1), 1)];

lb = [];
ub = [];

end

