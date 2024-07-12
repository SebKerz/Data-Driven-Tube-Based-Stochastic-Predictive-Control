function [out1] = compute_ce_matrices(sys, data, controller)
%Compute system matrices A, B, Bd, Acl using data & certainty equivalence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      DATA PREPROCESSING ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read necessary variables:
Hx = data.Hx;
Hu = data.Hu;
Hd = data.Hd;

%Controller:
K = controller.K;

%System parameters:
n = sys.n;
m = sys.m;
md = sys.md;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-9;

%Data-driven system identification (De Persis et al. 2020, Theorem 1):
BBdAcl_ce = Hx(n+1:2*n,:)*pinv([Hu(1:m,:)-K*Hx(1:n,:);Hd(1:md,:);Hx(1:n,:)]); BBdAcl_ce(abs(BBdAcl_ce)<tol) = 0;
B_ce = BBdAcl_ce(:,1:m);
Bd_ce = BBdAcl_ce(:,m+(1:md));
Acl_ce = BBdAcl_ce(:,end-n+1:end);
BBdA_ce = Hx(n+1:2*n,:)*pinv([Hu(1:m,:);Hd(1:md,:);Hx(1:n,:)]); BBdA_ce(abs(BBdA_ce)<tol) = 0;
A_ce = BBdA_ce(:,end-n+1:end);

out1.B_ce = B_ce;
out1.Bd_ce = Bd_ce;
out1.A_ce = A_ce;
out1.Acl_ce = Acl_ce;
end