function [out1] = compute_M_matrices(sys, data, params, controller)
%Data pre-processing incl. stabilizing statefeedback
%   input: sys-struct, hankel data struct
%   If sys-struct does not contain K already, then it is computed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      DATA PREPROCESSING ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read necessary variables:
A = sys.A;
B = sys.B;
Bd = sys.Bd;
Hx = data.Hx;
Hu = data.Hu;
Hd = data.Hd;
L = params.L;

%Blockdiagonal expansion:
K = controller.K;
Kcell = repmat({K}, 1, L+1);
BigK = blkdiag(Kcell{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 1e-9;

%Data-driven system identification (De Persis et al. 2020, Theorem 1):
% BAsi = Hx(n+1:2*n,:)*pinv([Hu(1:m,:)-K*Hx(1:n,:);Hd(1:md,:);Hx(1:n,:)]); BAsi(abs(BAsi)<tol) = 0;
% Bsi = BAsi(:,1:m);
% Asi = BAsi(:,end-n+1:end);
% Bdsi = BAsi(:,m+(1:md));

Acl = A + B*K;
n = sys.n;
m = sys.m;
Mzsi = zeros(n*L,m*L+n);
for i = 1:L
    if(i==1)
        Mzsi(1:n,1:m) = B;
    else
        Mzsi((n*(i-1)+1):n*i, 1:m) = Acl^(i-1)*B;
        Mzsi((n*(i-1)+1):n*i, (m+1):end-n) =  Mzsi((n*(i-2)+1):n*(i-1), 1:m*(L-1));
    end
    Mzsi((n*(i-1)+1):n*i, end-n+1:end) = Acl^i;  
end
%Extend by rows for x0:
Mzsi = [zeros(n,m*L),eye(n);Mzsi];
Mu = Mzsi(:,1:end-n);
Mz0 = Mzsi(:,end-n+1:end);

%Prediction matrix: e_K = Md*d + Me0*e0 (state-feedback dynamics)
Mesi = zeros(n*L,m*L+n);
for i = 1:L
    if(i==1)
        Mesi(1:n,1:m) = Bd;
    else
        Mesi((n*(i-1)+1):n*i, 1:m) = Acl^(i-1)*Bd;
        Mesi((n*(i-1)+1):n*i, (m+1):end-n) =  Mesi((n*(i-2)+1):n*(i-1), 1:m*(L-1));
    end
    Mesi((n*(i-1)+1):n*i, end-n+1:end) = Acl^i;  
end
%Extend by rows for x0:
Mesi = [zeros(n,m*L),eye(n);Mesi];
Md = Mesi(:,1:end-n);
Me0 = Mesi(:,end-n+1:end);

out1.Mu=Mu;
out1.Md=Md;
out1.Me0=Me0;
out1.Mz0=Mz0;
end