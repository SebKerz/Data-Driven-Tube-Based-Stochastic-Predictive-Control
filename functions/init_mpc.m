function mpc_init = init_mpc(sys,params, controller, data,x_ref_user)
Qcell = repmat({params.Q}, 1, params.L+1); %Weighting matrix - states
Qcell{params.L+1} = controller.P;  %terminal weights
BigQ = blkdiag(Qcell{:});

%Weighting matrix - inputs:
Rcell = repmat({params.R}, 1, params.L+1);
BigR = blkdiag(Rcell{:});

%Feedback Gain:
Kcell = repmat({controller.K}, 1, params.L+1);
BigK = blkdiag(Kcell{:});
KHx = BigK*data.Hx;

%Pseudo inverse
Pi = pinv([data.Hu-KHx;data.Hd;data.Hx(1:sys.n,:)])*[data.Hu-KHx;data.Hd;data.Hx(1:sys.n,:)];

%Weights H, f of cost function J = 0.5*alpha'*H*alpha + f'*alpha:
%Note: x'*BigQ*x + u'*BigR*u -> alpha'*data.Hx'*BigQ*data.Hx*alpha + alpha'*data.Hu'*BigR*data.Hu*alpha
H = blkdiag([BigQ+BigK'*BigR*BigK,BigK'*BigR;BigR*BigK, BigR],params.lambda*(Pi')*Pi);
H = (H+H')/2;       %for symmetry, bc of uncertainties smaller than 1e-9

mpc_init.H = H;
mpc_init.BigQ = BigQ;
mpc_init.BigR = BigR;
mpc_init.KHx = KHx;
mpc_init.BigK = BigK;

x_ref = [x_ref_user, repmat(x_ref_user(:,1),1,params.L+1)]; %Extend for last prediction horizon
mpc_init.x_ref = reshape(x_ref,[],1); %Reshape for OCP

%Number of MPC iterations:
mpc_init.sim_length = size(x_ref_user,2);
end