function [K,P] = get_KP_fromData(data,Q,R)
    
    md = size(data.D(:,1),1); %dimension of disturbance
    n = size(data.X(:,1),1); %dimension of state
    L = data.L;
    uKd = data.U(:,1:L+1);
    xKd = reshape(data.Hx*pinv([data.Hu;data.Hd;data.Hx(1:n,:)])*[uKd';zeros(md*(L+1),1);data.X(:,1)] , n, []);

    [K, P] = feedbackGainLMI(uKd, xKd, Q, R);
    % [Kref,Pref] = dlqr(A,B,Q,R); Kref = -Kref;
end