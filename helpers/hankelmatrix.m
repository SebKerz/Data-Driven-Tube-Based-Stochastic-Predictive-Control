function H = hankelmatrix(x, L)
%This function computes the Hankel matrix according to Berberich et al.
%"Data-Driven Model Predictive Control with Stability and Robustness
%Guarantees", Section 2.

    n = size(x,1);
    N = size(x,2);
    H = zeros(n*L, N-L+1);
    window = 0:(N-L);
    
    for i = 1:L
        H((n*(i-1)+1):n*i,:) = x(:,window + i);
    end

end
