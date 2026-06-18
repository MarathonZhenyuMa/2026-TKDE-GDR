function [alpha,alpha0] = updateAlpha_ALM2(A,GGG,Niter_ALM,rho,m,n_view,alpha)
% Obj_alpha = [];

alpha0 = alpha;

AA = [];
for v = 1:n_view
    AA = [AA,reshape(A{v},[m*m,1])]; % nm*n_view
end
M2 = GGG;
mm2 = reshape(M2,[m*m,1]); % nm*1
P = AA'*AA;  % n_view*n_view
q = 2*AA'*mm2; % n_view*1


Lambda = zeros(n_view,1);
mu = 0.0001;

for iter2 = 1:Niter_ALM
    % update beta
    beta = alpha+(1/mu)*(Lambda-P'*alpha);

    % update alpha
    r = (1/mu)*(mu*beta-P*beta+q-Lambda);
    alpha = EProjSimplex_new(r);

    % update Lambda mu
    Lambda=Lambda+mu*(alpha-beta);
    mu=rho*mu;

%     obj_alpha = trace(alpha'*P*alpha-alpha'*q);
%     Obj_alpha = [Obj_alpha obj_alpha];
end
% plot(Obj_alpha);
end

