function [B] = LearnB(B,X,Z,k,n,m,n_view)

for v = 1:n_view
    Bv = zeros(n,m);
    Dis = EuDist2(X{v},Z{v},0);
    [~,idx_m] = sort(Dis,2);
    for j = 1:n
        id = idx_m(j,1:(k+1));
        di = Dis(j,id);
        Bv(j,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);  % Eq.(35) in Nie et al. 2016
        %             gamma(j) = (k*di(k+1)-sum(di(1:k)))/2;
    end
    B{v} = Bv;
end