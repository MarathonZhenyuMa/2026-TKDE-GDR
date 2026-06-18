function [Z,ZZ,XX,m] = AnchorSelect(X,n,n_view,q)

Z = cell(1,n_view);  % anchor matrix m*dv
m = 2^q;
XX = [];
len = zeros(n_view,1);
for v = 1:n_view
    XX = [XX X{v}];         % 每个视图横向拼接为n*[\sum d^(v)]
    len(v) = size(X{v},2);  % 记录每个视图的拼接长度，即维度d^(v)
end
[~,locAnchor] = hKM(XX',[1:n],q,1); % q为bkhk层数
ZZ = locAnchor';
t1 = 1;
for v=1:n_view
    t2 = t1+len(v)-1;
    Z{v} = ZZ(:,t1:t2);
    t1 = t2+1;
end