%% Multi-view Clustering via Compact Graph Discrete Reconstruction (CGDR)
function [result,F,G,alpha,t,Obj,converge] = MvC_CGDR(X,label,k,h,init_label,is_normal,Maxiter)
% F: probabilistic sample indicator matrix n*c
% G: discrete anchor indicator matrix m*c
% t: running time
% Obj: objective functions
% converge: 1-converge 0-not converge

if nargin<7
    Maxiter = 30;
end

if nargin<6
    is_normal = 1;
end

if nargin<5
    init_label = 'N2HI';
end

n = size(X{1},1);

if nargin<4
    if n<=1024
        h = floor(log2(n));
    elseif 1024<n && n<=5e3
        h = 10;
    elseif 5e3<n && n<=1e4
        h = 12;
    elseif n>1e4
        h = 13;
    end
end

if nargin<3
    k = 20;
end


n_view = length(X);

for v = 1:n_view
    X{v} = double(X{v});
    X{v} = full(X{v});
end


c = length(unique(label));
B = cell(1,n_view);  % anchor graphs n*m*V
A = cell(1,n_view);  % compact graphs m*m*V

%% Normalization
if is_normal == 1
for v = 1:n_view
    for i = 1:n
        X{v} = full(X{v});
        X{v}(i,:) =(X{v}(i,:) - mean(X{v}(i,:)))/std(X{v}(i,:));
    end
end
end

%% Initialization
tic;
[Z,ZZ,~,m] = AnchorSelect(X,n,n_view,h);  % Select Anchors
B = LearnB(B,X,Z,k,n,m,n_view);          % Construct Anchor Graphs
alpha = ones(n_view,1)./n_view;          % initialize \alpha avergely
WA = zeros(m,m);
WB = zeros(n,m);
for v = 1:n_view
    WB = WB+alpha(v).*B{v};
    A{v} = B{v}'*B{v};
    WA = WA+alpha(v).*A{v};
end

% zero_diag = 0;
% if zero_diag == 1
%     for j = 1:m
%         WA(j,j) = 0;
%     end
% end

switch init_label
    case 'N2HI'
        g = n2hi(WA, c);  % initialize g via N2HI
    case 'KM'
        g = litekmeans(ZZ,c);
    case 'random'
        g = randi([1,c],1,m);g = g';
    case 'SVD+KM'
        [~,~,V] = svds(WB,c);
        g = litekmeans(V,c);
end

G = full(ind2vec(g')');
GGG = G*diag(sum(G))^(-1)*G'; % GGG = G*(G'G)^(-1)*G';
obj_init = norm(WA-GGG,'fro')^2;
objs = [];

converge = 1;
%% Optimization
for iter1 = 1:Maxiter
    % Update G via Fast CD
    [g,~] = fast_cd_BB(WA,g,m);
    G = full(ind2vec(g')');
    GGG = G*diag(sum(G))^(-1)*G';
    % Update alpha via ALM
    [alpha,~] = updateAlpha_ALM2(A,GGG,500,1.05,m,n_view,alpha);
    % Obj
    WA = zeros(m,m);
    for v = 1:n_view
        WA = WA+alpha(v).*A{v};
    end

%     if zero_diag == 1
%         for j = 1:m
%             WA(j,j) = 0;
%         end
%     end

    objs = [objs,norm(WA-GGG,'fro')^2];
    if iter1==1 && objs(1)>obj_init
        converge = 0;
    elseif iter1>1 && objs(iter1)>objs(iter1-1)
        converge = 0;
    end

    if iter1>1
        a = (objs(iter1-1)-objs(iter1))/objs(iter1-1);
    end
    if iter1>2 && a<1e-4
        break;
    end
end

%% Label Propagation from anchors to samples
F = zeros(n,c);
for v = 1:n_view
    Fv = alpha(v)*B{v}*G;
    F = F+Fv;
end
[~,preY] = max(F,[],2);
t = toc;
% result
result = ClusteringMeasure_All_mzy(label,preY);
Obj = [obj_init,objs];

end