%% Multi-view Clustering via Bipartite Graph Discrete Reconstruction (BGDR)
function [result,F,G,alpha,t,Obj,converge] = MvC_BGDR(X,label,k,h,init_label,is_normal,Maxiter)
% F: probabilistic sample indicator matrix n*c
% G: discrete anchor indicator matrix m*c
% t: running time
% Obj: objective functions
% converge: 1-converge 0-not converge

% Initialize f by label propagation from g

if nargin<7
    Maxiter = 30;
end

if nargin<6
    is_normal = 1;
end

if nargin<5
    init_label  = 'N2HI';
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
B = cell(1,n_view);  % anchor graph n*m
A1 = cell(1,n_view); % compact graph m*m
% A2 = cell(1,n_view); % big graph n*n


%% Normalization
% switch data_normal
%     case 'row'
%         for v = 1:n_view
%             for i = 1:n
%                 X{v} = full(X{v});
%                 X{v}(i,:) =(X{v}(i,:) - mean(X{v}(i,:)))/std(X{v}(i,:)); % 行归一化，强调样本之间差异减少
%             end
%         end
%     case 'column'
%         for v = 1:n_view
%             for j = 1:size(X{v},2)
%                 X{v} = full(X{v});
%                 X{v}(:,j) =(X{v}(:,j) - mean(X{v}(:,j)))/std(X{v}(:,j)); % 列归一化，强调特征尺度统一
%             end
%         end
% end

if is_normal == 1
    for v = 1:n_view
        for i = 1:n
            X{v} = full(X{v});
            X{v}(i,:) =(X{v}(i,:) - mean(X{v}(i,:)))/std(X{v}(i,:)); % 行归一化，强调样本之间差异减少
        end
    end
end

tic;

[Z,ZZ,XX,m] = AnchorSelect(X,n,n_view,h); % Select Anchors
B = LearnB(B,X,Z,k,n,m,n_view);           % Construct Anchor Graphs

%% Initialize
alpha = ones(n_view,1)./n_view; % initialize \alpha avergely
WA1 = zeros(m,m); % WA2 = zeros(n,n);
WB = zeros(n,m);
for v = 1:n_view
    A1{v} = B{v}'*B{v};  % m*m  O(nm^2)
%     A2{v} = B{v}*B{v}';  % n*n  O(mn^2)
%     WA2 = WA2+alpha(v).*A2{v};  % n*n O(Vn^2)
    WA1 = WA1+alpha(v).*A1{v};  % m*m O(Vm^2)
    WB = WB+alpha(v).*B{v}; % n*n
end


switch init_label
    case 'N2HI'
        g = n2hi(WA1, c);  % initialize g via N2HI
%         f = n2hi(WA2, c);  % initialize f via N2HI
        F_fuzzy = WB*ind2vec(g')'; % label propagation
        [~,f] = max(F_fuzzy,[],2);
    case 'KM'
        g = litekmeans(ZZ,c);
        f = litekmeans(XX,c);
%         F_fuzzy = WB*ind2vec(g')'; % label propagation
%         [~,f] = max(F_fuzzy,[],2);
    case 'random'
        g = randi([1,c],1,m);g = g';
        f = randi([1,c],1,n);f = f';
    case 'SVD+KM'
        [U1,~,V1] = svds(WB,c);
        f = litekmeans(U1,c);
        g = litekmeans(V1,c);
end

G = full(ind2vec(g')');

F = zeros(n,c);
for i = 1:n
    F(i,f(i)) = 1;
end
% F = full(ind2vec(f')');
FF = F*diag(sum(F))^(-0.5);
GG = G*diag(sum(G))^(-0.5);

obj_init = norm(WB-FF*GG','fro')^2;
objs = [];

converge = 1;
%% Optimization
for iter1 = 1:Maxiter
    % Update F
    V = WB*GG; % n*c
    [F,ff,~,~,~] = updateFF(V,F); % O(nc)
    FF = F*diag(ff)^(-0.5);

    % Update G
    U = WB'*FF;
    [G,gg,~,~,~] = updateGG(U,G);     % O(mct) t<10
    GG = G*diag(gg)^(-0.5);

    % Update alpha
    [alpha,~] = updateAlpha_ALM(B,FF,GG,500,1.05,n,m,n_view,alpha);


    % Compute Obj
    WB = zeros(n,m);
    for v = 1:n_view
        WB = WB+alpha(v).*B{v};
    end
    objs = [objs,norm(WB-FF*GG','fro')^2];
    if iter1==1 && objs(1)>obj_init
        converge = 0;
    elseif iter1>1 && objs(iter1)>objs(iter1-1)
        converge = 0;
    end

    if iter1>1
        a = (objs(iter1-1)-objs(iter1))/objs(iter1-1);
    end
    if iter1>2 && a<1e-3 % set threshold
        break;
    end
end
[~,preY] = max(F,[],2);
t = toc;
% result
result = ClusteringMeasure_All_mzy(label,preY);
% Fig
Obj = [obj_init,objs];    
% x = 0:1:(length(Obj)-1);
% plot(x,Obj);
% xlabel('# Iterations');
% ylabel('Objective Functions');
