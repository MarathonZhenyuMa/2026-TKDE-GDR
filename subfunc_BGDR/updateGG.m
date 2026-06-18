function [G,gg,obj_G,changed,converge] = updateGG(V,G) % O(mc)
% G=离散锚点簇指示矩阵 m*c
% Q=G加权矩阵/簇数c/最大迭代次数10
%% Preliminary
[m,c] = size(G);
obj_G = zeros(11,1);           % 每轮流交替更新一轮G后记录一次目标函数值，第一次是还没迭代时的目标函数

gg = sum(G)+eps*ones(1,c);     % tr(GTG) O(mc)
gv = sum(V.*G);                % tr(VTG) O(mc)

qv = zeros(1,c);
for cc = 1:c                        % O(c)
    qv(cc) = gv(cc)./sqrt(gg(cc));  % O(1)
end
obj_G(1) = sum(qv);            % objg O(c)

changed = zeros(10,1);
incre_G = zeros(1,c);
converge = true;
%% Update
for iterg = 1:10               % O(mct) t<10
    converged = true;
    for i = 1:m                           % O(mc)
        vi = V(i,:);
        [~,id0] = find(G(i,:)==1);
        for k = 1:c                       % O(c)
            if k == id0
                incre_G(k) = gv(k)/sqrt(gg(k)+eps) - (gv(k) - vi(k))/sqrt(gg(k)-1+eps);
            else
                incre_G(k) = (gv(k)+vi(k))/sqrt(gg(k)+1+eps) - gv(k)/sqrt(gg(k)+eps);
            end
        end

        [~,id] = max(incre_G);
        %         [~,id] = max(incre_g);     % 该行对应样本更新后的归属类别 1*1
        if id~=id0
            converged = false;               % not converge
            changed(iterg) = changed(iterg)+1; % change record
            G(i,id0) = 0;G(i,id) = 1;
            gg(id0) = gg(id0) - 1;           % id0 from 1 to 0, number -1
            gg(id)  = gg(id) + 1;            % id from 0 to 1, number +1
            gv(id0) = gv(id0) - vi(id0);     % id0 from 1 to 0, update gv
            gv(id)  = gv(id) + vi(id);       % id from 0 to 1, update gv
        end
    end
    if converged                             % m anchors traversal, false continue, true break
        break;
    end

    %% Obj tr(VT*Q)
    for cc = 1:c
        qv(cc) = gv(cc)/sqrt(gg(cc)+eps);
    end
    obj_G(iterg+1) = sum(qv);

    err_obj_g = obj_G(iterg+1)-obj_G(iterg);
    if err_obj_g < 0
        converge = false;
    end
end

end
    
    
                
            
