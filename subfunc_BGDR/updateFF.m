function [F,ff,obj_F,changed,converge_f] = updateFF(U,F) % O(nc)
%% Preliminary
[n,c] = size(F);
obj_F = zeros(11,1);           

ff = sum(F);                        % O(nc)
uf = sum(U.*F);                     % O(nc)

up = zeros(1,c);
for cc = 1:c                        % O(c)
    up(cc) = uf(cc)./sqrt(ff(cc));  % O(1)
end
obj_F(1) = sum(up);                 % objf

changed = zeros(10,1);
incre_F = zeros(1,c);
converge_f = true;
%% Update
for iterf = 1:10                    % O(nct) t<5
    converged = true;
    for i = 1:n
        ui = U(i,:);
        [~,id0] = find(F(i,:)==1);
        for k = 1:c                          % O(c)
            if k == id0
                incre_F(k) = uf(k)/sqrt(ff(k)+eps) - (uf(k) - ui(k))/sqrt(ff(k)-1+eps);
            else
                incre_F(k) = (uf(k)+ui(k))/sqrt(ff(k)+1+eps) - uf(k)/sqrt(ff(k)+eps);
            end
        end

        [~,id] = max(incre_F);
        if id~=id0                           % O(1)
            converged = false;               
            changed(iterf) = changed(iterf)+1;
            F(i,id0) = 0; F(i,id) = 1;
            ff(id0) = ff(id0) - 1;           % id0 from 1 to 0, number -1
            ff(id)  = ff(id) + 1;            % id from 0 to 1£¬number +1
            uf(id0) = uf(id0) - ui(id0);
            uf(id)  = uf(id) + ui(id);
        end
    end
    if converged
        break;
    end

    for cc = 1:c
        up(cc) = uf(cc)/sqrt(ff(cc)+eps);
    end
    obj_F(iterf+1) = sum(up);

    err_obj_f = obj_F(iterf+1)-obj_F(iterf);
    if err_obj_f < 0
        converge_f = false;
    end
end
end
    
    
                
            
