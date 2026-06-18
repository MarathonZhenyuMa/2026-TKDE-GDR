function [g_ind, obj] = fast_cd_BB(A,g_ind,m,cd_iter)
% 总复杂度 O(m^2 c + m^2 T)

  if nargin < 4
    cd_iter = 50;
  end

  G = ind2vec(g_ind')';    % 初始标签G
  n_clu = full(sum(G));    % 1*c向量 每簇个数 O(mc)

  H = A * G;               % O(m^2 c)
  D = ones(m,1);           % D = full(sum(H, 2));  
  yAy = full(sum(H .* G)); % O(mc)
  yDy = full(sum(G));      % yDy = full(sum(H)); 
  H = full(H);  % For efficient slicing

  obj(1) = sum(yAy ./ yDy);
  for iter = 1:cd_iter
    for u = 1:m
      p = g_ind(u);
      if n_clu(p) == 1
        continue;
      end

      a_Y = H(u, :);

      yAy_k = yAy + 2 * a_Y;
      yAy_k(p) = yAy(p);

      yDy_k = yDy + D(u);
      yDy_k(p) = yDy(p);

      yAy_0 = yAy;
      yAy_0(p) = yAy(p) - 2 * a_Y(p);

      yDy_0 = yDy;
      yDy_0(p) = yDy(p) - D(u);

      delta = yAy_k ./ yDy_k - yAy_0 ./ yDy_0;

      [~, r] = max(delta);
      if r ~= p
        yAy(p) = yAy_0(p);
        yDy(p) = yDy_0(p);
        yAy(r) = yAy_k(r);
        yDy(r) = yDy_k(r);

        A_m = A(:, u);            % O(m)
        H(:, r) = H(:, r) + A_m;  % O(m)
        H(:, p) = H(:, p) - A_m;  % O(m)

        g_ind(u) = r;
        n_clu(p) = n_clu(p) - 1;
        n_clu(r) = n_clu(r) + 1;
      end
    end
    obj(iter + 1) = sum(yAy ./ yDy);

    if iter > 2 && abs((obj(iter + 1) - obj(iter)) / obj(iter)) < 1e-9
      break;
    end
  end
end
