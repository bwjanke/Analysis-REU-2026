function [pcoef, qcoef, p, q] = aaaPolyPQ(zj, fj, wj)
    % uses convention r(z) = q(z) / p(z)

    zj = zj(:);
    wj = wj(:);
    fj = fj(:);
    m = length(zj);

    ellcoef = poly(zj);         % = coefficients of polynomial \prod_{t_j \in T_n} (z - t_j)

    pcoef = zeros(1, m);        % define vectors for coefficients of p, q,
    qcoef = zeros(1, m);        % 1 to m bc degree m-1 poly has m coefficients.

    for j = 1:m

        Ljcoef = deconv(ellcoef, [1, -zj(j)]);      % compute L_j; refer to preprint

        pcoef = pcoef + wj(j) * Ljcoef;
        qcoef = qcoef + wj(j) * Ljcoef * fj(j);     %compute individual terms in p, q

    end

    p = @(z) polyval(pcoef, z);
    q = @(z) polyval(qcoef, z);

end