function [num, den] = multipoint_pade(f_sample,x_sample)
    syms x;
    n = length(x_sample);
    g = zeros(n, n);
    g(1,:) = f_sample;

    for j=2:n
        for i=j:n %only the diagonals are relevant.
            g(j,i)=(g(j-1,j-1)-g(j-1,i))/((x_sample(i)-x_sample(j-1))*g(j-1,i));
        end
    end

    a = zeros(0,n);
    for j=1:n
        a(j)=g(j,j);
    end

    f=sym(1);
    for i=1:n-1
        f=a(n-i+1)*(x-x_sample(n-i))/f;
        f=f+1;
    end
    f=a(1)/f;
    [P,Q]=numden(simplifyFraction(f));
    P=coeffs(P,'all');
    Q=coeffs(Q,'all');
    normalizing=Q(length(Q));
    num=P/normalizing;
    den=Q/normalizing;
end