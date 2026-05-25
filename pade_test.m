clear 
close all

addpath("chebfun");

x_restricted0=2;
x_restricted1=3;
x_broad0=-100;
x_broad1=100;
n_restricted=1e5;
n_broad=1e6;

x_sample=linspace(x_restricted0,x_restricted1, n_restricted);
x_broad=linspace(x_broad0,x_broad1,n_broad);

syms z;
f_sym= -sqrt(z^2-1)-z;
f_num=matlabFunction(f_sym);
f_sample=f_num(x_sample);


[r,pol,res,zer,zj,fj,wj]=aaa(f_sample, x_sample,'tol',0,'mmax',100); 

l=length(wj);
vec=sym.empty(0,l);
for j=1:l
 term=1;
 for k=1:l
     if k~=j
         term=term*(z-zj(k));
     end
 end
 vec(j)=wj(j)*term;
end

aaa_Q=sum(vec);
for j=1:l
    vec(j)=vec(j)*fj(j);
end
aaa_P=sum(vec);

[aaa_P,aaa_Q]=numden(simplifyFraction(aaa_P/aaa_Q));

cf=coeffs(aaa_Q,z);
normalizing=cf(length(cf));
aaa_P=aaa_P/normalizing;
aaa_Q=aaa_Q/normalizing;

[pade_P,pade_Q]=multipoint_pade(fj,zj);

% ===================================
% ==== Symbolic Polynomials =========
% ===================================

r_aaa=root(aaa_Q);
r_aaa=vpa(r_aaa);
r_aaa=r_aaa(r_aaa>-100 & r_aaa<100);
r_pade=root(pade_Q);
r_pade=vpa(r_pade);
r_pade=r_pade(r_pade>-100 & r_pade<100);


diff=(pade_P/pade_Q)-(aaa_P/aaa_Q);
% for j=1:length(zj)
%     vpa(subs(diff,z,zj(j)))
% end

% =================================
% ===== Symbolic Coefficients =====
% =================================

% pade_P=coeffs(pade_P);
% pade_Q=coeffs(pade_Q);
% aaa_P=coeffs(aaa_P);
% aaa_Q=coeffs(aaa_Q);
% m=max([length(pade_P),length(pade_Q),length(aaa_P),length(aaa_Q)]);
% pade_P=paddata(pade_P,m,Side='leading');
% pade_Q=paddata(pade_Q,m,Side='leading');
% aaa_P=paddata(aaa_P,m,Side='leading');
% aaa_Q=paddata(aaa_Q,m,Side='leading');

% ==========================================
% ==== Numerical Coefficients ==============
% ==========================================

% pade_P=vpa(pade_P);
% pade_Q=vpa(pade_Q);
% aaa_P=vpa(aaa_P);
% aaa_Q=vpa(aaa_Q);

diff=matlabFunction(diff);
plot(x_broad,diff(x_broad)); hold on;

aaa_zeros=zeros(1,length(r_aaa));
pade_zeros=zeros(1,length(r_pade));
plot(r_aaa,aaa_zeros,'o'); hold on;
plot(r_pade,pade_zeros,'o'); hold off;
