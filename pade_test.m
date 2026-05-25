clear 
close all

addpath("chebfun");

x_restricted0=0;
x_restricted1=1;
x_broad0=-1;
x_broad1=2;
n_restricted=1e5;
n_broad=1e5;

x_sample=linspace(x_restricted0,x_restricted1, n_restricted);
x_broad=linspace(x_broad0,x_broad1,n_broad);

syms z;
f_sym= sin(z);
f_num=matlabFunction(f_sym);
f_sample=f_num(x_sample);


[r,pol,res,zer,wj,zj,fj]=aaa(f_sample, x_sample,'tol',0,'mmax',100);

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

Q=sum(vec);
for j=1:l
    vec(j)=vec(j)*fj(j);
end
P=sum(vec);

aaa_P=coeffs(P,z);
aaa_Q=coeffs(Q,z);
normalizing=aaa_Q(length(aaa_Q));
aaa_P=aaa_P/normalizing;
aaa_Q=aaa_Q/normalizing;

N=length(aaa_P)-1;

[pade_P,pade_Q]=multipoint_pade(fj,zj);
m=max([length(pade_P),length(pade_Q),length(aaa_P),length(aaa_Q)]);
pade_P=paddata(pade_P,m,Side='leading')
pade_Q=paddata(pade_Q,m,Side='leading')
aaa_P=paddata(aaa_P,m,Side='leading')
aaa_Q=paddata(aaa_Q,m,Side='leading')

