function ref=mc_reference_extrap_euler2(p,N_ref,M_ref,batchSize)
% referencna vrednost z ekstrapolirano Eulerjevo metodo
h=p.T/N_ref;
h2 = 2*h;
Nc=N_ref/2;

xInit=log(p.S0);

nDone=0;
sum_ref = 0;
sumsq_ref=0;

while nDone<M_ref

 m=min(batchSize,M_ref-nDone);

 Xf=xInit*ones(m,1);
 Xc=Xf;
 Xfa = Xf;
 Xca=Xf;

 for k=1:Nc

   dW1=sqrt(h)*randn(m,1);
   dW2 = sqrt(h)*randn(m,1);

   % drobna pot
   Xf=euler_step(Xf,dW1,h,p);
   Xf=euler_step(Xf,dW2,h,p);

   % antiteticna drobna pot
   Xfa=euler_step(Xfa,-dW1,h,p);
   Xfa=euler_step(Xfa,-dW2,h,p);

   % groba pot
   dWc=dW1+dW2;

   Xc=euler_step(Xc,dWc,h2,p);
   Xca = euler_step(Xca,-dWc,h2,p);

 end

 Gf=payoff_discounted(Xf,p);
 Gc = payoff_discounted(Xc,p);

 Gfa=payoff_discounted(Xfa,p);
 Gca=payoff_discounted(Xca,p);

 Gext=2*Gf-Gc;
 Gexta = 2*Gfa-Gca;

 Gref=0.5*(Gext+Gexta);

 sum_ref=sum_ref+sum(Gref);
 sumsq_ref = sumsq_ref+sum(Gref.^2);

 nDone=nDone+m;

end

[V_ref,SE_ref]=mean_se(sum_ref,sumsq_ref,M_ref);

ref.N_ref=N_ref;
ref.delta_ref = h;
ref.M_ref=M_ref;
ref.V_ref = V_ref;
ref.SE_ref=SE_ref;

end


function Xnew=euler_step(X,dW,h,p)

[a,b]=coeff_ab(X,p);

Xnew=X+a*h+b.*dW;

end


function [a,b]=coeff_ab(x,p)

u=x-log(p.S0)-p.x0;

b=p.f0+p.f1*tanh(u);
a = p.r-0.5*b.^2;

end


function G=payoff_discounted(X,p)

S=exp(X);
payoff=max(S-p.K,0);

G=exp(-p.r*p.T)*payoff;

end


function [mu,se]=mean_se(sumx,sumx2,n)

mu=sumx/n;

v=(sumx2-n*mu^2)/(n-1);
v=max(v,0);

se=sqrt(v/n);

end

