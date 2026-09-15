function out=mc_simplified_euler1(p,N,M,batchSize)
% poenostavljena sibka Eulerjeva metoda
h=p.T/N;
xInit = log(p.S0);

nDone=0;

sum_brez=0;
sumsq_brez = 0;

sum_ant=0;
sumsq_ant = 0;

while nDone<M

 m=min(batchSize,M-nDone);

 X=xInit*ones(m,1);
 Xa = X;

 for n=1:N

   xi=2*(rand(m,1)<0.5)-1;
   dW=sqrt(h)*xi;

   [a,b]=coeff_ab(X,p);
   X=X+a*h+b.*dW;

   [aa,ba] = coeff_ab(Xa,p);
   Xa = Xa+aa*h-ba.*dW;

 end

 G=payoff_discounted(X,p);
 Ga = payoff_discounted(Xa,p);

 Gant=0.5*(G+Ga);

 sum_brez=sum_brez+sum(G);
 sumsq_brez = sumsq_brez+sum(G.^2);

 sum_ant=sum_ant+sum(Gant);
 sumsq_ant = sumsq_ant+sum(Gant.^2);

 nDone=nDone+m;

end

out.N=N;
out.delta = h;
out.M=M;

[out.V_brez,out.SE_brez]=mean_se(sum_brez,sumsq_brez,M);
[out.V_ant,out.SE_ant] = mean_se(sum_ant,sumsq_ant,M);

end

function [a,b]=coeff_ab(x,p)

u=x-log(p.S0)-p.x0;

b=p.f0+p.f1*tanh(u);
a = p.r-0.5*b.^2;

end

function G=payoff_discounted(X,p)

S=exp(X);
payoff = max(S-p.K,0);

G=exp(-p.r*p.T)*payoff;

end

function [mu,se]=mean_se(sumx,sumx2,n)

mu=sumx/n;

v=(sumx2-n*mu^2)/(n-1);
v = max(v,0);

se=sqrt(v/n);

end