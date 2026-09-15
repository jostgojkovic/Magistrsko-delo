function out=mc_simplified_weak2(p,N,M,batchSize)
% poenostavljena sibka Taylorjeva metoda reda 2

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

   dW=simplified_increment_order2(m,h);

   X=step_sw2(X,dW,h,p);
   Xa = step_sw2(Xa,-dW,h,p);

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

function Xnew=step_sw2(X,dW,h,p)

[a,b,ap,bp,app,bpp]=coeffs(X,p);

term1=a*h+b.*dW;

term2=0.5*b.*bp.*(dW.^2-h);

term3=0.5*(ap.*b+a.*bp+0.5*bpp.*b.^2).*dW*h;

term4=0.5*(a.*ap+0.5*app.*b.^2)*h^2;

Xnew=X+term1+term2+term3+term4;

end

function dW=simplified_increment_order2(m,h)

U=rand(m,1);
dW=zeros(m,1);

val=sqrt(3*h);

dW(U<1/6)=val;
dW(U>5/6)=-val;

end

function [a,b,ap,bp,app,bpp]=coeffs(x,p)

u=x-log(p.S0)-p.x0;

th=tanh(u);
sech2 = 1-th.^2;

b=p.f0+p.f1*th;
bp = p.f1*sech2;
bpp=-2*p.f1*th.*sech2;

a=p.r-0.5*b.^2;
ap = -b.*bp;
app=-(bp.^2+b.*bpp);

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