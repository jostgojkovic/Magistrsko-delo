clear
clc
close all
rng(100)

lambda=2;
mu = 1;
X0=1;
T = 1;

N=2^8;
dt = T/N;

dW=sqrt(dt)*randn(1,N);
W = cumsum(dW);
t=dt:dt:T;

% tocna resitev
Xtocna = X0*exp((lambda-0.5*mu^2)*t + mu*W);

figure
plot(0:dt:T,[X0,Xtocna],'m-','LineWidth',1.2)
hold on

% Euler-Maruyamova metoda
R=4;
Dt = R*dt;
L=N/R;

Xem=zeros(1,L);
X = X0;

for j=1:L
    dWskupaj = sum(dW(R*(j-1)+1:R*j));
  X=X+lambda*X*Dt + mu*X*dWskupaj;
    Xem(j)=X;
end

plot(0:Dt:T,[X0,Xem],'g--*','LineWidth',1.2)
hold off

xlabel('t')
ylabel('X')
legend('Točna rešitev','Euler-Maruyamova aproksimacija','Location','northwest')
grid on


% mocna konvergenca
rng(100)

M=5e4;
Dtvals = 2.^(-(5:10));

mocnaNapaka=zeros(size(Dtvals));

for i=1:length(Dtvals)

 Dt=Dtvals(i);
 L = round(T/Dt);

 napake=zeros(M,1);

 for m=1:M

     dW=sqrt(Dt)*randn(1,L);
     WT = sum(dW);

     % tocna resitev v casu T
     XT=X0*exp((lambda-0.5*mu^2)*T+mu*WT);

     % Euler-Maruyamova aproksimacija
     X=X0;

     for j=1:L
       X = X+lambda*X*Dt+mu*X*dW(j);
     end

     napake(m)=abs(XT-X);
 end

 mocnaNapaka(i)=mean(napake);
end


% izpis rezultatov
fprintf('\nMočna konvergenca:\n')
fprintf('Dt napaka\n')

for i=1:length(Dtvals)
 fprintf('%8.6f    %10.6f\n',Dtvals(i),mocnaNapaka(i))
end


% ocena reda mocne konvergence
x = log(Dtvals);
y=log(mocnaNapaka);

p = polyfit(x,y,1);

q=p(1);
C = exp(p(2));

fprintf('\nOcenjeni red močne konvergence: %.4f\n',q)


%sibka konvergenca

rng(200)

lambda = 2;
mu=0.1;
X0=1;
T = 1;

M=5e4;
DtvalsWeak = 2.^(-(5:10));

sibkaNapaka=zeros(size(DtvalsWeak));

% tocna pricakovana vrednost
tocnoPovprecje=X0*exp(lambda*T);

for i=1:length(DtvalsWeak)

 Dt=DtvalsWeak(i);
 L = round(T/Dt);

 X=X0*ones(M,1);

 for j=1:L
    dW=sqrt(Dt)*randn(M,1);
   X = X+lambda*X*Dt+mu*X.*dW;
 end

 povprecje = mean(X);

 sibkaNapaka(i)=abs(tocnoPovprecje-povprecje);
end


% izpis
fprintf('\nŠibka konvergenca:\n')
fprintf(' Dt   napaka\n')

for i=1:length(DtvalsWeak)
  fprintf('%8.6f %10.6f\n',DtvalsWeak(i),sibkaNapaka(i))
end

% ocena reda sibke konvergence

x=log(DtvalsWeak);
y = log(sibkaNapaka);

p=polyfit(x,y,1);

qWeak = p(1);
CWeak=exp(p(2));

fprintf('\nOcenjeni red šibke konvergence: %.4f\n',qWeak)