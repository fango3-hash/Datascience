%% Spline arrunt eta M.B. erakoak kalkulatzeko

function [r,y] = splinekubik (x,f,opt,x0,xn)
  if isa(f,'function_handle')
    F=f(x);
  else
    F=f;
  end
   plot(x,F,'-o')
   hold on
 if opt==1 && nargin==3
   M=arrunt(F,x);
  
 else
    M=mugalde(F,x,x0,xn);
   
   
   
 end
 h=diff(x);
 yy=diff(F)./h;
 n3=length(h);
 MM=diff(M)./h(1:n3);
 yy=yy(1:n3);

  m=yy-((MM./6).*(h.^2))-(M(1:n3).*h(1:n3)./2);
r=[min(x):0.005:max(x)];
if opt==1
plot(r,horner(r,x,F,m,M,MM),'-g')
grid on
title(' f(x) funtzioari dagokion hurbilketa s(x) arruntarekin');
xlabel('nodoak');
ylabel('hurbilketa');
legend('f(x) funtzioa','p(x) hurbilketa nodoetan')
else
plot(r,horner(r,x,F,m,M,MM),'-g')
grid on
title(' f(x) funtzioari dagokion hurbilketa mugalde baldintzazko s(x)rekin');
xlabel('nodoak');
ylabel('hurbilketa');
legend('f(x) funtzioa','p(x) hurbilketa nodoetan')
end
r=r;
y=horner(r,x,F,m,M,MM);
end

 function y=arrunt(F,x)
 n=length(x);
 h=diff(x);
 n1=length(h);
  hh=h(1:n1-1)+h(2:n1);
 w=2*hh;
FF=diff(F)./h;
n2=length(FF);
d=6*(FF(2:n2)-FF(1:n2-1));
T1=diag(w(1:n-2));
T2=diag(h(1:n-3),1);
T3=diag(h(1:n-3),-1);
T=T1+T2+T3;
M(1)=0;
M(2:n-1)=crout(T,d);
M(n)=0;

 
y=M;
end
function y=mugalde(F,x,x0,xn)
  n=length(x);
  h=diff(x);
 w(1)=2*(h(1));
 hh=h(2:n-1)+h(1:n-2);
 w(2:n-1)=2*hh;
 w(n)=2*h(n-1);
 FF=diff(F)./h;
 n2=length(FF);
 d(1)=6*(FF(1)-x0);
 d(2:n-1)=6*(FF(2:n2)-FF(1:n2-1));
 d(n)=6*(xn-FF(n-1));
T1=diag(w(1:n));
T2=diag(h(1:n-1),1);
T3=diag(h(1:n-1),-1);
T=T1+T2+T3;
M=crout(T,d);

y=M;
end
function y=horner(r,x,F,m,M,MM)
y=zeros(size(r));
n=length(x);
 for i=1:n-1
    ind1=r<= x(i+1) & r>= x(i);
    r1=r(ind1)-x(i);
    r2=(r(ind1)-x(i)).^2;
    r3=(r(ind1)-x(i)).^3;
  
   y(ind1)=F(i)+(m(i)*r1)+((M(i)/2)*r2)+((MM(i))*r3./6);
  end
  y=y;
end
function y= crout(T,d)
  a0=diag(T);
  b=diag(T,-1);
 c0=diag(T,1);
 length(c0);
  n=size(T)-1;
  
  for i=1:n
    c(i)=c0(i)/a0(i);
  end
  a(1)=a0(1)-(b(1)*c0(1));
  for i=1:n
    a(i+1)=a0(i+1)-(b(i)*c0(i));
  end
  d(n+1)=d(n)/a(n);
  for i=2:n
    d(i)=(d(i)-b(i-1)*d(i-1))/a(i);
  end
  for i=n:-1:1
    d(i)=d(i)-c(i)*d(i+1);
  end
  y=d;
 end
