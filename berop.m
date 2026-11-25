%% Bero potentziaren kalkuluafunction y=berop(a,b,Sp,Sn,Rp,Rn,Tdif,N,Ap,An,Lp,Ln,Ep,En,Ti,j)  n=300;  x=linspace(a,b,n+1);       f=@(x) (x.*(Sp+Sn)*((Tdif)+Ti)-0.5*x.^2*(Lp*Rp/Ap+Ln*Rn/An)+(Ep*Ap/Lp+En*An/Ln)*(Tdif))*N;    if isa(f,'function_handle')==1      F=f(x);      else      F=f;           end            y=F;
endfunction
