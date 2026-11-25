%% Potentzia elektrikoaren kalkuluafunction y=welek(a,b,Sp,Sn,Rp,Rn,Tdif,N,Ap,An,Lp,Ln,j)  n=300;  x=linspace(a,b,n+1);        f=@(x) (x.*(Sp+Sn)*(Tdif)-x.^2*(Lp*Rp/Ap+Ln*Rn/An))*N;    if isa(f,'function_handle')==1      F=f(x);    else       F=f;            end      y=F;
endfunction
