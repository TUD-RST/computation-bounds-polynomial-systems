% Chua system with cubic nonlinearity
% structure of the Lyapunov candidate function
load_package(redlog);
rlset r;
on ofsfvs;
on rlqevarseltry;
on rlsimpl;
off rational;

% condition on the parameters
bedpar:=G>0 and C1>0 and C2>0 and L>0 and a<0 and b>0;

% equations of the Chua system
f1:=(G/C1)*(x2-x1)-(1/C1)*(a*x1+b*x1^3);
f2:=(G/C2)*(x1-x2)+(1/C2)*x3;
f3:=-x2/L;
f:={f1,f2,f3};
x:={x1,x2,x3};

% quadratic Lyapunov canidate function
V:=p33*x3^2+2*p23*x2*x3+2*p13*x1*x3+p22*x2^2+2*p12*x1*x2+p11*x1^2;

% Lie derivative
procedure lieder(f,h,x);
	begin scalar n,i,l;
		n:=length(x);
		l:=for i:=1:n sum df(h,part(x,i))*part(f,i);
	return(l);
	end;

% Lie derivative of the Lyapunov canidate function
write("Lie derivative:");
dV:=lieder(f,V,x);

% remove the constant denominator
write("denominator:");
den(dV);

write("numerator:");
dV:=num(dV);

% conditions for negative definite derivative
procedure conddv(dVi,xi);
	begin scalar u2,u3,u4,phi,psi;
		u2:=coeffn(dVi,xi,2);
		u3:=coeffn(dVi,xi,3);
		u4:=coeffn(dVi,xi,4);
		phi:=all({lambda1,lambda2},bedpar and ((-1<=lambda1 and lambda1<=1 and -1<=lambda2 and lambda2<=1) impl (u4<=0)));
		psi:=rlqe(phi);
		return(psi);
	end;

% direction x1
write("direction x1:");
dV1:=sub(x2=lambda1*x1,x3=lambda2*x1,dV);
d:=deg(dV1,x1);
write("degree in x1: ",d," : ",coeffn(dV1,x1,d));

% direction x2
write("direction x2:");
dV2:=sub(x1=lambda1*x2,x3=lambda2*x2,dV);
d:=deg(dV2,x2);
write("degree in x2: ",d," : ",coeffn(dV2,x2,d));

% directions x3
write("direction x3:");
dV3:=sub(x1=lambda1*x3,x2=lambda2*x3,dV);
d:=deg(dV3,x3);
write("degree in x3: ",d," : ",coeffn(dV3,x3,d));

% generate necessary conditions
psi1:=conddv(dV1,x1);
psi2:=conddv(dV2,x2);
psi3:=conddv(dV3,x3);

% combine conditions
psi:=psi1 and psi2 and psi3;

% store results
off nat$ 
out "c.ineq"$
write psi;
shut "c.ineq"$
on nat;
