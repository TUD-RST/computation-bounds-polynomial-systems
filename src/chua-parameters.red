% Chua system with cubic nonlinearity
% structure of the Lyapunov candidate function
load_package(redlog);
rlset r;
on ofsfvs;
on rlqevarseltry;
on rlsimpl;
off rational;

% parameters
G:=0.7;
C1:=1/9;
C2:=1;
L:=1/7;
a:=-4183/5120;
b:=189/4096;

% equations of the Chua system
f1:=(G/C1)*(x2-x1)-(1/C1)*(a*x1+b*x1^3);
f2:=(G/C2)*(x1-x2)+(1/C2)*x3;
f3:=-x2/L;
f:={f1,f2,f3};
x:={x1,x2,x3};

% quadratic Lyapunov canidate function
V0:=p33*x3^2+2*p23*x2*x3+2*p13*x1*x3+p22*x2^2+2*p12*x1*x2+p11*x1^2;

% simplification according to the structure, normalization p11=0
V:=sub(p12=0,p13=0,p11=1,V0);

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

% numerator
write("numerator:");
dV:=num(dV);

% monomials of degree 2 and 4 in the Lie derivative
a1:=df(dV,x1,4); 
a2:=df(dV,x2,2);
a3:=df(dV,x3,2);
a4:=df(dV,x1,1,x2,1);
a5:=df(dV,x1,1,x3,1);
a6:=df(dV,x2,1,x3,1);

write("conditions:");
bed:=a1<0 and a2<0 and a3<0 and a6=0 and p22>0 and p33>0 and p22-p23^2>0;

% possibility to eliminate some variables
psi:=rlqe(ex({},bed));

% store results
off nat$ 
out "c0.ineq"$
write psi;
shut "c0.ineq"$
on nat;

