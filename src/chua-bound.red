% Chua system with cubic nonlinearity
% structure of the Lyapunov candidate function
load_package(redlog);
rlset r;
on ofsfvs;
on rlqevarseltry;
on rlsimpl;
off rational;

% equations of the Chua system
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

% quadratic Lyapunov candidate function
V0:=p33*x3^2+2*p23*x2*x3+2*p13*x1*x3+p22*x2^2+2*p12*x1*x2+p11*x1^2;

% use computed parameter values
V:=sub(p12=0,p13=0,p33=1,p11=1,p22=1449/214,p23=-35/107,V0);

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

% value for gamma
g:=2021.851;

% expression with quantifiers
phi:=ex(q,all({x1,x2,x3},
   q>0 and dV<=-q*(V-g)));
psi:=rlqe(phi);

% print result
write("Result: ",psi);

