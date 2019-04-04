function [sys,x0,str,ts] = tdesmc(t,x,u,flag)

switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 2,
    sys=mdlUpdate(t,x,u);
case 3,
    sys=mdlOutputs(t,x,u);
case {4,9}
    sys=[];
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
sizes = simsizes;
sizes.NumOutputs     = 5;
sizes.NumInputs      = 2;
sizes.NumDiscStates = 33;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = zeros(1,33);
str = [];
ts  = [0.001 0];

function sys=mdlUpdate(t,x,u) %【离散情况下的x(n+1)】
persistent ka
if t==0
    ka=0;
end
a=5
M=0.001;
L=0.001;
k=1;   %小则可以消除抖振，但是误差会增大
fn=1;
lamba=0.005;
%FNN参数

alfa=0.05;
xite=0.35;
%FNN参数

bi=[2 2 2 2 2];
% 	ci=0.1*ones(3,5);   
ci=[-3 -1 0 1 3;-3 -1 0 1 3];

qd=u(1);
dqd=15*pi*cos(pi*t-pi/2);
ddqd=-15*pi*pi*sin(pi*t-pi/2);
p=3;
l=5;
fi=2;
the1=1; %调节后误差变化明显
the2=0.1; %调节后误差变化明显
enta=6;
T=4;
ex=0.01;

q=u(2);
dq=(q-x(1))/L;
ddq=(dq-x(2))/L;
% ddq=(u(2)-2*q_xx+q_x)/L^2;

e=q-qd;
de=dq-dqd;
dde=ddq-ddqd
ee=x(33)+e;

% s=de+lamba*e
% slaw=-fn*sign(s)-k*s;  %指数趋近
% s=de+lamba*abs(de)^(l/p) ;%FSMC

% s=e+the1*e^fi+the2*de^(l/p)
s=e+the1*abs(e)^fi*sign(e)+the2*abs(de)^(l/p)*sign(de)%NFTSM滑膜面
H=x(4)-M*x(3);
%FNN
xi=[x(3),x(4)]';
FS1=0;
for j=1:1:5
   u1(j)=exp(-norm(xi(1)-ci(1,j))^2/(2*bi(j)*bi(j)));
%     gs1=-[(xi(1)-(l1-4))/2]^2;
% 	u1(l1)=exp(gs1);
end
for j=1:1:5
%      u2(j)=exp(-norm(x2-ci_1(2,j)^2/(2*bi_1(j)*bi_1(j)))
     gs2=-norm(xi(2)-ci(2,j))^2/(2*bi(j)*bi(j));
     u2(j)=exp(gs2);
%     gs2=-[(xi(2)-(l2-4))/2]^2;
% 	u2(l2)=exp(gs2);
end
% for j=1:1:5
%     u3(j)=exp(-norm(xi(3)-ci(3,j))^2/(2*bi(j)*bi(j)));
% end
% 

for l1=1:1:5
	for l2=1:1:5
%         for l3=1:1:5
            FS2(5*(l1-1)+l2)=u1(l1)*u2(l2);
            FS1=FS1+u1(l1)*u2(l2);
%         end
	end
end
FS=FS2/(FS1+0.001);
kk=0.005;
gama=2;
% % S=gama*s*FS
SS=(1/M)*s*the2*(l/p)*abs(de)^(l/p-1)*FS*kk;
for i=1:1:5*5
    thta(i,1)=x(i+4);
end
unn=-thta'*FS'

   
%RBF
% KD=0.04;KP=0.03
% KD=0.001;KP=0.01  %有干扰时
% s=de+KD*e+KP*ee;
% tu=H+M*(ddqd+KD*de+KP*e); %conventional formulation of the TDC
% KP=0.05;KI=0.0001;
% tu=H+M*(ddqd+KI*ee+KP*e); %iPI of the TDC
% KP=0.02;KI=0.01;KD=0.05;
% tu=H+M*(ddqd+KP*e+KI*ee+KD*de);
% tu=H+M*(ddqd+k*de+slaw);  %指数趋近滑模控制
% ueq=H+M*(ddq+lamba*(l/p)*abs(de)^(l/p-1)*de)
ccc=1
% ueq=-M*1/the2*p/l*(de^(2-l/p)+the1*fi*e^(fi-1)*de^(2-l/p))+H+M*ddqd;
T1=abs(de)^(2-l/p)*sign(de);
ueq=-M*1/the2*p/l*(T1+the1*fi*e^(fi-1)*T1)+H+M*ddqd
% ucor=-M*k*sign(s);
if s< -ex
    sat=-1;
elseif s>ex
        sat=1;
else sat=s/ex;
end

% ucor=-M*enta/T*sat; %where η> 0 and T> 0 are the convergence factor and the boundary layer thickness, respectively
% unn

% tu=ueq+ucor+unn;  %NFTSM
 aa=9;
bb=20;
cc=2;
ka=ka+x(6)*L
% ucor=-M*enta/T*sat; %where η> 0 and T> 0 are the convergence factor and the boundary layer thickness, respectively
%  ucor=-M*(aa)*sat  %无自适应率
 ucor=-M*(ka+bb)*sat; %有自适应率
% ucor=-M*enta/T*sat; %where η> 0 and T> 0 are the convergence factor and the boundary layer thickness, respectively
dka=s*the2*(l/p)*fi*abs(de)^(l/p-1)*sign(s)  %鲁棒性比普通的要好 
% tu=ueq+unn+ucor;  %NFTSM
 tu=ueq+ucor;

%  if t==2
%     tu=tu+3;
% end 
% tu=H+M*(ddqd+KI*ee+KP*e)-ueq-ucor;

% if t==3 %&& t< 3.002
%     tt=0.01*q+0.1*dq^2+0.5*cos(q)
% %     tt=0.1*q+2*dq^2+5*cos(q)
% %     td=td+0.001;
%     kkk=1
% end 
% % tu=tt+tu;
% td1=8*q+10*dq^2+5*cos(q) % t=4s
% % % td1=0.15*q+2*dq^2+0.5*cos(q)  %t=3s
% % td
% if t==4
%     tu=td1+tu;
% end 
% 
% if t>1.5&&t<1.53
%     tu=0.5*tu;
% end 



sys(1)=u(2);
sys(2)=dq;
sys(3)=ddq;
sys(4)=tu;
for i=5:1:29
    sys(i)=SS(i-4);
end
sys(30)=dqd;
sys(31)=ddqd;
sys(32)=s;
sys(33)=ee;


function sys=mdlOutputs(t,x,u)

sys(1)=x(4);
sys(2)=x(2);
sys(3)=x(32);
sys(4)=x(30);
sys(5)=x(31);