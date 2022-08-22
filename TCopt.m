%先运行 ARLCompute.m，计算出ARL1,ARL0
%已经把在本程序里调用ARLCompute.m了,不需要先运行它了
%计算ECT
%tic
syms x y q %x,y,q是时间t的等价变量
cycleN=10;
%决策变量Tp ,n3，H,L,n1, n2可变化
syms n1 n2 n3 n4  %n3 代表h, n4代表Tp
% ECT=zeros(1,cycleN);
% 
H=6; %控制限
L=4;%控制图的警戒限
% n1=2;
% n2=6;
% n4=100;
% n3=5;%抽样间隔时间
%%%%%%%%%%%
fx=(1.7/52.3)*((x/52.3)^(1.7-1))*exp(-(x/52.3)^1.7); %f(t)系统异常时间服从WEIBUL分布，即系统发生了偏移，从受控变为失控
fy=(1.7/52.3)*((y/52.3)^(1.7-1))*exp(-(y/52.3)^1.7);
Fx=1-exp(-(x/52.3)^1.7);% 系统发生异常时间的概率累积分布函数
Ftp=1-exp(-(n4/52.3)^1.7);%当时间到达TP时的失效的累计概率
g0x=(1.4/138.9)*((x/138.9)^(1.4-1))*exp(-(x/138.9)^1.4);%g0(t)，受控时，设备失效时间服从Weibul分布
g0y=(1.4/138.9)*((y/138.9)^(1.4-1))*exp(-(y/138.9)^1.4);
g0q=(1.4/138.9)*((q/138.9)^(1.4-1))*exp(-(q/138.9)^1.4);
g1x=(1.4/84.7)*((x/84.7)^(1.4-1))*exp(-(x/84.7)^1.4);%g1(t)，系统失控时，设备失效时间服从Weibul分布
g1y=(1.4/84.7)*((y/84.7)^(1.4-1))*exp(-(y/84.7)^1.4);
g1q=(1.4/84.7)*((q/84.7)^(1.4-1))*exp(-(q/84.7)^1.4);   
G0x=1-exp(-(x/138.9)^1.4);
G0y=1-exp(-(y/138.9)^1.4);
G0tp=1-exp(-(n4/138.9)^1.4);
G1x=1-exp(-(x/84.7)^1.4);
G1y=1-exp(-(y/84.7)^1.4);
G1tp=1-exp(-(n4/84.7)^1.4);

[ARL0,ARL1]=ARLComputeopt( n1,n2);

 ARL1=real(ARL1);
 ARL1=abs(ARL1);
mm=1-1/ARL1;%论文中的 4.1节的ββ
m1=(y-x)/n3; %时间y 和x之间的采样次数
m2=(n4-x)/n3;
pxy=1-mm^m1; %x,y时间段内检出过程质量失控的概率
pxtp=1-mm^m2;
pp1=0.6612;%落在安全区的概率，利用ARL0.M程序计算出来
pp2=0.3388;%落在非安全警戒区的概率
ps1=int(g0x*(1-Fx),x,0,n4);
%ps1=double(ps1);
ps2=(1-Ftp)*(1-G0tp)*(1-G1tp);
%ps2=double(ps2);
ps3=int(fx*int(g1y*(1-pxy)*int(g0q,q,y,inf),y,x,n4),x,0,n4);
%ps3=double(ps3);
ps4=int(fx*int(pxy*int(g0q*g1q,q,y,inf),y,x,n4),x,0,n4);
%ps4=double(ps4);
ps5=int(fx*(1-pxtp),x,0,n4)*(1-G0tp)*(1-G1tp);
%ps5=double(ps5);
PZ=ps1+ps2+ps3+ps4+ps5;
P1=ps1/PZ;
P2=ps2/PZ;
P3=ps3/PZ;
P4=ps4/PZ;
P5=ps5/PZ;
%n3=1时，可计算出如下概率值
%P1=0.2176;
%P2=0.0110;
%P3=0.0790;
%P4=0.6845;
%P5=0.0080;
Ta=1;
T1=4;
T2=5;
T3=6;
Ts=0.1;
Cq=5;
Cl1=50;
Cl2=300;
Cd=200;
Cma=20;
Cm1=200;
Cm2=1000;
Cm3=1100;
Cms=2;%初始化成本参数
En=n1*pp1+n2*pp2;%样本容量的期望
T1c1=int(x*g0x*(1-Fx),x,0,n4)/P1;
%T1c1=double(T1c1);
Ts1=T1c1+T3+T1c1/(n3*ARL0)*(Ta+Ts);
%Ts1=double(Ts1);
C1=T1c1/n3*En*Cq+T1c1*Cl1+Cm3+T1c1/(n3*ARL0)*(Cma+Cms)+T3*Cd;%情景1的成本C（S1)
Ts2=n4+T1+n4/(n3*ARL0)*(Ta+Ts);
%Ts2=double(Ts2);
C2=n4/n3*En*Cq+n4*Cl1+Cm1+n4/(n3*ARL0)*(Cma+Cms)+T1*Cd;%情景2
T1c3=int(x*fx*(1-G0x)*int(g1y*(1-pxy)/((1-G1x)*P3),y,x,n4),x,0,n4);
%T1c3=double(T1c3);
T0c3=int(y*g1y*(1-G0y)*int(fx*(1-pxy)/((1-G1x)*P3),x,0,y),y,0,n4)-T1c3;
%T0c3=double(T0c3);
Ts3=T1c3+T0c3+T3+T1c3/(n3*ARL0)*(Ta+Ts);
%Ts3=double(Ts3);
C3=(T1c3+T0c3)/n3*En*Cq+T1c3*Cl1+T0c3*Cl2+Cm3+T1c3/(n3*ARL0)*(Cma+Cms)+T3*Cd;%情景3
T1c4=int(x*fx*(1-G0x)/((1-G1x)*P4)*(int(g1y*pxy,y,x,n4)+int(g1y*pxtp,y,n4,inf)),x,0,n4);
%T1c4=double(T1c4);
T0c4=ARL1*n3;
Ts4=T1c4+T0c4+T2+T1c4/(n3*ARL0)*Ta+(T1c4/(n3*ARL0)+1)*Ts;
C4=(T1c4+n3*ARL1)/n3*En*Cq+T1c4*Cl1+n3*ARL1*Cl2+Cm2+T1c4/(n3*ARL0)*Cma+(T1c4/(n3*ARL0)+1)*Cms+T2*Cd;%情景4
T1c5=int(x*fx*(1-G0x)*int(g1y*(1-pxtp)/((1-G1x)*P5),y,n4,inf),x,0,n4);
%T1c5=double(T1c5);
Ts5=n4+T2+T1c5/(n3*ARL0)*(Ta+Ts);
C5=n4/n3*En*Cq+T1c5*Cl1+(n4-T1c5)*Cl2+Cm2+T1c5/(n3*ARL0)*(Cma+Cms)+T2*Cd;%情景5
%ECT(cycleN)=(C1*P1+C2*P2+C3*P3+C4*P4+C5*P5)/(Ts1*P1+Ts2*P2+Ts3*P3+Ts4*P4+Ts5*P5);
ECT=(C1*P1+C2*P2+C3*P3+C4*P4+C5*P5)/(Ts1*P1+Ts2*P2+Ts3*P3+Ts4*P4+Ts5*P5);

%plot(ECT)
%toc  m