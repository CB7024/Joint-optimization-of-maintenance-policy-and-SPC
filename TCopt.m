%������ ARLCompute.m�������ARL1,ARL0
%�Ѿ����ڱ����������ARLCompute.m��,����Ҫ����������
%����ECT
%tic
syms x y q %x,y,q��ʱ��t�ĵȼ۱���
cycleN=10;
%���߱���Tp ,n3��H,L,n1, n2�ɱ仯
syms n1 n2 n3 n4  %n3 ����h, n4����Tp
% ECT=zeros(1,cycleN);
% 
H=6; %������
L=4;%����ͼ�ľ�����
% n1=2;
% n2=6;
% n4=100;
% n3=5;%�������ʱ��
%%%%%%%%%%%
fx=(1.7/52.3)*((x/52.3)^(1.7-1))*exp(-(x/52.3)^1.7); %f(t)ϵͳ�쳣ʱ�����WEIBUL�ֲ�����ϵͳ������ƫ�ƣ����ܿر�Ϊʧ��
fy=(1.7/52.3)*((y/52.3)^(1.7-1))*exp(-(y/52.3)^1.7);
Fx=1-exp(-(x/52.3)^1.7);% ϵͳ�����쳣ʱ��ĸ����ۻ��ֲ�����
Ftp=1-exp(-(n4/52.3)^1.7);%��ʱ�䵽��TPʱ��ʧЧ���ۼƸ���
g0x=(1.4/138.9)*((x/138.9)^(1.4-1))*exp(-(x/138.9)^1.4);%g0(t)���ܿ�ʱ���豸ʧЧʱ�����Weibul�ֲ�
g0y=(1.4/138.9)*((y/138.9)^(1.4-1))*exp(-(y/138.9)^1.4);
g0q=(1.4/138.9)*((q/138.9)^(1.4-1))*exp(-(q/138.9)^1.4);
g1x=(1.4/84.7)*((x/84.7)^(1.4-1))*exp(-(x/84.7)^1.4);%g1(t)��ϵͳʧ��ʱ���豸ʧЧʱ�����Weibul�ֲ�
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
mm=1-1/ARL1;%�����е� 4.1�ڵĦ¦�
m1=(y-x)/n3; %ʱ��y ��x֮��Ĳ�������
m2=(n4-x)/n3;
pxy=1-mm^m1; %x,yʱ����ڼ����������ʧ�صĸ���
pxtp=1-mm^m2;
pp1=0.6612;%���ڰ�ȫ���ĸ��ʣ�����ARL0.M����������
pp2=0.3388;%���ڷǰ�ȫ�������ĸ���
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
%n3=1ʱ���ɼ�������¸���ֵ
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
Cms=2;%��ʼ���ɱ�����
En=n1*pp1+n2*pp2;%��������������
T1c1=int(x*g0x*(1-Fx),x,0,n4)/P1;
%T1c1=double(T1c1);
Ts1=T1c1+T3+T1c1/(n3*ARL0)*(Ta+Ts);
%Ts1=double(Ts1);
C1=T1c1/n3*En*Cq+T1c1*Cl1+Cm3+T1c1/(n3*ARL0)*(Cma+Cms)+T3*Cd;%�龰1�ĳɱ�C��S1)
Ts2=n4+T1+n4/(n3*ARL0)*(Ta+Ts);
%Ts2=double(Ts2);
C2=n4/n3*En*Cq+n4*Cl1+Cm1+n4/(n3*ARL0)*(Cma+Cms)+T1*Cd;%�龰2
T1c3=int(x*fx*(1-G0x)*int(g1y*(1-pxy)/((1-G1x)*P3),y,x,n4),x,0,n4);
%T1c3=double(T1c3);
T0c3=int(y*g1y*(1-G0y)*int(fx*(1-pxy)/((1-G1x)*P3),x,0,y),y,0,n4)-T1c3;
%T0c3=double(T0c3);
Ts3=T1c3+T0c3+T3+T1c3/(n3*ARL0)*(Ta+Ts);
%Ts3=double(Ts3);
C3=(T1c3+T0c3)/n3*En*Cq+T1c3*Cl1+T0c3*Cl2+Cm3+T1c3/(n3*ARL0)*(Cma+Cms)+T3*Cd;%�龰3
T1c4=int(x*fx*(1-G0x)/((1-G1x)*P4)*(int(g1y*pxy,y,x,n4)+int(g1y*pxtp,y,n4,inf)),x,0,n4);
%T1c4=double(T1c4);
T0c4=ARL1*n3;
Ts4=T1c4+T0c4+T2+T1c4/(n3*ARL0)*Ta+(T1c4/(n3*ARL0)+1)*Ts;
C4=(T1c4+n3*ARL1)/n3*En*Cq+T1c4*Cl1+n3*ARL1*Cl2+Cm2+T1c4/(n3*ARL0)*Cma+(T1c4/(n3*ARL0)+1)*Cms+T2*Cd;%�龰4
T1c5=int(x*fx*(1-G0x)*int(g1y*(1-pxtp)/((1-G1x)*P5),y,n4,inf),x,0,n4);
%T1c5=double(T1c5);
Ts5=n4+T2+T1c5/(n3*ARL0)*(Ta+Ts);
C5=n4/n3*En*Cq+T1c5*Cl1+(n4-T1c5)*Cl2+Cm2+T1c5/(n3*ARL0)*(Cma+Cms)+T2*Cd;%�龰5
%ECT(cycleN)=(C1*P1+C2*P2+C3*P3+C4*P4+C5*P5)/(Ts1*P1+Ts2*P2+Ts3*P3+Ts4*P4+Ts5*P5);
ECT=(C1*P1+C2*P2+C3*P3+C4*P4+C5*P5)/(Ts1*P1+Ts2*P2+Ts3*P3+Ts4*P4+Ts5*P5);

%plot(ECT)
%toc  m