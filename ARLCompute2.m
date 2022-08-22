%马尔可夫链法计算平均运行链长
%计算出ARO和ARL1
%以下是四个待优化变量
syms n1 n2 
H=6; %控制限
L=4;%控制图的警戒限
%%%%%%%%%%%%%%
m=10;%将受控区域[0 H]分成10份
k=0.5;
g=2*H/(2*m+1);
%h=6;
%syms n1 n2 h
p=[];
for i=0:m 
    c=(i*g)^2;
    ncx2=@(x)ncx2pdf(x,3,c);%非中心卡方分布，自由度为3，偏向参数为c
    for j=0:m
        p(i+1,j+1)=quadgk(ncx2,(k+(j-0.5)*g)^2,(k+(j+0.5)*g)^2);%求出非中心卡方分布概率
    end
end
r1=[1 0 0 0 0 0 0 0 0 0 0];
I=eye(m+1);%单位矩阵
NN=[n1 n1 n1 n1 n1 n1 n1 n2 n2 n2 n2];
ARL0=r1*(I-p)^(-1)*NN';
for jx=-m:m+1
    for ix=-m:m+1
        if jx<0
            h1(ix+m+1,jx+m+1)=normcdf(-k+(jx-ix+0.5)*g,1,1)-normcdf(-k+(jx-ix-0.5)*g,1,1);
        elseif jx==0
            h1(ix+m+1,jx+m+1)=normcdf(-k+(jx-ix+0.5)*g,1,1)-normcdf(k+(jx-ix-0.5)*g,1,1);
        elseif jx>0
            h1(ix+m+1,jx+m+1)=normcdf(k+(jx-ix+0.5)*g,1,1)-normcdf(k+(jx-ix-0.5)*g,1,1);
        end
    end
end
for iy=-m:m+1
    c0=(iy*g)^2;
    ncx3=@(x)ncx2pdf(x,3,c0);%非中心卡方分布，自由度为3，偏向参数为c
    for jy=-m:m+1
        p0(iy+m+1,jy+m+1)=quadgk(ncx3,(k+(jy-0.5)*g)^2,(k+(jy+0.5)*g)^2);%求出非中心卡方分布概率
    end
end
Q=kron(h1,p0);
for jx=-m:m+1
    for ix=-m:m+1
        if jx<0
            h0(ix+m+1,jx+m+1)=normcdf(-k+(jx-ix+0.5)*g,0,1)-normcdf(-k+(jx-ix-0.5)*g,0,1);
        elseif jx==0
            h0(ix+m+1,jx+m+1)=normcdf(-k+(jx-ix+0.5)*g,0,1)-normcdf(k+(jx-ix-0.5)*g,0,1);
        elseif jx>0
            h0(ix+m+1,jx+m+1)=normcdf(k+(jx-ix+0.5)*g,0,1)-normcdf(k+(jx-ix-0.5)*g,0,1);
        end
    end
end
Q0=kron(h0,p0);
[x,y]=eig(Q0);
b=x(:,48);
I0=eye(484);
s=ones(484,1);
ARL1=b'*(I0-Q)^(-1)*s*n2;
