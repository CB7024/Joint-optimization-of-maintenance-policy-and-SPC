%����ɷ���������ƽ����������
%�����ARO��ARL1
%�������ĸ����Ż�����
syms n1 n2 
H=6; %������
L=4;%����ͼ�ľ�����
%%%%%%%%%%%%%%
m=10;%���ܿ�����[0 H]�ֳ�10��
k=0.5;
g=2*H/(2*m+1);
%h=6;
%syms n1 n2 h
p=[];
for i=0:m 
    c=(i*g)^2;
    ncx2=@(x)ncx2pdf(x,3,c);%�����Ŀ����ֲ������ɶ�Ϊ3��ƫ�����Ϊc
    for j=0:m
        p(i+1,j+1)=quadgk(ncx2,(k+(j-0.5)*g)^2,(k+(j+0.5)*g)^2);%��������Ŀ����ֲ�����
    end
end
r1=[1 0 0 0 0 0 0 0 0 0 0];
I=eye(m+1);%��λ����
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
    ncx3=@(x)ncx2pdf(x,3,c0);%�����Ŀ����ֲ������ɶ�Ϊ3��ƫ�����Ϊc
    for jy=-m:m+1
        p0(iy+m+1,jy+m+1)=quadgk(ncx3,(k+(jy-0.5)*g)^2,(k+(jy+0.5)*g)^2);%��������Ŀ����ֲ�����
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
