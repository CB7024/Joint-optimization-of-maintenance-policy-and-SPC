%��������ؿ����������龰����
a=52.3;
a0=138.9;
a1=84.7;
b=1.7;
b0=1.4;
b1=1.4;
Tp=100;%�Զ�����߲���Tp=100
h=1;
m=0;%��ʼ������
q1=0;
q2=0;
q3=0;
q4=0;
q5=0;
ARL1=10;%��ʱ��p����ϵ
for i=1:500000%1000һ��>2m����Ȼ���m�޷�����Ҫ��
    if m<200000%m�������õ����������ķ������ݸ���
    t1=wblrnd(a,b);%�����쳣����ʱ��
    t2=wblrnd(a0,b0);%�ܿ�ʱ�豸ʧЧ����ʱ��
    t3=wblrnd(a1,b1);%ʧ��ʱ�豸ʧЧ����ʱ��
            if t2<t1 && t2<Tp && t3 >t2
                q1=q1+1;
                else if t1>Tp && t2>Tp && t3>Tp
                    q2=q2+1;
                    else if t2>t3 && t3<Tp && t3<t1+ARL1*h && t3>t1
                        q3=q3+1;
                        else if t2>t1+ARL1*h && t1+ARL1*h<t3 && t1+ARL1*h<Tp 
                            q4=q4+1;
                            else if t2>Tp && t3>Tp && t1+ARL1*h>Tp && t1<Tp 
                                q5=q5+1;
                                end
                            end
                        end
                    end
             end
            m=m+1;
        end
end