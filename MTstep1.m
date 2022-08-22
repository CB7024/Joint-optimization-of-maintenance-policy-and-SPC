%计算机蒙特卡洛仿真计算情景概率
a=52.3;
a0=138.9;
a1=84.7;
b=1.7;
b0=1.4;
b1=1.4;
Tp=100;%自定义决策参数Tp=100
h=1;
m=0;%初始化参数
q1=0;
q2=0;
q3=0;
q4=0;
q5=0;
ARL1=10;%暂时跟p的联系
for i=1:500000%1000一般>2m，不然最后m无法满足要求
    if m<200000%m是我想获得的满足条件的仿真数据个数
    t1=wblrnd(a,b);%过程异常发生时间
    t2=wblrnd(a0,b0);%受控时设备失效发生时间
    t3=wblrnd(a1,b1);%失控时设备失效发生时间
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