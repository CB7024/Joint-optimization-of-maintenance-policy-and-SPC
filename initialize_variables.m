function f = initialize_variables(pop, M, V, min_range, max_range)%f��һ������Ⱥ������ɵľ���
min = min_range;
max = max_range;
K = M + V;%%K���������Ԫ�ظ�����Ϊ�˱��ڼ��㣬���߱�����Ŀ�꺯������һ���γ�һ�����顣  
%���ڽ���ͱ��죬����Ŀ������Ծ��߱�������ѡ��
for i = 1 : pop
    for j = 1 : V
        f(i,j) = round(min(j) + (max(j) - min(j))*rand(1));%f(i j)��ʾ������Ⱥ�е�i�������еĵ�j�����߱�����
                                                    %���д���Ϊÿ����������о��߱�����Լ�����������ȡֵ
    end
    f(i, K) = evaluate_objective(f(i,1:7)); % M��Ŀ�꺯������ V�Ǿ��߱�������
                                                    %Ϊ�˼򻯼��㽫��Ӧ��Ŀ�꺯��ֵ������Ⱦɫ���V + 1 �� K��λ�á�
end
