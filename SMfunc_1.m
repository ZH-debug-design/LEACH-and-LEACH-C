function result=SMfunc_1(CM,CH,node)

     %%%%%% ���ھ������� %%%%%% 
     f2=0;%���ھ�������
     for i=1:length(CM)
         d_min=inf;
         for j=1:length(CH)
             d=sqrt((node(CM(i)).location(1)-node(CH(j)).location(1))^2+((node(CM(i)).location(2)-node(CH(j)).location(2)))^2);
             if d<d_min
                d_min=d; 
             end
         end
         f2=f2+d_min^2;%��¼�����Ӧ��ֵ
     end
    
     %%%%%% ��Ӧ�Ⱥ��� %%%%%%
     F=f2;
     result=F;
     
end



