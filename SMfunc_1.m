function result=SMfunc_1(CM,CH,node)

     %%%%%% 簇内距离因子 %%%%%% 
     f2=0;%簇内距离因子
     for i=1:length(CM)
         d_min=inf;
         for j=1:length(CH)
             d=sqrt((node(CM(i)).location(1)-node(CH(j)).location(1))^2+((node(CM(i)).location(2)-node(CH(j)).location(2)))^2);
             if d<d_min
                d_min=d; 
             end
         end
         f2=f2+d_min^2;%记录最佳适应度值
     end
    
     %%%%%% 适应度函数 %%%%%%
     F=f2;
     result=F;
     
end



