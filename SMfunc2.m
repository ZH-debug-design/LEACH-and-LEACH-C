function result = SMfunc2(CH,CH_candidate,node,bounder,node_HE)
     for i=1:length(CH)%����Ŷ�
         Next_CH_location(i,1)=node(CH(i)).location(1)+(-bounder)+rand*(2*bounder);
         Next_CH_location(i,2)=node(CH(i)).location(2)+(-bounder)+rand*(2*bounder);
     end
     for i=1:length(CH)%ѡ�������ڵ�
         d_min=inf;
         for j=1:length(CH_candidate)
             d=sqrt((Next_CH_location(i,1)-node(CH_candidate(j)).location(1))^2+(Next_CH_location(i,2)-node(CH_candidate(j)).location(2))^2);
             if d_min>d
                d_min=d;
                NextCH(i)=CH_candidate(j);
             end
         end
     end
     if length(unique(NextCH))~=length(NextCH)%�����ظ��Ĵ�ͷ,�����ظ���ͷ����Ӹ������ڵ��������ѡ��ͷ����
        UniqueCH=unique(NextCH);
        CH_candidate=node_HE(~ismember(node_HE,UniqueCH));
        count=(length(NextCH)-length(UniqueCH));
        RandSet=randperm(length(CH_candidate),count);
        UniqueCH=[UniqueCH CH_candidate(RandSet)];
        NextCH=UniqueCH;
     end
     result=NextCH;%������
end



