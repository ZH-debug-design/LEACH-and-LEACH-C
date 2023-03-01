% close all;
clear
% clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1.一阶无线电能耗模型（First Order Radio Model）%%%%%%%%%%%%%%%%%%%
for zh=1:1
    N=100;%节点个数
    sinkx=50;sinky=175;%基站位置
    xm=100;ym=100;%监测范围
    p=0.05;%簇头节点比例
    Eo=0.1;%初始能量
    Eelec=50*10^(-9);%电路能耗（Joules/bit）
    ETx=50*10^(-9);%传送能耗（Joules/bit）
    ERx=50*10^(-9);%接收能耗（Joules/bit）
    Efs=10*10^(-12);%自由空间模型能耗系数
    Emp=0.0013*10^(-12);%多径衰落模型能耗系数
    d0=sqrt(Efs/Emp);%临界值
    EDA=5*10^(-9);%数据聚合能耗
    Lc=100;%控制包大小（bit）
    Ld=4000;%数据包大小（bit）
    AD_Range=sqrt(xm^2+ym^2);%广播范围
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2.网络初始化 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for zh=1:1
    load location.mat;%载入位置
    for i=1:N
        node(i).energy=Eo;%剩余能量
%         node(i).location=[rand*xm rand*ym];%节点位置
        node(i).location=location(i,:);%节点位置
        node(i).status=0;%表明节点是否成为或成为过簇头：0为簇成员节点，1为簇头节点
        node(i).CH=0;% 表明节点的簇头
        node(i).CM=[];% 表明节点的簇成员
        node(i).CM_count=0;% 表明簇内的簇成员节点数
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3.网络运行  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for zh=1:1
    
    %%%%%%性能指标参数设置（包括存活节点数量、网络剩余能量、控制开销）%%%%%%
    for z=1:1
        r=0;%轮数
        number_AliveNodes=N;%存活节点数量
        residual_Energy=N*Eo;%网络剩余能量
        contorl_Overhead=0;%控制开销
        AliveNodes=1:1:N;%存活节点
        Throughput=0;%吞吐量
        %         figure(1)%画图
        %         for i=1:N
        %             plot(node(i).location(1),node(i).location(2),'bo','linewidth',1);
        %             hold on;
        %         end
        %         plot(sinkx,sinky,'rd','linewidth',2);hold on;
    end
    
    while number_AliveNodes>0%当节点存活数量>0时，网络正常运行
        %%%%%%%%%%%%%%%%%%%%%%%%%%% ①簇建立阶段 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%节点判断自己是否成为簇头，节点选择簇头加入%%%%%%%%%
        for z=1:1
            %%%%判断是否所有节点都已成为簇头%%%%
            if mod(r,round(1/p))==0%若轮换簇头完成一轮，则更新节点的簇头状态
                for i=1:N
                    node(i).status=0;
                end
            end
            
            %%%%节点判断自己是否为簇头%%%%
            CH=[];%本轮选择的簇头集合
            count_CH=0;%簇头个数
            for i=1:N
                if node(i).energy>0 && node(i).status==0 %在非死亡及没有成为过簇头的节点中选择簇头
                    if rand<=(p/(1-p*mod(r,round(1/p))))
                        node(i).status=1;%更改节点状态
                        CH=[CH i];%将节点加入簇头集
                        count_CH=count_CH+1;%统计簇头数量
                    end
                end
            end
            CM=AliveNodes(~ismember(AliveNodes,CH));%簇成员集合
            
            %%%%%%%%% 画图（簇头与簇成员节点连线）%%%%%%%%%
            CM=AliveNodes(~ismember(AliveNodes,CH));%簇成员集合
            %             figure(2)
            %             if count_CH>0
            %                 for i=1:length(CM)
            %                      d_min=inf;
            %                      for j=1:length(CH)
            %                          d=sqrt((node(CM(i)).location(1)-node(CH(j)).location(1))^2+((node(CM(i)).location(2)-node(CH(j)).location(2)))^2);
            %                          if d<d_min
            %                             d_min=d;
            %                             node(CM(i)).CH=CH(j);
            %                          end
            %                      end
            %                      plot([node(CM(i)).location(1),node(node(CM(i)).CH).location(1)],[node(CM(i)).location(2),node(node(CM(i)).CH).location(2)],'b');hold on
            %                 end
            %                 pause(0.1);
            %                 hold off;
            %             end
            
            %%%%簇头广播信息消耗能量（包含能量消耗）%%%%
            for i=1:count_CH
                d=AD_Range;
                if d>d0
                    node(CH(i)).energy=node(CH(i)).energy-Lc*(Eelec+Emp*d^4);
                    contorl_Overhead=contorl_Overhead+Lc*(Eelec+Emp*d^4);%控制开销
                else
                    node(CH(i)).energy=node(CH(i)).energy-Lc*(Eelec+Efs*d^2);
                    contorl_Overhead=contorl_Overhead+Lc*(Eelec+Efs*d^2);%控制开销
                end
            end
            
            %%%%节点根据距离选择簇头加入（包含能量消耗）%%%%
            for i=1:length(AliveNodes)
                d_min=inf;%判断加入的条件
                F=find(CH==AliveNodes(i));
                if length(F)==0%如果节点不是是簇头
                    for j=1:count_CH
                        d=sqrt((node(AliveNodes(i)).location(1)-node(CH(j)).location(1))^2+(node(AliveNodes(i)).location(2)-node(CH(j)).location(2))^2);
                        if d<=AD_Range %如果节点在簇头广播范围内
                            node(AliveNodes(i)).energy=node(AliveNodes(i)).energy-Lc*Eelec;%节点接收信息的能耗
                            contorl_Overhead=contorl_Overhead+Lc*Eelec;%控制开销
                            if d<d_min %判断簇头和节点距离是否最近
                                d_min=d;
                                node(AliveNodes(i)).CH=CH(j);%为节点选择簇头
                            end
                        end
                    end
                    if  node(AliveNodes(i)).CH~=0 %如果节点拥有簇头
                        node(node(AliveNodes(i)).CH).CM=[node(node(AliveNodes(i)).CH).CM  AliveNodes(i)];% 将节点添加到簇头当中
                        node(node(AliveNodes(i)).CH).CM_count=node(node(AliveNodes(i)).CH).CM_count+1;% 对应簇头的簇成员数量+1
                        % 画图（簇头与簇成员节点连线） %
                        %                     plot([node(i).location(1),node(node(i).CH).location(1)],[node(i).location(2),node(node(i).CH).location(2)],'b');hold on
                    end
                end
            end
            
            %%%%节点返回注册信息（包含能量消耗） %%%%
            for i=1:length(AliveNodes)
                F=find(CH==AliveNodes(i));
                if length(F)==0 && node(AliveNodes(i)).CH~=0 %如果节点不是是簇头、拥有簇头
                    d=sqrt(((node(AliveNodes(i)).location(1)-node(node(AliveNodes(i)).CH).location(1)))^2+((node(AliveNodes(i)).location(2)-node(node(AliveNodes(i)).CH).location(2)))^2);
                    if d>d0
                        node(AliveNodes(i)).energy=node(AliveNodes(i)).energy-Lc*(Eelec+Emp*d^4);
                        contorl_Overhead=contorl_Overhead+Lc*(Eelec+Emp*d^4);%控制开销
                    else
                        node(AliveNodes(i)).energy=node(AliveNodes(i)).energy-Lc*(Eelec+Efs*d^2);
                        contorl_Overhead=contorl_Overhead+Lc*(Eelec+Efs*d^2);%控制开销
                    end
                end
            end
            
            %%%%簇头发送TDMA时间表（包含能量消耗）%%%%
            for i=1:count_CH%找到距离最远的节点，以此距离广播TDMA时间表
                d_Max=-inf;
                for j=1:node(CH(i)).CM_count%找到距离最远的节点
                    id=node(CH(i)).CM(j);
                    d=sqrt(((node(id).location(1)-node(CH(i)).location(1)))^2+((node(id).location(2)-node(CH(i)).location(2)))^2);
                    if d>d_Max
                        d_Max=d;
                    end
                end
                if node(CH(i)).CM_count>0%如果簇内只有簇头，则不广播
                    if d_Max>d0
                        node(CH(i)).energy=node(CH(i)).energy-Lc*(Eelec+Emp*d_Max^4);
                        contorl_Overhead=contorl_Overhead+Lc*(Eelec+Emp*d_Max^4);%控制开销
                    else
                        node(CH(i)).energy=node(CH(i)).energy-Lc*(Eelec+Efs*d_Max^2);
                        contorl_Overhead=contorl_Overhead+Lc*(Eelec+Efs*d_Max^2);%控制开销
                    end
                end
            end
            
            %%%%节点接收TDMA时间表（包含能量消耗）%%%%
            for i=1:length(AliveNodes)
                F=find(CH==AliveNodes(i));
                if length(F)==0  && node(AliveNodes(i)).CH~=0%如果节点不是簇头、能量大于0、拥有簇头
                    node(AliveNodes(i)).energy=node(AliveNodes(i)).energy-Lc*Eelec;
                    contorl_Overhead=contorl_Overhead+Lc*Eelec;%控制开销
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%% ②稳态数据传输阶段 %%%%%%%%%%%%%%%%%%%%%%%%%%%
        for z=1:1
            %%%%%%%%%%%% 数据传送(消耗能量) %%%%%%%%%%%%
            for i=1:count_CH
                for j=1:node(CH(i)).CM_count%簇内节点传送数据
                    %%%%%%%%%%%%对应节点发送数据消耗能量%%%%%%%%%%%%
                    id=node(CH(i)).CM(j);%节点序号
                    d=sqrt(((node(id).location(1)-node(CH(i)).location(1)))^2+((node(id).location(2)-node(CH(i)).location(2)))^2);
                    if d>d0
                        node(id).energy=node(id).energy-Ld*(Eelec+Emp*d^4);
                    else
                        node(id).energy=node(id).energy-Ld*(Eelec+Efs*d^2);
                    end
                    %%%%%%%%%%%%对应簇头接收数据能耗%%%%%%%%%%%%
                    node(CH(i)).energy=node(CH(i)).energy-Ld*Eelec;
                end
            end
            %%%%%%%%%%%%簇头将簇内数据融合并将数据送往基站（能量消耗）%%%%%%%%%%%%
            for i=1:count_CH
                node(CH(i)).energy=node(CH(i)).energy-node(CH(i)).CM_count*Ld*EDA;%数据融合
                d=sqrt((node(CH(i)).location(1)-sinkx)^2+((node(CH(i)).location(2)-sinky)^2));%将数据传送至基站
                if d>d0
                    node(CH(i)).energy=node(CH(i)).energy-Ld*(Eelec+Emp*d^4);
                else
                    node(CH(i)).energy=node(CH(i)).energy-Ld*(Eelec+Efs*d^2);
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%更新参数%%%%%%%%%%%%%%%%%%%%%%%
        for z=1:1
            %%%%%%%%%%%%更新节点参数%%%%%%%%%%%%
            for i=1:N
                node(i).CH=0;% 表明节点的簇头
                node(i).CM=[];% 表明节点的簇成员
                node(i).CM_count=0;% 表明簇内的簇成员节点数
            end
            
            %%%%%%%%%%%%更新性能指标参数%%%%%%%%%%%%
            r=r+1%轮数+1
            
            %%%吞吐量%%%
            Throughput=Throughput+count_CH*Ld;%吞吐量
            trace_Throughput(r)=Throughput;
            
            %%%网络生命周期%%%
            number_AliveNodes=N;
            AliveNodes=1:1:N;
            for i=1:N
                if node(i).energy<=0
                    number_AliveNodes=number_AliveNodes-1;%跟新节点数
                    AliveNodes=AliveNodes(~ismember(AliveNodes,i));%删除死亡节点
                end
            end
            network_Lifetime(r)=number_AliveNodes;%存活节点数量
            
            %%%网络剩余能量%%%
            E=0;
            for i=1:N
                if node(i).energy>0
                    E=E+node(i).energy;
                end
            end
            network_ResidualEnergy(r)=E;%网络剩余能量
            
            %%%网络时延%%%
            %             network_Delay(r)=delay;%网络时延
            
            %%%控制开销%%%
            network_Control_Overhead(r)=contorl_Overhead;%控制开销
            
            %%%簇头个数%%%
            network_CH_Count(r)=count_CH;%簇头个数
            
            %%%参数更新%%%
            %             delay=0;%网络时延初始化
            %             contorl_Overhead=0;%控制开销初始化
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4.实验结果  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for zh=1:1
%     figure(3)
%     plot(network_Lifetime,'-o','MarkerIndices',maker_idx,'linewidth',1.25) ;hold on
%     figure(4)
%     plot(network_ResidualEnergy,'-o','MarkerIndices',maker_idx,'linewidth',1.25) ;hold on
%     figure(5)
%     plot(trace_Throughput,'-o','MarkerIndices',maker_idx,'linewidth',1.25) ;hold on
%     figure(6)
%     plot(network_Control_Overhead,'-o','MarkerIndices',maker_idx,'linewidth',1.25) ;hold on

%     figure(3)
%     plot(network_Lifetime,':k','linewidth',1) ;hold on
%     figure(4)
%     plot(network_ResidualEnergy,':k','linewidth',1) ;hold on
%     figure(5)
%     plot(trace_Throughput,':k','linewidth',1) ;hold on
%     figure(6)
%     plot(network_Control_Overhead,':k','linewidth',1) ;hold on
    
    figure(3)
    plot(network_Lifetime,':k','linewidth',1.25) ;hold on
    figure(4)
    plot(network_ResidualEnergy,':k','linewidth',1.25) ;hold on
%     figure(5)
%     plot(trace_Throughput,':k','linewidth',1.25) ;hold on
%     figure(6)
%     plot(network_Control_Overhead,':k','linewidth',1.25) ;hold on
%     
%      figure(3)
%     plot(network_Lifetime,'.k','linewidth',0.5) ;hold on
%     figure(4)
%     plot(network_ResidualEnergy,'.k','linewidth',0.5) ;hold on
%     figure(5)
%     plot(trace_Throughput,'.k','linewidth',0.5) ;hold on
%     figure(6)
%     plot(network_Control_Overhead,'.k','linewidth',0.5) ;hold on
end



