clear;clc;close all
Function_name='F5'; % 使用方程的名字，对应 Functions_details 文件
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);  %得到具体的方程即目标函数，变量的维度，变量的上下限
pop_num=40;  % Number of search agents 种群数量
Max_iter=501;    % Maximum numbef of iterations 最大迭代次数
%以下是各种优化算法的比较
Time_compare=[];      %算法的运行时间比较
Fival_compare=[];       %算法的最终目标比较
curve_compare=[];     %算法的过程函数比较
%麻雀搜索算法
name_all=[];     %算法的名称记录
%%
iter=1;
 % 灰狼优化算法 
t1=clock;
[fMin_GWO,bestX_GWO,GWO_curve]=GWO(pop_num,Max_iter,lb,ub,dim,fobj);      % 灰狼优化算法 
t2=clock;
time_GWO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_GWO];
Time_compare=[Time_compare,time_GWO(end)];
curve_compare=[curve_compare;GWO_curve];
name_all{1,iter}='GWO';
iter=iter+1;


% 改进灰狼优化算法 
t1=clock;
[fMin_IGWO,bestX_IGWO,IGWO_curve]=IGWO(pop_num,Max_iter,lb,ub,dim,fobj);  
t2=clock;
time_IGWO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_IGWO];
Time_compare=[Time_compare,time_IGWO(end)];
curve_compare=[curve_compare;IGWO_curve];
name_all{1,iter}='IGWO';
iter=iter+1;
%鲸鱼优化算法
t1=clock;
[fMin_WOA,bestX_WOA,WOA_curve]=WOA(pop_num,Max_iter,lb,ub,dim,fobj); % 鲸鱼优化算法    
t2=clock;
time_WOA=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_WOA];
Time_compare=[Time_compare,time_WOA(end)];
curve_compare=[curve_compare;WOA_curve];
name_all{1,iter}='WOA';
iter=iter+1;
%%%%%多版本优化器  Multi-Verse Optimizer
t1=clock;
[fMin_MVO,bestX_MVO,MVO_curve]=MVO(pop_num,Max_iter,lb,ub,dim,fobj);     
t2=clock;
time_MVO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_MVO];
Time_compare=[Time_compare,time_MVO(end)];
curve_compare=[curve_compare;MVO_curve];
name_all{1,iter}='MVO';
iter=iter+1;
%北方苍鹰优化算法Northern_Goshawk_Optimization_A_New_Swarm-Based_Algorithm
t1=clock;
[fMin_NGO,bestX_NGO,NGO_curve]=NGO(pop_num,Max_iter,lb,ub,dim,fobj);     
t2=clock;
time_NGO=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_NGO];
Time_compare=[Time_compare,time_NGO(end)];
curve_compare=[curve_compare;NGO_curve];
name_all{1,iter}='NGO';
iter=iter+1;
% 阿狸八八和四十大盗优化算法 Ali baba and the Forty Thieves (AFT) algorithm 
t1=clock;
[fMin_AFT,bestX_AFT,AFT_curve]=AFT(pop_num,Max_iter,lb,ub,dim,fobj);     
t2=clock;
time_AFT=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
Fival_compare=[Fival_compare,fMin_AFT];
Time_compare=[Time_compare,time_AFT(end)];
curve_compare=[curve_compare;AFT_curve];
name_all{1,iter}='AFT';
iter=iter+1;


























%%
% load('color_list')
% figure(2)
% color_all=color_list(randperm(length(color_list)),:);
% %画迭代过程曲线
% for N=1:length(Fival_compare)
%      semilogy(curve_compare(N,:),'Color',color_all(N,:),'LineWidth',2)
%      hold on
% end
% xlabel('迭代次数');
% ylabel('目标函数值');
% grid on
% box on
% legend(name_all)
% % %%  运行值和最终目标函数比较
% figure(3)
% color=color_list(randperm(length(color_list)),:);
% width=0.7; %柱状图宽度
% for  i=1:length(Fival_compare) 
%    set(bar(i,Fival_compare(i),width),'FaceColor',color(i,:),'EdgeColor',[0,0,0],'LineWidth',2)
%    hold on
%     
%    %在柱状图 x,y 基础上 绘制误差 ,low为下误差，high为上误差，LineStyle 误差图样式，'Color' 误差图颜色  
%    % 'LineWidth', 线宽,'CapSize',误差标注大小
% %    errorbar(i, y(i), low(i), high(i), 'LineStyle', 'none', 'Color', color(i+3,:), 'LineWidth', 1.5,'CapSize',18);
% end
% ylabel('obj-value')
% ax=gca;
% ax.XTick = 1:1:length(Fival_compare);
% set(gca,'XTickLabel',name_all,"LineWidth",2);
% set(gca,"FontName","Times New Roman","FontSize",12,"LineWidth",2)
% %% 运行时间比较
% figure(4)
% color=color_list(randperm(length(color_list)),:);
% width=0.7; %柱状图宽度
% for  i=1:length(Fival_compare) 
%    set(bar(i,Time_compare(i),width),'FaceColor',color(i,:),'EdgeColor',[0,0,0],'LineWidth',2)
%    hold on
%    %在柱状图 x,y 基础上 绘制误差 ,low为下误差，high为上误差，LineStyle 误差图样式，'Color' 误差图颜色  
%    % 'LineWidth', 线宽,'CapSize',误差标注大小
% %    errorbar(i, y(i), low(i), high(i), 'LineStyle', 'none', 'Color', color(i+3,:), 'LineWidth', 1.5,'CapSize',18);
% end
% ylabel('Time')
% ax=gca;
% ax.XTick = 1:1:length(Fival_compare);
% set(gca,'XTickLabel',name_all,"LineWidth",2);
% set(gca,"FontName","Times New Roman","FontSize",12,"LineWidth",2)
% %%
% % 画出目标函数的图示,只能画到三维，目标函数dim的设置可能很多维
% figure(1)
% x=lb/2:0.1:ub/2;  %x轴,y轴范围
% y=x;
% L_num=length(x);
% f_value=[];      %对应x,y的函数值
% for i=1:L_num
%     for j=1:L_num
%         f_value(i,j)=fobj([x(i),y(j)]);
%     end
% end
% surfc(x,y,f_value,'LineStyle','none');