% 给每一个运行结果都画一张xz图片
clear
figure
mkdir('../pictures')
d = dir('../orbit_results');
isub = [d(:).isdir]; % returns logical vector
name_folds = {d(isub).name}';
name_folds(ismember(name_folds,{'.','..'})) = [];
var1_num = 2; % 第一个自变量的数目, 这里是磁面
var2_num = 6; % 第二个自变量的数目, 这里是磁面

% poin_dir = './poincare.plt';
% poincare = load(poin_dir);
% poin_x=poincare(:,3);
% poin_z=poincare(:,4);

for i=1:var1_num
    % 每一个起始磁面循环
    lost_list = zeros(1,var2_num);
    diff_list = zeros(1,var2_num);
    for j = 1:var2_num
        hold on
        
        %% 画器壁
        path = ['../orbit_results/',name_folds{(i-1)*var2_num+j},'/orbit_results/wall.plt']; 
        m=importdata(path,' ',2);
        wall=m.data;
        xw=wall(:,1);
        zw=wall(:,2);
        plot(xw,zw)
        
%         %% 画轴线
%         path = ['../orbit_results/',name_folds{(i-1)*var2_num+j},'/orbit_results/axis.plt']; 
%         m=importdata(path,' ',2);
%         axis=m.data;
%         x1=axis(:,1);
%         y1=axis(:,2);
%         x2=axis(:,3);
%         y2=axis(:,4);
%         plot(x1,y1)
%         plot(x2,y2)

        %% 画轨迹
        path = ['../orbit_results/',name_folds{(i-1)*var2_num+j},'/orbit_results/traj1.plt']; 
        m=importdata(path,' ',2);
        traj1=m.data;
        t=traj1(:,1);
        x=traj1(:,2);
        z=traj1(:,3);
        pol=traj1(:,4);
        thet=traj1(:,5);
        zet=traj1(:,6);
        plot(x,z,'r-')
        
%        %% 画磁面
%         path = ['../orbit_results/',name_folds{(i-1)*var2_num+j},'/orbit_results/surf.plt']; 
%         m=importdata(path,' ',1);
%         surf=m.data;
%         x01=surf(:,1);
%         z01=surf(:,2);
%         x02=surf(:,3);
%         z02=surf(:,4);
%         x03=surf(:,5);
%         z03=surf(:,6);
%         x04=surf(:,7);
%         z04=surf(:,8);
%         x05=surf(:,9);
%         z05=surf(:,10);
%         x06=surf(:,11);
%         z06=surf(:,12);
%         x07=surf(:,13);
%         z07=surf(:,14);
%         x08=surf(:,15);
%         z08=surf(:,16);
%         plot(x01,z01)
%         plot(x02,z02)
%         plot(x03,z03)
%         plot(x04,z04)
%         plot(x05,z05)
%         plot(x06,z06)
%         plot(x07,z07)
%         plot(x08,z08)
%         
%         %% 画spokes
%         path = ['../orbit_results/',name_folds{(i-1)*var2_num+j},'/orbit_results/spokes.plt'];         
%         m=importdata(path,' ');
%         spokes=m;
%         xx1=spokes(:,1);
%         zz1=spokes(:,2);
%         xx2=spokes(:,3);
%         zz2=spokes(:,4);
%         xx3=spokes(:,5);
%         zz3=spokes(:,6);
%         xx4=spokes(:,7);
%         zz4=spokes(:,8);
%         xx5=spokes(:,9);
%         zz5=spokes(:,10);
%         xx6=spokes(:,11);
%         zz6=spokes(:,12);
%         xx7=spokes(:,13);
%         zz7=spokes(:,14);
%         xx8=spokes(:,15);
%         zz8=spokes(:,16);
%         xx9=spokes(:,17);
%         zz9=spokes(:,18);
%         xx10=spokes(:,19);
%         zz10=spokes(:,20);
%         xx11=spokes(:,21);
%         zz11=spokes(:,22);
%         xx12=spokes(:,23);
%         zz12=spokes(:,24);
%         xx13=spokes(:,25);
%         zz13=spokes(:,26);
%         xx14=spokes(:,27);
%         zz14=spokes(:,28);
%         xx15=spokes(:,29);
%         zz15=spokes(:,30);
%         xx16=spokes(:,31);
%         zz16=spokes(:,32);
%         plot(xx1,zz1)
%         plot(xx2,zz2)
%         plot(xx3,zz3)
%         plot(xx4,zz4)
%         plot(xx5,zz5)
%         plot(xx6,zz6)
%         plot(xx7,zz7)
%         plot(xx8,zz8)
%         plot(xx9,zz9)
%         plot(xx10,zz10)
%         plot(xx11,zz11)
%         plot(xx12,zz12)
%         plot(xx13,zz13)
%         plot(xx14,zz14)
%         plot(xx15,zz15)
%         plot(xx16,zz16)
%         
       %% 画庞加莱
%         plot(poin_x,poin_z,'.','MarkerSize',1.5)
        hold off
        xlabel('$X/cm$','interpreter','latex')
        ylabel('$Z/cm$','interpreter','latex')
        path = ['../pictures/',name_folds{(i-1)*var2_num+j}]; 
        saveas(gcf,[path,'.fig'])
        saveas(gcf,[path,'.png'])
        close gcf
    end
end 





