mkdir('../pictures')
d = dir('../orbit_results');
isub = [d(:).isdir]; % returns logical vector
name_folds = {d(isub).name}';
name_folds(ismember(name_folds,{'.','..'})) = [];
row = 4;
col = 5;
figure
set(gcf,'Position',[10 10 2000 2000]);
for i=1:20
    % 每一个起始磁面循环
    lost_list = zeros(1,10);
    diff_list = zeros(1,10);
    for j = 1:10
        % 统计10个俯仰角对应的信息
        diff = ['../orbit_results/',name_folds{20*j-20+i},'/orbit_results/diffusion.plt'];
        diff_mean = diff_read(diff);
        diff_list(j) = diff_mean * 1e6;
    end
    subplot(row,col,i)
    plot(linspace(0,1,10),diff_list,'color',[255,140,105]/255);
    xlabel('lambda');
    ylabel('$< ( \delta \psi _p ) ^2 > /10^{-6}$','interpreter','latex')
    polo = linspace(0,1,20);
    polo = polo(i);
    polo = sprintf('%4.3f',polo);
    title(['polo=',polo]);
end 
saveas(gcf,'../pictures/diff_rip_NTM.fig')
saveas(gcf,'../pictures/diff_rip_NTM.png')