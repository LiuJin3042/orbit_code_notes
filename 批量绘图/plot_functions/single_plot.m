mkdir('../pictures')
d = dir('../orbit_results');
isub = [d(:).isdir]; % returns logical vector
name_folds = {d(isub).name}';
name_folds(ismember(name_folds,{'.','..'})) = [];
hold off
for i=1:11
    % 每一个起始磁面循环
    lost_list = zeros(1,10);
    diff_list = zeros(1,10);
    for j = 1:21
        % 统计10个俯仰角对应的信息
        % read orbit.out and substract lost particles
        orbout = ['../orbit_results/',name_folds{20*j-20+i},'/orbit_results/orbit.out']; 
        orbit_out_file = fileread(orbout);
        expr = 'total lost\s*(\d*)';
        lost_number = regexp(orbit_out_file,expr,'tokens');
        % 把匹配到的损失粒子数放到lost_list
        lost_list(j) = str2num(lost_number{1}{1});
        diff = ['../orbit_results/',name_folds{20*j-20+i},'/orbit_results/diffusion.plt'];
        diff_mean = diff_read(diff);
        diff_list(j) = diff_mean * 1e6;
    end
    yyaxis left
    plot(linspace(0,1,10),lost_list)
    ylabel('lost particles')
    xlabel('lambda');
    yyaxis right
    plot(linspace(0,1,10),diff_list);
    ylabel('$< ( \delta \psi _p ) ^2 > /10^{-6}$','interpreter','latex')
    polo = linspace(0,1,20);
    polo = polo(i);
    polo = sprintf('%4.3f',polo);
    title(['polo=',polo]);
    saveas(gcf,['../pictures/',name_folds{20*j-20+i}(1:36),'.png'])
    saveas(gcf,['../pictures/',name_folds{20*j-20+i}(1:36),'.png'])
    close gcf
end 
