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
    % ÿһ����ʼ����ѭ��
    lost_list = zeros(1,10);
    for j = 1:10
        % ͳ��10�������Ƕ�Ӧ����Ϣ
        % read orbit.out and substract lost particles
        orbout = ['../orbit_results/',name_folds{20*j-20+i},'/orbit_results/orbit.out']; 
        orbit_out_file = fileread(orbout);
        expr = 'total lost\s*(\d*)';
        lost_number = regexp(orbit_out_file,expr,'tokens');
        % ��ƥ�䵽����ʧ�������ŵ�lost_list
        lost_list(j) = str2num(lost_number{1}{1});
    end
    subplot(row,col,i)
    plot(linspace(0,1,10),lost_list)
    ylabel('lost particles')
    xlabel('lambda');
    polo = linspace(0,1,20);
    polo = polo(i);
    polo = sprintf('%4.3f',polo);
    title(['polo=',polo]);
end 
saveas(gcf,'../pictures/lost_rip_NTM.fig')
saveas(gcf,'../pictures/lost_rip_NTM.png')