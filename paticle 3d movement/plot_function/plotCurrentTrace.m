function plotCurrentTrace(pos, t)
    global trace x_min x_max y_min y_max z_min z_max;
    if t ~= 0
%         D: ά����
        [~,D,~] = size(pos);
        hold on
        axis equal
        box on
        grid on
        set(gca, 'linewidth', 1.5)
%         ����trace������ʼ����,���trace(֮ǰ�Ĳ���)����t(��ǰ����),����Ϊ0
        start = max(t-trace,1);
        for j = 1:D
            eval(['d',num2str(j),'(:)=pos(1,j,start:t);'])
        end
        if D == 3
            axis([x_min x_max y_min y_max z_min z_max])
            plot3(d1, d2, d3, 'linewidth', 1.5)
        elseif D == 2
            axis([x_min x_max y_min y_max])
            plot(d1, d2, 'linewidth', 1.5)
        end
    end
end