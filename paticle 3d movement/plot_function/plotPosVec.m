function plotPosVec(current_pos, t, pos_all)
%current_pos, 是当前粒子位置,当前时间
%pos_all, 是粒子所有位置,也就是主程序的pos
    global dt time_start duration x_min x_max y_min y_max z_min z_max x_label y_label z_label;
    fprintf('present time is %fs/%fs)\n', time_start + t*dt, time_start + duration)
    [~,D] = size(current_pos);
    figure(1)
    if D == 2
        scatter(current_pos(:,1), current_pos(:,2), 'or', 'filled')
    elseif D == 3
        scatter3(current_pos(:,1), current_pos(:,2), current_pos(:,3), 'or', 'filled')
    end
    axis equal
    box on
    grid on
    set(gca, 'linewidth', 1.5, 'xtick', floor(linspace(x_min, x_max, 11)), 'ytick', floor(linspace(y_min, y_max, 11)), 'ztick', floor(linspace(z_min, z_max, 11)))
    hold on
%     画历史轨迹
    plotCurrentTrace(pos_all, t)
    text(x_max*13/25, y_min*20/25, z_min*20/25, 'Moving Trace of Single Particle', 'horiz', 'center', 'color', 'r')
    axis([x_min x_max y_min y_max z_min z_max])
    xlabel(x_label)
    ylabel(y_label)
    zlabel(z_label)
    hold off
end