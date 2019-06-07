function plotPosVec(current_pos, t, pos_all)
%current_pos, 是当前粒子位置,当前时间
%pos_all, 是粒子所有位置,也就是主程序的pos
    global x_min x_max y_min y_max z_min z_max;
    figure(1)
    scatter3(current_pos(:,1), current_pos(:,2), current_pos(:,3), 'ok', 'filled')
    axis equal
    box on
    grid on
    set(gca, 'linewidth', 1.5, 'xtick', floor(linspace(x_min, x_max, 11)), 'ytick', floor(linspace(y_min, y_max, 11)), 'ztick', floor(linspace(z_min, z_max, 11)))
    hold on
%     画速度矢量的程序,省略了
%     for i = 1:length(vx)
%         line([pos(i,1) pos(i,1)+vx(i)/2], [pos(i,2), pos(i,2)+vy(i)/2], 'linewidth', 1.2)
%     end
%     画历史轨迹
    plotCurrentTrace(pos_all, t)
    text(x_max*13/25, y_min*20/25, z_min*20/25, 'Moving Trace of Single Particle', 'horiz', 'center', 'color', 'r')
    axis([x_min x_max y_min y_max z_min z_max])
    hold off
end