function plotCurrentTrace(pos, t)
    global x_min x_max y_min y_max z_min z_max;
    if t ~= 0
%         D: 维度数
        [~,D,~] = size(pos);
%         a: 粒子数
        [a, ~, ~] = size(pos);
        hold on
        axis equal
        box on
        grid on
        set(gca, 'linewidth', 1.5)
        if D == 3
            axis([x_min x_max y_min y_max z_min z_max])
            for i = 1:a
                x = zeros(1, t);
                y = zeros(1, t);
                z = zeros(1, t);
                for j = 1:t
                    x(j) = pos(i, 1, j);
                    y(j) = pos(i, 2, j);
                    z(j) = pos(i, 3, j);
                end
                plot3(x, y, z, 'linewidth', 1.5)
            end
        elseif D == 2
            axis([x_min x_max y_min y_max])
            for i = 1:a
                x = zeros(1, t);
                y = zeros(1, t);
                for j = 1:t
                    x(j) = pos(i, 1, j);
                    y(j) = pos(i, 2, j);
                end
                plot(x, y, 'linewidth', 1.5)
            end
        end
    end
end