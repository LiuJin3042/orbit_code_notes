function plotPosition(pos, time)
%画出粒子运动的函数
    global dt isOutVideo;
    figure(1)
    if isOutVideo == true
        %打开一个视频对象
        writerObj = VideoWriter('test_plt.avi');
        writerObj.FrameRate = floor(1/dt);
        open(writerObj);
    end
    for t = 1:length(time)
        %画出粒子当前位置,包括粒子位置,之前的轨迹
        plotPosVec(pos(:,:,t), t, pos)
        if isOutVideo == true
            %保存为视频
            frame = getframe;
            frame.cdata = imresize(frame.cdata, [1000,1000]);
            writeVideo(writerObj, frame);
        end
    end
    if isOutVideo == true
        close(writerObj);
    end
end