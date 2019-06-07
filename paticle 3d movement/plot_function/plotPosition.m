function plotPosition(pos, time)
    global isOutVideo;
    figure(1)
    if isOutVideo == true
        writerObj = VideoWriter('test_plt.avi');
        open(writerObj);
        myMovie(1:length(time)) = struct('cdata', [], 'colormap', []);
    end
    for t = 1:length(time)
        plotPosVec(pos(:,:,t), t, pos)
        if isOutVideo == true
            frame = getframe;
            frame.cdata = imresize(frame.cdata, [685, 685]);
            writeVideo(writerObj, frame);
        end
    end
    if isOutVideo == true
        close(writerObj);
    end
end