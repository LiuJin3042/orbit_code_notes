function plotPosition(pos, time)
%���������˶��ĺ���
    global isOutVideo;
    figure(1)
    if isOutVideo == true
        %��һ����Ƶ����
        writerObj = VideoWriter('test_plt.avi');
        open(writerObj);
    end
    for t = 1:length(time)
        %�������ӵ�ǰλ��,��������λ��,֮ǰ�Ĺ켣
        plotPosVec(pos(:,:,t), t, pos)
        if isOutVideo == true
            %����Ϊ��Ƶ
            frame = getframe;
            frame.cdata = imresize(frame.cdata, [685, 685]);
            writeVideo(writerObj, frame);
        end
    end
    if isOutVideo == true
        close(writerObj);
    end
end