time = linspace(0,100,1);
y = sin(time);
writerObj = VideoWriter('test_plt.avi');
open(writerObj);
myMovie(1:length(time)) = struct('cdata', [], 'colormap', []);
for t = 1:length(time)
    plot(time(1:t),y(1:t));
    frame = getframe;
    frame.cdata = imresize(frame.cdata, [685, 685]);
    writeVideo(writerObj, frame);
end
close(writerObj)