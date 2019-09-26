mkdir('../pictures')
d = dir('../orbit_results');
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
row = 5;
col = 6;
for i=1:length(nameFolds)
    poin_dir = ['../orbit_results/',nameFolds{i},'/orbit_results/poincare.plt'];
    wall_dir = ['../orbit_results/',nameFolds{i},'/orbit_results/wall.plt'];
    name = nameFolds{i};
    name = name(10:end);
    poinxz(row,col,i,wall_dir,poin_dir,name);
end 
saveas(gcf,'../pictures/poinxz.fig')
saveas(gcf,'../pictures/poinxz.png')