function diff_mean = diff_read(diff)
diff_file = importdata(diff,' ',2);
if iscell(diff_file)
    diff_mean = -1;
elseif isstruct(diff_file)
    diff_mean = mean(diff_file.data(:,2));
end