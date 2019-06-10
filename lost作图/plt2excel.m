function []=plt2excel(filename)
data = importdata([filename,'.plt']);
xlswrite(['./',filename,'.xlsx'],data.data)
end