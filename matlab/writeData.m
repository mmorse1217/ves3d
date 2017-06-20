function count = writeData(fileName, suffix, data)
% WRITEDATA(fileName, data) - Project specific interface to write binary files

if(isstr(suffix))
  fileName = ['../data' filesep fileName suffix '.bin'];
else
  fileName = ['../data' filesep fileName num2str(suffix) '.bin'];
end

if(~exist('../data', 'dir'))
  mkdir('../data')
end

fid = fopen(fileName, 'w');
if(fid>0)
	count = fwrite(fid, data, 'double');
	fclose(fid);
else
	count = fid;
end
