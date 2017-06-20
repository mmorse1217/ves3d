function [val rFlag] = readData(fileInit, suffix, dim)
% READDATA(fileInit, suffix, dim) - Project specific interface to read binary files
  
  nd = prod(dim); % number of entries to be read
  if(isstr(suffix))
    fileName = ['../data' filesep fileInit suffix '.bin'];
  else
    fileName = ['../data' filesep fileInit num2str(suffix) '.bin'];
  end
  rFlag = exist(fileName,'file'); 
  if(rFlag)
	  fid = fopen(fileName,'r');
	  val = fread(fid, nd,'double');
	  fclose(fid);
  
	  if(size(dim,2)>1)
		  val = reshape(val, dim);
	  end
  else
	  val = [];
  end
