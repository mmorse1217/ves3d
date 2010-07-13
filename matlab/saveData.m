function count = writeData(fileName, data, precision)
  
  if(~strcmp(precision,'double')), data = single(data);end
    
  fid = fopen(fileName, 'w');
  if(fid>0)
    count = fwrite(fid, data, precision);
    fclose(fid);
  else
    count = fid;
  end