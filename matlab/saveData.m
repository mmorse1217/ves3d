function count = writeData(fileName, data, precision)
  
  fid = fopen(fileName, 'w');
  if(fid>0)
    count = fwrite(fid, data, precision);
    fclose(fid);
  else
    count = fid;
  end