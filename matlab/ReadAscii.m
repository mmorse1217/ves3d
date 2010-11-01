function out = ReadAscii(name,count)

if ( nargin < 2 )
  count = {};
else
  count = {count};
end

fid = fopen(name,'r');
out = fscanf(fid,'%g',count{:});
fclose(fid);
