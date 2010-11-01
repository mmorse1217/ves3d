function out = ReadBinary(name, precision)

if ( nargin < 2 )
  precision = 'double';
end

fid = fopen(name,'r');
out = fread(fid,precision);
fclose(fid);
