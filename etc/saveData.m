function saveData(fileName,data,precision,extenIn)
  
  exten = 'bin';
  if(nargin>3), exten = extenIn;end
    
  fileName = [fileName '_' precision '.' exten];
  switch exten
   case 'bin'
    fId = fopen(fileName,'w'); 
    fwrite(fId, data, precision); 
    fclose(fId);
    
   case 'txt'
    if(strcmp(precision,'single'))
      precision = [];
    else
      precision = ['-' precision];
    end
      
    save(fileName,'data',['-ascii' precision]);
    
   otherwise
    disp('The format can not be deduced for the file extension');
  end
