template<typename T>
T ParserLite(const char *filename)
{
    ifstream in_file(filename);
    T class_instance;
    
    if ( in_file.is_open() ) {
        
        string line, var_name;
        
        while ( getline ( in_file, line ) ) {
            string::size_type i = line.find_first_not_of ( " \t\n\v" );
            
            if ( i == string::npos || (i !=string::npos && line[i] == '#') )
                continue;
            
            char *cstr = new char [line.size()+1];
                        
            i = line.find_first_of ( "=" );
            var_name = line.substr(0,i);
            line = line.substr(++i);
            
            sscanf(var_name.c_str(),"%s ",cstr);
            class_instance.SetMember(cstr,line);

            delete cstr;
        }   
    }
    return class_instance;
}

template <typename T>
T String2Num(T &num, const string &s)
{
    std::istringstream iss(s);
    isss >> num;
    return (num);
}
