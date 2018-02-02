function mkdirIfNotExist( dirpath )
    
    if ( dirpath(end) ~= '/' )
        dirpath = [dirpath '/'];
    end
    
    if ( exist(dirpath, 'dir') == 0 )
        mkdir(dirpath);
    else
        rmdir(dirpath, 's');
        mkdir(dirpath);
    end
    
end