% This script adds the absolute location of
% the shared functions to the path



%warning('off','MATLAB:rmpath:DirNotFound');

path1 = pwd;
rmpath(genpath(path1));

%cd ..

path2 = pwd;
%if ( isunix )
%    addpath(genpath(path2));
%else
addpath(genpath(path2));
%end

[ver, date] = version;     
%year=str2double(date(end-1:end));

if 0
dir1='library_vectorization';             %original library
dir2='library_vectorization_faster';      %improved library

if (str2double(ver(1:3))>=9.1)       %since version 9.1 there is an aritmetic expansion
    cd(dir1); rmpath(genpath(pwd));
else
    cd(dir2); rmpath(genpath(pwd));
end

end

cd(path1)



