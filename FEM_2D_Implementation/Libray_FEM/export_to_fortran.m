% https://stackoverflow.com/questions/48396639/read-matlab-2d-array-to-intel-fortran-and-write-from-intel-fortran-to-matlab-fil

cd matrices_fortran/

fileID = fopen('elems2nodes.bin','w');
fwrite(fileID,mesh.elems2nodes,'double');
fclose(fileID);

fileID = fopen('nodes2coord.bin','w');
fwrite(fileID,mesh.nodes2coord,'double');
fclose(fileID);

cd ..


