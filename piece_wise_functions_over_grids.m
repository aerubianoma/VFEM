x = -2:0.8:2;
y = -2:0.8:2;
[X,Y] = meshgrid(x,y);
Z = zeros(6);
Z(:,1:3) = exponential(X(:,1:3),Y(:,1:3));
Z(:,4:6) = combination(X(:,4:6),Y(:,4:6));
T = delaunay(X,Y);
trimesh(T,X,Y,Z)

function piece_wise_one = exponential(X,Y)
    
    piece_wise_one = zeros(6,3);
    piece_wise_one = (exp(X)+exp(Y))*0+1;

end

function piece_wise_two = combination(X,Y)
    
    piece_wise_two = zeros(6,3);
    piece_wise_two = (X + 10*Y)*0 + 1;

end