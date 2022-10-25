function h=show_nodal_scalar_frame(nodalValue,coordinates,elements,m)

nodalDisplacement=0*coordinates;

if size(coordinates,2)==2
    X=coordinates(:,1)+nodalDisplacement(:,1);
    Y=coordinates(:,2)+nodalDisplacement(:,2);
    
    h=trisurf(elements,X,Y,nodalValue,'FaceColor','flat','LineWidth',0.1,'EdgeColor','w');
    switch m
        case 1
            title('Solution')
        case 2
            title('Solution')
        case 3
            title('Solution')
        case 4
            title('Solution')
    end
    
elseif size(coordinates,2)==3  %3D only testing
    
    X=coordinates(:,1)+nodalDisplacement(:,1);
    Y=coordinates(:,2)+nodalDisplacement(:,2);
    Z=coordinates(:,3)+nodalDisplacement(:,3);
    
    elements=[elements(:,[1 2 3]); elements(:,[1 2 4]); elements(:,[2 3 4]); elements(:,[3 1 4])];
    
    h=trisurf(elements,X,Y,Z,nodalValue,'FaceColor','interp');
end

%set(gcf,'renderer','zbuffer');