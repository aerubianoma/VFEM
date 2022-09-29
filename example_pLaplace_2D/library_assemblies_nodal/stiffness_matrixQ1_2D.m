function [K,areas]=stiffness_matrixQ1_2D(elements,coordinates,coeffs)
%coeffs can be either P0 (elementwise constant) or P1 (elementwise nodal) function 
%represented by a collumn vector with size(elements,1) or size(coordinates,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally

NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension

%particular part for a given element in a given dimension
NLB=4; %number of local basic functions, it must be known!
coord=zeros(DIM,NLB,NE);
for d=1:DIM
    for i=1:NLB
        coord(d,i,:)=coordinates(elements(:,i),d);
    end
end   

%area
areas1=areas_triangle(elements,coordinates);
areas2=areas_triangle(elements(:,[1 3 4]),coordinates);
areas=areas1+areas2;

ips=intrec(2);
dphi = phider(coord,ips.coord','Q1'); 
dphi = squeeze(dphi); 

product=zeros(NLB,NLB,NE);
for i=1:ips.nip
    product=product+ips.w(i)*amtam(squeeze(dphi(:,:,i,:)),squeeze(dphi(:,:,i,:)));
end 

if (nargin<3)
    Z=astam(areas,product); 
else
    if numel(coeffs)==size(coordinates,1)  %P1->P0 averaging
        coeffs=evaluate_average_point(elements,coeffs);
    end  
    Z=astam((areas.*coeffs),product);
end
Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);

%copy this part for a creation of a new element
X=permute(Y,[2 1 3]);
K=sparse(X(:),Y(:),Z(:));  
      
    function areas=areas_triangle(elements,coordinates)
        a=coordinates(elements(:,1),1);
        b=coordinates(elements(:,2),1);
        c=coordinates(elements(:,3),1);
        d=coordinates(elements(:,1),2);
        e=coordinates(elements(:,2),2);
        f=coordinates(elements(:,3),2);
        areas=abs(b.*f - c.*e - a.*f + c.*d + a.*e - b.*d)/2; 
    end

    end
