function M=mass_matrixQ1_2D(elements,areas,coeffs)
%coeffs can be either P0 (elementwise constant) or P1 (elementwise nodal) function 
%represented by a collumn vector with size(elements,1) or size(coordinates,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally

NE=size(elements,1); %number of elements
%DIM=size(coordinates,2); %problem dimension

%particular part for a given element in a given dimension
NLB=4; %number of local basic functions, it must be known!
%coord=zeros(DIM,NLB,NE);
%for d=1:DIM
%    for i=1:NLB
%        coord(d,i,:)=coordinates(elements(:,i),d);
%    end
%end   

ips=intrec(3);
phi=shapefun (ips.coord','Q1');
Y=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);

%copy this part for a creation of a new element
M_local=zeros(NLB);
for i=1:ips.nip
    M_local=M_local+ips.w(i)*phi(:,i)*phi(:,i)';
end

if (nargin<3)
    Z=astam(areas,reshape(repmat(M_local,1,NE),NLB,NLB,NE));
else
    if numel(coeffs)~=size(elements,1)  %P1->P0 averaging
        coeffs=evaluate_average_point(elements,coeffs);
    end      
    Z=astam(areas.*coeffs,reshape(repmat(M_local,1,NE),NLB,NLB,NE));
end
X=permute(Y,[2 1 3]);  
M=sparse(X(:),Y(:),Z(:));






    end
