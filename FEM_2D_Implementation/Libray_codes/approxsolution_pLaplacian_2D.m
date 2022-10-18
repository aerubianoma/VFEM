function [u_approx, Ju, time, iters] = approxsolution_pLaplacian_2D(exanum,mesh,u_init,all_methods_on,params)

nodesMinim=mesh.nodesFree;  %%%%%

Hstr=mesh.Hstr(nodesMinim,nodesMinim);     %Hessian pattern in free nodes

nnMinim = numel(nodesMinim);

f_constant = params.f_constant;
difference_type = params.difference_type;
eps = params.eps;
p = params.p;

elems2nodes = mesh.elems2nodes;
nodes2elems = mesh.nodes2elems;
dphi_x = mesh.Grad1_elems;
dphi_y = mesh.Grad2_elems;

areas = mesh.areas;

M=mass_matrixP1_2D(mesh.elems2nodes,areas);

f=f_constant*ones(size(u_init));
b = M*f;

if exanum == 1
    fun_g = @energy_exact;
else
    fun_g = @energy_approx;
end

%additional postprocessing
elems_j = cell(nnMinim,1);
elems2nodes_j = cell(nnMinim,1);
Grad_elems_j = cell(nnMinim,2);
for ii=1:nnMinim %components of gradient
    j = nodesMinim(ii);
    elems=nonzeros(nodes2elems(j,:))';
    elems_j{ii} = elems;
    elems2nodes_j{ii}=elems2nodes(elems,:);
    Grad_elems_j{ii,1}=dphi_x(elems,:);
    Grad_elems_j{ii,2}=dphi_y(elems,:);   
end

Dintgrds2 = zeros(nnMinim,1);
for ii=1:nnMinim
    elems = elems_j{ii};
    e2n = elems2nodes_j{ii};
    intgrds2 = f_constant*ones(size(e2n,1),1)/size(e2n,2);
    Dintgrds2(ii) = sum(areas(elems).*intgrds2);
end

[u_nodesMinim, Ju, time, iters] = start_minimizers(fun_g,Hstr,u_init(nodesMinim),all_methods_on);
u_approx = u_init;
u_approx(nodesMinim) = u_nodesMinim;

function [e, g] = energy_approx(u_nodesMinim)
    u = u_init;
    u(nodesMinim) = u_nodesMinim;
    e = energy(u);

    if nargout == 2  % numeric computation of gradient   
        epsilon = u_nodesMinim*eps; 
        epsilon(epsilon==0) = eps;   
        
        switch difference_type
            case 1
                es_plus = energiesLocal(u,+epsilon);
                es_minus = energiesLocal(u,-epsilon);
                d = 2;
            case 2
                
            case 3
                
            otherwise
                error('Bad type of numerical difference')
        end   
        g = (es_plus - es_minus)./epsilon/d;  % g IS COMPUTED DIRECTLY WITHOUT FOR-LOOP
        g=g-Dintgrds2;
    end
end

function e=energy(v)
    v_elems=v(elems2nodes);                           %nodal values
    v_x_elems=sum(dphi_x.*v_elems,2);                 %local x-gradient
    v_y_elems=sum(dphi_y.*v_elems,2);                 %local y-gradient
    intgrds1=(1/p)*sum(abs([v_x_elems v_y_elems]).^p,2);
    e=sum(areas.*intgrds1) - b'*v;
end

function es = energiesLocal(v,epsilon)   
    es = zeros(nnMinim,1);    
    for i=1:nnMinim                                       %loop over internal points
        j=nodesMinim(i);
        v(j)=v(j)+epsilon(i);                             %local elements
        v_elems=v(elems2nodes_j{i});                      %local nodal values
        
        v_x_elems=sum(Grad_elems_j{i,1}.*v_elems,2);      %local x-gradient
        v_y_elems=sum(Grad_elems_j{i,2}.*v_elems,2);      %local y-gradient
           
        intgrds1=(1/p)*(abs(v_x_elems).^p+abs(v_y_elems).^p);
        es(i) = sum(areas(elems_j{i}).*intgrds1);
        v(j)=v(j)-epsilon(i);
    end 
end

function [e, g] = energy_exact(u)
    
    v = u_init;
    v(nodesMinim) = u;
    
    v_elems = v(elems2nodes);
    v_x_elems = sum(dphi_x.*v_elems,2);
    v_y_elems = sum(dphi_y.*v_elems,2);
    intgrds1=(1/p)*sum(abs([v_x_elems v_y_elems]).^p,2);
    intgrds2=f_constant*sum(v_elems,2)/size(v_elems,2);
    e=sum(areas.*(intgrds1-intgrds2));
    
    if nargout == 2
        g = zeros(nnMinim,1);
        for i=1:nnMinim
            j = nodesMinim(i);
            elems = elems_j{i};
            e2n = elems2nodes_j{i};
            v_elems = v(e2n);
            v_x_elems = sum(Grad_elems_j{i,1}.*v_elems,2);
            v_y_elems = sum(Grad_elems_j{i,2}.*v_elems,2);
            
            D_DvDx1_elems = sum((e2n==j).*Grad_elems_j{i,1},2);
            D_DvDx2_elems = sum((e2n==j).*Grad_elems_j{i,2},2);
            
            intgrds1 = ( D_DvDx1_elems.*sign(v_x_elems).*abs(v_x_elems).^(p-1) + ...
                                          D_DvDx2_elems.*sign(v_y_elems).*abs(v_y_elems).^(p-1) );
            g(i) = sum(areas(elems).*intgrds1);
        end
        g=g-Dintgrds2;
    end
    
end

end




