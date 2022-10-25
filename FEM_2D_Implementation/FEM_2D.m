%{
Copyright (c) 2021, Alexej Moskovka
Copyright (c) 2016, Jan Valdman
Copyright (c) 2015, Jan Valdman
Copyright (c) 2014, Immanuel Anjam
Copyright (c) 2012, Immanuel Anjam
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
* Neither the name of  nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
* Neither the name of University of South Bohemia &  Institute of Information Theory and Automation, Czech Republic nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
* Neither the name of University of Jyväskylä nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%}

% This is a modified version by Andres Rubiano


clear variables
close all
clc

add_paths

%=======================================================
all_levels = 5:5;                   % problem size 
all_methods_on = [1 1 0 0];       % options
RUN = kron(all_methods_on,ones(numel(all_levels),1));
RUN(4:end,3) = 0;                 % exclude some evaluations for option 3
RUN(5:end,4) = 0;                 % exclude some evaluations for option 4  

difference_type = 1;              % 1 ... central, 2 ... forward, 3 ... backward
eps = 1e-5;                       % difference parameter
 
mesh_plot = true;                % plot mesh
graphs = true;                    % graphs of solutions
   
% setup for pLaplacian
approx_solution=@approxsolution_pLaplacian_2D;
params.f_constant = -1;                 % right-hand-side
params.p = 2;                            % power in p-Laplacian
params.pConj = params.p/(params.p-1); 
fprintf('power= %f, conjugate power = %f: \n \n',params.p,params.pConj);
nbfn=1;                                 % number of basis functions in one node 
%=======================================================

all_times = zeros(length(all_levels),4);        % running times of all options
all_iters = zeros(length(all_levels),4);        % number of iters of all options

params.difference_type = difference_type;
params.eps = eps;

u_approx_C = cell(length(all_levels),4);
Ju_approx = zeros(length(all_levels),4);    
all_ndofs = zeros(length(all_levels),1);

separate1 = '=========================================================================';

for le=1:length(all_levels)
    
    level = all_levels(le);
    mesh = mesh_construct(level,mesh_plot);
    
    ndofs = numel(nbfn*mesh.nodesInternal);       % number of degrees of freedom 
    all_ndofs(le)=ndofs;
     
    nn=size(mesh.nodes2coord,1);                  % number of all nodes
    %u_init = ones(nbfn*nn,1);
    u_init = zeros(nbfn*nn,1);
    
    fprintf(separate1); fprintf('\n');
    
    for m=1:4
        if all_methods_on(m) && RUN(le,m)
            fprintf(separate1); fprintf('\n\n');
            switch m
                case 1
                    fprintf('nn = %d,  OPTION 1 : exact gradient, known hessian pattern \n',ndofs);
                case 2
                    fprintf('nn = %d,  OPTION 2 : approx. gradient, known hessian pattern \n',ndofs);
                case 3
                    fprintf('nn = %d,  OPTION 3 : approx. gradient, no hessian pattern \n',ndofs);
                case 4
                    fprintf('nn = %d,  OPTION 4 : no gradient, no hessian pattern \n',ndofs);    
            end
            if m==1
                exanum = 1;
            else
                exanum = 2;
            end
            [u_approx_C{le,m}, Ju_approx(le,m), time, iters] = approx_solution(exanum,mesh,u_init,circshift([0 0 0 1],m),params);
            all_times(le,m) = time;
            all_iters(le,m) = iters;

            if graphs
               figure; 
               show_nodal_scalar_frame(u_approx_C{le,m},mesh.nodes2coord,mesh.elems2nodes,m); axis square; view(3);
               %subplot(1,2,2); tricontour(mesh.nodes2coord,mesh.elems2nodes,u_approx_C{le,m},8); axis square;
               switch m
                   case 1
                       title('')
                   case 2
                       title('Contours')
                   case 3
                       title('Contours')
                   case 4
                       title('Contours')
                end
            end
        end
    end          
end
report_console(all_levels,all_ndofs,all_times,all_iters,RUN);

 if all_methods_on(1) && all_methods_on(2)
     n = length(all_levels);
     err = zeros(2,n);
     for i=1:n
         err(1,i) = max(abs(u_approx_C{i,1} - u_approx_C{i,2}));
         err(2,i) = sum(abs(u_approx_C{i,1} - u_approx_C{i,2}))/length(u_approx_C{i,1});
     end
     format long
     err
     format short
 end

generate_table_paper_latex(all_levels,all_ndofs,all_times,all_iters,RUN)

figure(10)
X = mesh.nodes2coord(:,1);
Y = mesh.nodes2coord(:,2);
Z = zeros(size(X));
one_vec =  ones(size(Y));
for n=1:200
    Z = Z + 2*(4*(sech((pi*n)/2))*(sin((pi*n)/2)^2)*(sin(X*(n*pi))).*(sinh(((Y-one_vec)*(n*pi))./2)).*(sinh((Y*(n*pi))./2)))./(pi*n)^3;;
end

h=trisurf(mesh.elems2nodes,X,Y,Z,'FaceColor','flat','LineWidth',0.1,'EdgeColor','w');

norm(u_approx_C{1,1}-Z)

