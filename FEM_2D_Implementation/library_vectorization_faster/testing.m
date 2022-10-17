clear all

nx=2; ny=2; nz=1e7;

vv=2*ones(nx,nz);
AAA=ones(nx,ny,nz);
BBB=2*ones(nx,ny,nz);

tic
DDD=avtam(vv,AAA);
toc
size(DDD)

% tic
% q = bsxfun(@times,DDD,DDD); 
% toc
% size(q)
% 
% tic
% bsxfun(@times,AAA(1,:,:),BBB(:,1,:));
% toc
% 
% 
% A3d = rand(200,300,400);
% A1 = A3d; A2 = A3d; A3 = A3d;
% 
% tic
% m = mean(A1,3);
% for i=1:size(A1,3)
%    A1(:,:,i) = A1(:,:,i) - m;
% end
% nonvec = toc
% 
% tic
% A2 = bsxfun(@minus,A2,mean(A2,3));
% vec = toc
% 
% 
% tic  
% bsxfun(@times,ones(2,2,1000),[1 2; 3 4]);
% 
% toc
