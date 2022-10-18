function [u, Ju, time, iters] = start_minimizers(fun_g_approx,Hstr,u_init,all_methods_on)

disp='iter';
% disp='final';
% disp='off';

gradient='GradObj'; yes='on';

if all_methods_on(1)
    % option 1 - exact gradient, known hessian pattern
    options1 = optimoptions(@fminunc,'Algorithm','trust-region',gradient,yes,'HessPattern',Hstr,'display',disp);
    tic
    [u, Ju, ~, output] = fminunc(fun_g_approx,u_init,options1);        
    time = toc;
    iters = output.iterations;
end

if all_methods_on(2)
    % option 2 - approx. gradient, known hessian pattern
    options2 = optimoptions(@fminunc,'Algorithm','trust-region',gradient,yes,'HessPattern',Hstr,'display',disp);
    tic
    [u, Ju, ~, output] = fminunc(fun_g_approx,u_init,options2);
    time = toc;
    iters = output.iterations;
end

if all_methods_on(3)
    % option 3 - approx. gradient, no hessian pattern
    options3 = optimoptions(@fminunc,'Algorithm','trust-region',gradient,yes,'display',disp);
    tic      
    [u, Ju, ~, output] = fminunc(fun_g_approx,u_init,options3);
    time = toc;
    iters = output.iterations;
end

if all_methods_on(4)
    n = numel(u_init);
    % options4 = optimoptions(@fminunc,'display',disp,'MaxFunctionEvaluations',n*1e5,'MaxIterations',n*10);
    options4 = optimoptions(@fminunc,'display',disp,'MaxFunEvals',n*1e5,'MaxIter',n*10);
    tic
    [u, Ju, ~, output] = fminunc(fun_g_approx,u_init,options4);
    time = toc;
    iters = output.iterations;
end

end

