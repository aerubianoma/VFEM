function report_console(all_levels,all_nn,all_times,all_iters,RUN)

separate2 = '|------|----------|----------|----------|----------|';

fprintf('=========================\n');
fprintf('P-LAPLACE IN 2D: RESULTS \n');
fprintf('=========================\n\n');

for k=1:length(all_levels)
    for j=1:4
        if RUN(k,j) == 0
            all_times(k,j) = inf;
            all_iters(k,j) = inf;
        end
    end
end

% TIMES
fprintf('Times\n');
fprintf(separate2); fprintf('\n');
fprintf('|  nn  | option 1 | option 2 | option 3 | option 4 | \n');
fprintf(separate2); fprintf('\n');
for k=1:length(all_levels)
    nn = all_nn(k);
    time1 = all_times(k,1);  time2 = all_times(k,2);  time3 = all_times(k,3);  time4 = all_times(k,4);
    if nn>=10 && nn<=99
        space = "   ";
    elseif nn>=100 && nn<=999
        space = "  ";
    elseif nn>=1000 && nn<=9999
        space = " ";
    else
        space = "";
    end
    fprintf('|');  fprintf(space);
    fprintf('%d |%9.2f |%9.2f |%9.2f |%9.2f | \n',nn,time1,time2,time3,time4);
    fprintf(separate2); fprintf('\n');
end

% ITERATIONS
fprintf('\n');
fprintf('Iterations\n');
fprintf(separate2); fprintf('\n');
fprintf('|  nn  | option 1 | option 2 | option 3 | option 4 | \n');
fprintf(separate2); fprintf('\n');
for k=1:length(all_levels)
    nn = all_nn(k);
    iter1 = all_iters(k,1);  iter2 = all_iters(k,2);  iter3 = all_iters(k,3);  iter4 = all_iters(k,4);
    if nn>=10 && nn<=99
        space = "   ";
    elseif nn>=100 && nn<=999
        space = "  ";
    elseif nn>=1000 && nn<=9999
        space = " ";
    else
        space = "";
    end
    fprintf('|');  fprintf(space);
    fprintf('%d |%9d |%9d |%9d |%9d | \n',nn,iter1,iter2,iter3,iter4);
    fprintf(separate2); fprintf('\n');
end
