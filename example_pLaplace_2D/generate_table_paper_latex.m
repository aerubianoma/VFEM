function generate_table_paper_latex(all_levels,all_nn,all_times,all_iters,RUN)

separate2 = '|------|----------|----------|----------|----------|';

for k=1:length(all_nn)
    for j=1:4
        if RUN(k,j) == 0
            all_times(k,j) = 0;
            all_iters(k,j) = 0;
        end
    end
end

% TIMES
fprintf('\n');
fprintf('Times\n');
fprintf('|  n  | option 1 | option 2 | option 3 | option 4 | \n');
fprintf(separate2); fprintf('\n');
for k=1:length(all_nn)
    nn = all_nn(k);  [nm,ne] = format_own(nn);
    t1 = all_times(k,1);  t2 = all_times(k,2); t3 = all_times(k,3);  t4 = all_times(k,4);
    fprintf(' %2.2fe%d & %9.2f & %9.2f & %9.2f & %9.2f \\\\ \n',nm,ne,t1,t2,t3,t4);
    %fprintf(separate2); fprintf('\n');
end

% ITERATIONS
fprintf('\n');
fprintf('Iterations\n');
fprintf('|  n  | option 1 | option 2 | option 3 | option 4 | \n');
fprintf(separate2); fprintf('\n');
for k=1:length(all_nn)
    nn = all_nn(k);  [nm,ne] = format_own(nn);
    i1 = all_iters(k,1);  i2 = all_iters(k,2); i3 = all_iters(k,3);  i4 = all_iters(k,4);
    fprintf(' %2.2fe%d & %9d & %9d & %9d & %9d \\\\ \n',nm,ne,i1,i2,i3,i4);
    %fprintf(separate2); fprintf('\n');
end

% BOTH TIMES AND ITERATIONS
fprintf('\n');
fprintf('Times + Iterations\n');
fprintf('|  n  | option 1 | option 2 | option 3 | option 4 | \n');
fprintf(separate2); fprintf('\n');
for k=1:length(all_nn)
    nn = all_nn(k);  
    [nm,ne] = format_own(nn);
    t1 = all_times(k,1);  i1 = all_iters(k,1);  
    t2 = all_times(k,2);  i2 = all_iters(k,2);
    t3 = all_times(k,3);  i3 = all_iters(k,3);  
    t4 = all_times(k,4);  i4 = all_iters(k,4);
       
    txt1=string_out(t1,i1);
    txt2=string_out(t2,i2);
    txt3=string_out(t3,i3);
    txt4=string_out(t4,i4);
    
    fprintf(' %d & %s & %s & %s & %s \\\\ \n',nn,txt1,txt2, txt3, txt4);
    %fprintf(separate2); fprintf('\n');
end

function text_out=string_out(t,i)
if i==0
   text_out=sprintf(' - & - ');
else
   text_out=sprintf(' %5.2f & %4d ',t,i); 
end
end

end
