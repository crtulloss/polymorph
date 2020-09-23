close all
clear all

polymorph_pixel_tasks = [63,59,56,56,56,56,57,57,57,59,57,57,54,60,61,...
    57,58,58,57,56,57,55,53,49,49,49,48,47,47,42,42,42,42,42,42,42,42,...
    42,42,42,42,45,45,45,45,45,44,43,44,44,44,44,44,44,44,49,46,47,47,...
    47,47,46]';
mindfuzz_design_tasks = [36,40,40,44,43,41,39,33,35,35,35,35,35,35,35,...
    35,35,35,35,34,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,...
    33,33,33,33,33,33,33,33,33,33,33,37,33,40,39,36,36,36,36,36,36,36,...
    36,36,36]';
total_design_tasks = polymorph_pixel_tasks + mindfuzz_design_tasks;

l = length(polymorph_pixel_tasks);

x = linspace(0, l-1, l);

polymorph_deltas = polymorph_pixel_tasks(2:l) -...
    polymorph_pixel_tasks(1:l-1);
mindfuzz_deltas = mindfuzz_design_tasks(2:l) -...
    mindfuzz_design_tasks(1:l-1);
total_deltas = polymorph_deltas + mindfuzz_deltas;

figure
plot(x,polymorph_pixel_tasks);
hold on
plot(x,mindfuzz_design_tasks);
plot(x,total_design_tasks);
title('total task lists');
legend('polymorph','mindfuzz','total');

figure
plot(-1*polymorph_deltas);
hold on
plot(-1*mindfuzz_deltas);
plot(-1*total_deltas);
title('daily progress');
legend('polymorph', 'mindfuzz', 'total');