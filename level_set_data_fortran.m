clear all
close all
clc

axis([0 1 0 1 0 1]);

max_dim = input('maximum dimension is : ');
size = 2^max_dim;
for j=0:1:50
    data4 = load(['tree_3dlevel' num2str(j) '.txt']);
    data4 = reshape(data4,size+1,size+1,size+1);
    
    l = linspace(0,1,size+1);
    lp = linspace(1,0,size+1);
    
    plot3([1 2],[1 2],[1,2]);
    hold on;
    figure(1), isosurface(l,lp,l,data4, 0), axis([0 1 0 1 0 1]);
    hold off;

    drawnow();
end