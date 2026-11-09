function helper_plot_S(S) 

[M,K] = size(S);

figure(55); clf;

for m = 1:M
    subplot(M,1,m);
    plot(S(m,:));
end
