N = 2.^(5:8);
Rmax = [10 30 100 300];
betas = 0.1:0.3:6;
for n=N 
    for r=Rmax
        main(r, n+1, n+1, betas);
    end
end
