Nr = [32 64 128];
Rmax = 30;
betas = [1.6 2.4 3.2];
for n=Nr 
    for r=Rmax
        main(r, n, betas)
    end
end
