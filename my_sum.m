A = rand(150,150);
B = rand(150,150);
f = @() sum(A.'.*B, 1);
timeit(f)
