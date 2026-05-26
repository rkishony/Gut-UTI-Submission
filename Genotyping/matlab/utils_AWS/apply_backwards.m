function rA = apply_backwards(func, A, d)
rA = flip(func(flip(A,d),d),d);
end
