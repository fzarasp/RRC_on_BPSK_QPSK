function [rand_signal] = BB_signal(M , N)
rand_signal = randi(M , [1 , N]) - 1;
end