function [noise] = noise_generator(n , N0)
noise1 = (randn(1 , n) +1i * (randn(1,n)));
noise = sqrt(N0/2) * noise1 / sqrt(var(noise1));
end