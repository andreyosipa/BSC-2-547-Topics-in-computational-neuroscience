function [ spikes ] = bernoulli_spikes( time_limit, sp_rate )
delta = 0.001;
N = time_limit/delta;
P = sp_rate * delta;
spikes = binornd(1, P, 1, N);
end

