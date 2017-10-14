function [ spikes ] = exp_isi( time_limit, sp_rate )
% genrate spikes with bernoulli distribution with
% given spiking rate for given time interval.

lambda = 1/sp_rate; % mean ISI
time = 0;
delta = 0.001;
spikes = zeros(1, length(0:delta:time_limit));
while time < time_limit
    isi = exprnd(lambda,1,1); % generate next ISI
    time = time + isi;
    if time < time_limit
        spikes(ceil(time/delta)) = 1; % fill appropriate bin
    end
end
end

