time_bound = 30000;
delta = 0.001;
sp_rate = 20;
b_spikes = bernoulli_spikes(time_bound, sp_rate);
e_spikes = exp_isi(time_bound, sp_rate);
b_sp_counts = zeros(1, time_bound);
e_sp_counts = zeros(1, time_bound);

disp(length(b_spikes));
disp(length(e_spikes));

for i=1:time_bound
   for idx = (1+(i-1)/delta) : (i/delta)
       if b_spikes(idx) == 1
           b_sp_counts(i) = b_sp_counts(i) + 1;
       end
       if e_spikes(idx) == 1
           e_sp_counts(i) = e_sp_counts(i) + 1;
       end
   end
end

p_sp_counts = poissrnd(sp_rate, 1, time_bound);

figure
hold on
ksdensity(p_sp_counts);
ksdensity(b_sp_counts);
ksdensity(e_sp_counts);
legend("Poisson distribution", "Binomial spikes", "Exponential ISI");
hold off

figure
hold on
histogram(p_sp_counts, 'FaceAlpha', 0.5);
histogram(b_sp_counts,'FaceAlpha', 0.5);
histogram(e_sp_counts, 'FaceAlpha', 0.5);
legend("Poisson distribution", "Binomial spikes", "Exponential ISI");
hold off