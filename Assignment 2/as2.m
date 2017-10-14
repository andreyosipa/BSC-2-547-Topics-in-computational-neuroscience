%Discovering RF
% given parameters
sx = 1;
sy = 2;
k = 1/0.56/180*pi;
phi = pi/2;

nonlinearity = @(x) max(0,x)^2;

range = -5:0.2:5;
[X,Y] = meshgrid(range);
rf_values = zeros(length(range));
for ix=1:length(range)
   for iy=1:length(range)
      rf_values(ix,iy) = RF( X(ix, iy), Y(ix, iy), sx, sy, phi, k);
   end
end
figure
contour3(X,Y,rf_values,20);
title("RF");



%Simple cell simulation

image_fun = @(x,y,alpha) sin(k.*x./180.*pi - alpha);

% number of different phases to simulate
n_alphas = 200;
out_rates = zeros(1,n_alphas);
alphas = linspace(-pi,pi,n_alphas);
%constant, fitted to have maen spiking rate of 50Hz
adj_c_s = 50/4.6778e-07;
for alpha_idx = 1:n_alphas
    image = zeros(length(range));
    for ix=1:length(range)
        for iy=1:length(range)
            image(ix,iy) = image_fun( X(ix, iy), Y(ix, iy), alphas(alpha_idx));
        end
    end
    out_rates(alpha_idx) = simple_cell(image, range, phi, nonlinearity, adj_c_s);
end

figure
hold on
plot(alphas, out_rates);
set(gca,'fontsize',18);
title("Spiking rate for simple cell");
xlabel("alpha, rad");
ylabel("spike rate");
hold off



%Complex cell simulation

n_alphas = 200;
out_rates = zeros(1,n_alphas);
alphas = linspace(-pi,pi,n_alphas);
for alpha_idx = 1:n_alphas
    image = zeros(length(range));
    for ix=1:length(range)
        for iy=1:length(range)
            image(ix,iy) = image_fun( X(ix, iy), Y(ix, iy), alphas(alpha_idx));
        end
    end
    out_rates(alpha_idx) = complex_cell(image, range);
end

figure
hold on
plot(alphas, out_rates);
ylim([0,60]);
set(gca,'fontsize',18);
xlabel("alpha, rad");
ylabel("spike rate");
title("Spiking rate for complex cell");
hold off



%Tuning curve for complex cell

image_fun_r = @(x,y,alpha,beta) sin(k.*((1-beta).*x./180.*pi + beta.*y/180*pi) - alpha);

n_betas = 200;
out_rates = zeros(1,n_betas);
% grating rotations
betas = [linspace(0,1,n_betas/2), linspace(1,0,n_betas/2)];
for beta_idx = 1:n_betas
    image = zeros(length(range));
    for ix=1:length(range)
        for iy=1:length(range)
            image(ix,iy) = image_fun_r( X(ix, iy), Y(ix, iy), 0, betas(beta_idx));
        end
    end
    out_rates(beta_idx) = complex_cell(image, range);
end

figure
hold on
plot(out_rates);
set(gca,'XTick',linspace(1, n_betas,21));
set(gca, 'XTickLabel', [linspace(0,1,11), linspace(0.9,0,10)]);
set(gca,'fontsize',16);
xlabel("beta");
ylabel("spike rate");
title("Tuning curve for complex cell");
hold off



%Reverse corellation for simple cell

n_images = 10000;
images = normrnd(0,1,51,51,n_images);
sp_counts_simple = zeros(1,n_images);
%constant to make spikes "noticable, despite smnall values.
adj_const = 1e+6;
for idx=1:n_images
    sp_counts_simple(idx) = poissrnd(adj_const * simple_cell(images(:,:,idx), range, phi, nonlinearity, 1));
end

figure
sp_counts_simple = round(0.2 .* sp_counts_simple ./ mean(sp_counts_simple));
disp("mean spiking rate SC: " + mean(sp_counts_simple));
plot(sp_counts_simple);
title("Responses to the noise, Simple Cell");

%compute avrage image
average = zeros(51);
count = 0;
for idx=1:n_images
   if sp_counts_simple(idx) >= 1 
      average = average + images(:,:,idx); 
      count = count + 1;
   end
end

average = average ./ count;
figure
surf(average);
title("RF estimation, Simple Cell");

disp("Correlation matrix with true RF");
disp(corrcoef(reshape(average, 1, 51^2), reshape(rf_values, 1, 51^2)));



%Reverse corellation for complex cell

n_images = 10000;
images = normrnd(0,3,51,51,n_images);
sp_counts_complex = zeros(1,n_images);
for idx=1:n_images
    sp_counts_complex(idx) = poissrnd( 1e-3 * complex_cell(images(:,:,idx), range));
end

figure
sp_counts_complex = round(0.3 .* sp_counts_complex ./ mean(sp_counts_complex));
disp("Mean spiking rate CC: " + mean(sp_counts_complex));
plot(sp_counts_complex);
title("Responses to the noise, Complex Cell");

%compute average
average = zeros(51);
count = 0;
for idx=1:n_images
   if sp_counts_complex(idx) >= 1 
      average = average + images(:,:,idx); 
      count = count + 1;
   end
end
average = average ./ count;
figure
surf(average); 
title("RF estimation, Complex cell");



%STA simple cell
range = -5:1:5;
n_images = 10000;
images = normrnd(0,1,11,11,n_images);
out_rates = zeros(1,n_images);
for idx = 1:n_images
    out_rates(idx) = poissrnd(simple_cell(images(:,:,idx), range, phi, nonlinearity, adj_c_s));
end
out_rates = round(0.2 .* out_rates ./ mean(out_rates));

% find all images that triggered spikes
count = 0;
triggers = [];
for idx=1:n_images
   if out_rates(idx) >= 1 
      count = count + 1;
      triggers = [triggers, reshape(images(:,:,idx), 1, 11^2)];
   end
end

triggers = reshape(triggers, count, 11^2);
covar = cov(triggers);
[V,D] = eig(covar);
figure
scatter(1:11^2, sort(diag(D))); 
title("Covariance matrix eigenvalues (SC)");
figure
surf(reshape(V(:,11^2),11,11))
title("Highest eigenvalue filter (SC)");
figure
surf(reshape(V(:,11^2 -1),11,11))
title("Second highest eigenvalue filter (SC)");



%STA complex cell
range = -5:1:5;
n_images = 10000;
images = normrnd(0,1,11,11,n_images);
out_rates = zeros(1,n_images);
for idx = 1:n_images
    out_rates(idx) = poissrnd(complex_cell(images(:,:,idx), range));
end
out_rates = round(0.2 .* out_rates ./ mean(out_rates));

count = 0;
triggers = [];
for idx=1:n_images
   if out_rates(idx) >= 1 
      count = count + 1;
      triggers = [triggers, reshape(images(:,:,idx), 1, 11^2)];
   end
end

triggers = reshape(triggers, count, 11^2);
covar = cov(triggers);
[V,D] = eig(covar);
figure
scatter(1:11^2, sort(diag(D))); 
title("Covariance matrix eigenvalues (CC)");
figure
surf(reshape(V(:,11^2),11,11))
title("Highest eigenvalue filter (CC)");
figure
surf(reshape(V(:,11^2 -1),11,11))
title("Second highest eigenvalue filter (CC)");