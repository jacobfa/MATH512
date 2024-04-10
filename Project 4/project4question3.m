%Part A
% Parameters
mu = 2;
sigma = 0.1;
X0 = 3;
T = 1;

% Array of different Delta_t values
Delta_ts = linspace(0.001, 0.1, 100); % From 10^-4 to 10^-1, 20 points
%Delta_ts = logspace(-4, -2, 100); % From 10^-4 to 10^-1, 20 points
% Initialize array to hold errors
errors_weak = zeros(size(Delta_ts));

for i = 1:length(Delta_ts)
    Delta_t = Delta_ts(i);
    N = round(T / Delta_t);
    numPaths = 100000;
    finalValues = zeros(numPaths, 1);
    
    for j = 1:numPaths
        W = sqrt(Delta_t) * randn(1, N); % Brownian increments
        X_EM = X0;
        for k = 1:N
            X_EM = X_EM + mu * X_EM * Delta_t+ sigma * X_EM * W(k);
        end
        finalValues(j) = X_EM;
    end
    
    % Expected value of the Euler-Maruyama approximation
    E_X_EM = mean(finalValues);
    
    % Exact solution's expected value
    E_X_exact = X0 * exp(mu * T);
    
    % Weak error estimation
    errors_weak(i) = abs(E_X_EM - E_X_exact);
end

% % Plotting
% Assuming errors_weak and Delta_ts are your data vectors for weak convergence
% Use polyfit to perform linear regression
coeffs_weak = polyfit(Delta_ts, errors_weak, 1);

% Use polyval to evaluate the polynomial coefficients at Delta_ts values for plotting
fit_weak = polyval(coeffs_weak, Delta_ts);

% Calculate R^2 for the weak convergence fit
% Step 1: Total sum of squares (SST)
SST_weak = sum((errors_weak - mean(errors_weak)).^2);
% Step 2: Residual sum of squares (SSR)
SSR_weak = sum((errors_weak - fit_weak).^2);
% Step 3: R^2 calculation
R_squared_weak = 1 - (SSR_weak / SST_weak);

% Display R^2
disp(['R^2 for the weak convergence fit: ', num2str(R_squared_weak)]);

% Plotting for visualization
figure;
plot(Delta_ts, errors_weak, 'o', Delta_ts, fit_weak, '-');
xlabel('\Delta t');
ylabel('Weak Error');
title(['Weak Convergence with Linear Fit, R^2=', num2str(R_squared_weak)]);
legend('Data', 'Linear Fit', 'Location', 'Best');
grid on;



% Part B

% Initialize array to hold strong errors
errors_strong = zeros(size(Delta_ts));

for i = 1:length(Delta_ts)
    Delta_t = Delta_ts(i);
    N = round(T / Delta_t);
    numPaths = 100000;
    errors = zeros(numPaths, 1);
    
    for j = 1:numPaths
        W = sqrt(Delta_t) * randn(1, N); % Brownian increments
        X_EM = X0;
        for k = 1:N
            X_EM = X_EM + mu * X_EM * Delta_t + sigma * X_EM * W(k);
        end
        
        % Exact solution at T
        X_exact = X0 * exp(mu * T + (sigma^2 / 2) * T);
        
        % Strong error for each path
        errors(j) = abs(X_EM - X_exact);
    end
    
    % Average strong error
    errors_strong(i) = mean(errors);
end

% % Plotting
% Assuming errors_strong and Delta_ts are already defined
% Transform Delta_ts to represent Delta t^0.5 for the strong convergence case
Delta_ts_sqrt = sqrt(Delta_ts);

% Linear regression on the transformed Delta_ts and errors_strong
coeffs_strong = polyfit(Delta_ts_sqrt, errors_strong, 1);

% Evaluate the polynomial at the transformed Delta_ts values
fit_strong = polyval(coeffs_strong, Delta_ts_sqrt);

% Calculate R^2
% Step 1: Total sum of squares (SST)
SST = sum((errors_strong - mean(errors_strong)).^2);
% Step 2: Residual sum of squares (SSR)
SSR = sum((errors_strong - fit_strong).^2);
% Step 3: R^2 calculation
R_squared = 1 - (SSR / SST);

% Display R^2
disp(['R^2 for the strong convergence fit: ', num2str(R_squared)]);

% Plotting for visualization
figure;
plot(Delta_ts_sqrt, errors_strong, 'o', Delta_ts_sqrt, fit_strong, '-');
xlabel('\Delta t^{0.5}');
ylabel('Strong Error');
title(['Strong Convergence with Linear Fit, R^2=', num2str(R_squared)]);
legend('Data', 'Linear Fit', 'Location', 'Best');
grid on;

figure;
plot(Delta_ts, errors_strong, 'o');
hold on;
plot(Delta_ts, errors_weak, '^');
hold off;


