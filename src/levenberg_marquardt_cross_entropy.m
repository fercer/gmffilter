function x_opt = levenberg_marquardt_cross_entropy(x_ini, y_targets, fun, dfun, x_min, x_max, max_iters, tolerance, delta_trust, nu_val)
% LEVENBERG_MARQUARTD searches for the optimal solution of the optimization
% problem in fun.
%
% Author: Fernando Cervantes Sanchez
% Date: June 2019
% e-mail: iie.fercer@gmail.com


    n_vars = length(x_ini);
    if nargin < 10
        nu_val = 0.2;
    end
        
    if nargin < 9
        compute_delta_trust = true;
    else
        compute_delta_trust = false;
    end
    
    if nargin < 8
        tolerance = 1e-6;
    end
    
    if nargin < 7
        max_iters = 50;
    end

    if nargin < 6
        x_max = Inf*ones(n_vars, 1);
        upper_unbounded = true;
    else
        upper_unbounded = false;
    end
    
    if nargin < 5
        x_min = -Inf*ones(n_vars, 1);
        lower_unbounded = true;
    else
        lower_unbounded = false;
    end
    
    % Initialize a trust region max region
    if compute_delta_trust
        if ~(upper_unbounded && lower_unbounded)
            delta_max = min([x_max - x_ini; x_ini - x_min]);
        elseif upper_unbounded
            delta_max = min(x_max - x_ini);
        elseif lower_unbounded
            delta_max = min(x_ini - x_min);
        else
            delta_max = 5*min(x_ini);
        end
        
        delta_trust = 0.25 * delta_max;
    end
    
    x_opt = main_loop(x_ini, y_targets, fun, dfun, x_min, x_max, max_iters, tolerance, delta_trust, nu_val);    
end


function x_opt = main_loop(x_ini, y_targets, fun, dfun, x_min, x_max, max_iters, tolerance, delta_trust, nu_val)
% Levenberg-Marquardt algorithm main loop
    n_vars = length(x_ini);
    iters = 1;
    x_opt = x_ini;
    fx = feval(fun, x_ini);
    residuals = sqrt(-y_targets.*log(fx+1e-12) - (1-y_targets).*log(1 - fx + 1e-12));
    Loss = 0.5*(residuals'*residuals);
    exit_flag = false;
    while iters <= max_iters && ~exit_flag
        delta_trust = min([x_opt - x_min; x_max - x_opt; delta_trust]);
        
        % Compute the Jacobian as the derivative of the residuals function
        dfx = feval(dfun, x_opt);
        J = bsxfun(@rdivide, bsxfun(@times, -y_targets./(fx+1e-12), dfx) + bsxfun(@times, (1-y_targets)./(1-fx+1e-12), dfx), residuals + 1e-12);
        JJ = J'*J;
        Jr = J'*residuals;

        p = -JJ\Jr;
        norm_p = norm(p);
        
        % Subiteration to find a betters solution reducing the trust region
        % radius delta_trust
        for sub_iters = 1:10            
            p_test = p;
            if norm_p >= delta_trust
                lambda = compute_lambda(Jr, JJ, delta_trust);
                H = JJ + lambda*eye(n_vars);
                p_test = -H\Jr;
                norm_p = norm(p_test);
            end
            x_test = x_opt + p_test;
                        
            
            % Test if the vector p satsified the descend conditions
            fx_test = feval(fun, x_test);
            residuals_test = sqrt(-y_targets.*log(fx_test+1e-12) - (1-y_targets).*log(1 - fx_test + 1e-12));
            Loss_test = 0.5*(residuals_test'*residuals_test);

            m_sub_p = Loss + p_test'*Jr + 0.5 * p_test'*JJ*p_test;

            if abs(Loss - m_sub_p) < tolerance
                %display('The subproblem cannot be minimized further than it is now ... (Stopping)')
                exit_flag = true;
                break
            end

            rho_k = (Loss - Loss_test)/(Loss - m_sub_p);

            % Update the delta_trust
            if rho_k < 0.25
                delta_trust = 0.25 * norm_p;
            else
                if rho_k > 0.75 && abs(norm_p - delta_trust) < tolerance
                    delta_trust_max = min([x_opt - x_min; x_max - x_opt]);
                    delta_trust = min([delta_trust_max, 2*delta_trust]);
                end
            end

            if rho_k > nu_val
                Loss = Loss_test;
                residuals = residuals_test;
                x_opt = x_test
                break
            end
        end
        
        display(['Iteration ', num2str(iters), ', Loss: ', num2str(Loss), ', ||Jr||: ', num2str(norm(Jr))]);
        iters = iters + 1;        
    end
end


function lambda = compute_lambda(Jr, JJ, delta_trust)
% Computes a value of lambda that satisfies the constrain of p having a
% norm lower or equal to delta_trust
    delta_trust_2 = delta_trust^2;
    [Q, L, ~] = svd(JJ);
    Q_Jr_2 = (Q'*Jr).^2;
        
    p_norm_lambda = @(lambda) sum(bsxfun(@rdivide, Q_Jr_2, bsxfun(@plus, diag(L), lambda).^2)) - delta_trust_2;
    lambda_search_limits = [-min(diag(L))+1e-16, max(diag(L))];
    
    while true        
        p_norm_left = p_norm_lambda(lambda_search_limits(1));
        
        if ~isfinite(p_norm_left)
            lambda_search_limits(1) = lambda_search_limits(1) + 1e-6;
            continue
        end
        
        p_norm_right = p_norm_lambda(lambda_search_limits(2));
        if sign(p_norm_left)*sign(p_norm_right) < 0
            break
        else
            lambda_search_limits(2) = lambda_search_limits(2) * 10;
        end
    end
    
    lambda = fzero(p_norm_lambda, lambda_search_limits);
end

