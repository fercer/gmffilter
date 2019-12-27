/***********************************************************************************************************************************
CENTRO DE INVESTIGACION EN MATEMATICAS
DOCTORADO EN CIENCIAS DE LA COMPUTACION
FERNANDO CERVANTES SANCHEZ

FILE NAME : ec_algorithm.h

PURPOSE : Declares common functions for Evolutionary Computation (EC) algorithms

FILE REFERENCES :
Name        I / O        Description
None----       ----------

ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES :
None
************************************************************************************************************************************/
function varargout = de(ec_params)
% DE defines the functions required to optimize a numerical optimization
% problem using the Differential Evolution algorithm
    if ~isfield(ec_params, 'mu')
        mu = 4+floor(3*log(ec_params.N));  % population size, offspring number
        ec_params.mu = mu;
    end
    
    if ~isfield(ec_params, 'crossover_rate')
        ec_params.crossover_rate = 0.9;
    end
    
    if ~isfield(ec_params, 'mutation_factor')
        ec_params.mutation_factor = 0.9;
    end
        
    if ~isfield(ec_params, 'stopfitness')
        ec_params.stopfitness = 1e-10;  % stop if fitness < stopfitness (minimization)
    end
    
    if ~isfield(ec_params, 'stopeval')
        ec_params.stopeval = 1e3*ec_params.N^2;   % stop after stopeval number of function evaluations
    end
    
    if ~isfield(ec_params, 'display_info')
        ec_params.display_info = 0;
    end
    
    ec_params.counteval = 0;        
    ec_params.generation = 0;  
    
    %% Main loop:
    
    ec_params = initializePopulation(ec_params);
    
    while ~checkStoppingCriteria(ec_params)
        ec_params = mutateIndividuals(ec_params);
        ec_params = crossoverIndividuals(ec_params);
        ec_params = selectParents(ec_params);
        ec_params.generation = ec_params.generation + 1;
    end
    
    varargout{1} = ec_params.best_fitness;
    
    if nargout >= 2
      varargout{2} = ec_params.elite;
    end
    
    if nargout >= 3
      varargout{3} = ec_params;
    end
end


%% Generate the initial population:
function ec_params = initializePopulation(ec_params)
    ec_params.population = zeros(ec_params.N, ec_params.mu);
    ec_params.fitness = zeros(ec_params.mu, 1);
    ec_params.best_fitness = inf;

    for k=1:ec_params.mu,
        ec_params.population(:, k) = rand(ec_params.N,1) .* (ec_params.upper_bound - ec_params.lower_bound) + ec_params.lower_bound;

        % Evaluate the new individial:       
        ec_params.fitness(k) = feval(ec_params.evaluationFcn, ec_params.population(:, k), ec_params);
        ec_params.counteval = ec_params.counteval+1;
                
        if ec_params.fitness(k) < ec_params.best_fitness
            ec_params.best_fitness = ec_params.fitness(k);
            ec_params.elite = ec_params.population(:,k);
        end
    end
    
    % additional arrays required by the algorithm:
    ec_params.mutation_population = zeros(ec_params.N, ec_params.mu);
    ec_params.crossover_population = zeros(ec_params.N, ec_params.mu);
    
    ec_params.crossover_fitness = zeros(ec_params.mu, 1);
end


%% Mutation process
% x_i = x_r1 + F * (x_r2 - x_r3)
function ec_params = mutateIndividuals(ec_params)
    for k=1:ec_params.mu
        r_indices = randsample(setdiff(1:ec_params.mu, k),3);
        ec_params.mutation_population(:, k) = ec_params.population(:, r_indices(1)) + ec_params.mutation_factor * ...
            (ec_params.population(:, r_indices(2)) - ec_params.population(:, r_indices(3)));        
    end
        
    % Verify that population is within boundaries    
    while 1
      d2u = bsxfun(@minus, ec_params.mutation_population, ec_params.upper_bound);
      ubi = (d2u > 0.0);
      d2u = bsxfun(@minus, ec_params.upper_bound, d2u).* ubi;
    
      d2l = bsxfun(@minus, ec_params.lower_bound, ec_params.mutation_population);
      lbi = (d2l > 0.0);
      d2l = bsxfun(@plus, d2l, ec_params.lower_bound).* lbi;
      
      ec_params.mutation_population = and(~ubi, ~lbi).*ec_params.mutation_population + d2u + d2l;
      if (sum(ubi) + sum(lbi)) < 1
        break
      end
    end

end


%% Crossover process
function ec_params = crossoverIndividuals(ec_params)
    for k = 1:ec_params.mu
        ec_params.crossover_population(:, k) = ec_params.population(:, k);
        
        crossover_mutation_idx = union(find(rand(ec_params.N, 1) <= ec_params.crossover_rate), randi(ec_params.N));
        ec_params.crossover_population(crossover_mutation_idx, k) = ec_params.mutation_population(crossover_mutation_idx, k);

        ec_params.crossover_fitness(k) = feval(ec_params.evaluationFcn, ec_params.crossover_population(:, k), ec_params);

        ec_params.counteval = ec_params.counteval+1;
    end
end


%% Selection process
function ec_params = selectParents(ec_params)
    selection_idx = find(ec_params.crossover_fitness <= ec_params.fitness);
    ec_params.population(:, selection_idx) = ec_params.crossover_population(:, selection_idx);
    ec_params.fitness(selection_idx) = ec_params.crossover_fitness(selection_idx);
    
    [generation_best_fitness, best_individual_idx] = min(ec_params.fitness);
    if generation_best_fitness < ec_params.best_fitness
        ec_params.best_fitness = generation_best_fitness;
        ec_params.elite = ec_params.population(:, best_individual_idx);
    end
end


%% Verify stopping criteria:
function stop_computing = checkStoppingCriteria(ec_params)
    stop_computing = 0;
    
    % Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable 
    if ec_params.best_fitness <= ec_params.stopfitness
      stop_computing = 1;
    end

    if ec_params.counteval > ec_params.stopeval
        stop_computing = 1;
    end

    if ec_params.generation >= ec_params.max_generations
        stop_computing = 1;
    end
    
    % Output 
    if ec_params.display_info
        disp(['[' num2str(ec_params.generation) '] ' num2str(ec_params.counteval) ': ' ...
          num2str(min(ec_params.fitness)) ' (best:' num2str(ec_params.best_fitness) ', worst: ' num2str(max(ec_params.fitness)) ')']);
    end
end
