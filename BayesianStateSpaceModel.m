function [samples, stats, structArray] = BayesianStateSpaceModel(aimingError, rotation, showCursor)
    NChains = 1;
    NBurnin = 2*10^4;
    NSamples = 5*10^4;
    doparallel = 1;
    TrialStart = 1;
    TrialEnd = 900;

    subjects = 1:size(aimingError,2);
    rotation = rotation(TrialStart:TrialEnd,subjects);
    showCursor = showCursor(TrialStart:TrialEnd,subjects);
    aimingError = aimingError(TrialStart:TrialEnd,subjects);

    indices = abs(aimingError+rotation)>30;
    aimingError(indices) = NaN;

    dataStruct = struct('y', aimingError, ...
        'p', rotation,...
        'v', showCursor,...
        'NSubjects', size(aimingError,2),...
        'TrialStart', TrialStart,...
        'TrialEnd', TrialEnd);

    for i=1:NChains
        S(i).A1mu = 0.5;
        S(i).A1prec = 0.5;
        S(i).A1 = ones(size(aimingError,2),1);
        S(i).B1mu = 0.5;
        S(i).B1prec = 0.5;
        S(i).B1 = ones(size(aimingError,2),1);
        S(i).x = [aimingError; aimingError(900,:)];
        S(i).etaprec = ones(size(aimingError,2),1);
        S(i).epsilonprec = ones(size(aimingError,2),1);
    end;  

    [samples, stats, structArray] = matjags( ...
    dataStruct, ...                     % Observed data   
    fullfile(pwd, 'BayesianStateSpaceFitReachingTotal.txt'), ...    % File that contains model definition
    S, ...                          % Initial values for latent variables
    'doparallel' , doparallel, ...      % Parallelization flag
    'nchains', NChains,...              % Number of MCMC chains
    'nburnin', NBurnin,...              % Number of burnin steps
    'nsamples', NSamples, ...           % Number of samples to extract
    'thin', 1, ...                      % Thinning parameter
    'dic', 0, ...                       % Do the DIC?
    'monitorparams', {'A','B','q','r'}, ...   % List of latent variables to monitor
    'savejagsoutput' , 1 , ...          % Save command line output produced by JAGS?
    'verbosity' , 2 , ...               % 0=do not produce any output; 1=minimal text output; 2=maximum text output
    'cleanup' , 1);                    % clean up of temporary files?    