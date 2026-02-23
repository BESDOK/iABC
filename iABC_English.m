%{

Improved Artificial Bee Colony Algorithm (iABC)
Besdok, E. 2026, Kayseri, TÜRKİYE


iABC is a variant of the ABC algorithm developed to solve real-valued numerical optimization problems.
iABC is designed as a single-objective, bounded, swarm-based, iterative population evolution-based evolutionary global minimizer algorithm.
iABC utilizes three distinct phases while searching for the global minimum of a given problem: Employed-Bee search, Onlooker-Bee search, and Scout-Bee search phases.

%}


function out = iABC_English(ObjFnc, userdata, N, D, lb, ub, MaxCycle)

% Setting of Initial Population:
% This process initializes the set of Flowers (i.e., FlowerField) from which honeybees can collect nectar. Technically, FlowerField is an
% [N D]-dimensional matrix. Each row of this matrix corresponds to a D-dimensional Flowers vector.
% That is; FlowerField=[Flowers_i] such that i={1,2,3,...,N}, containing N potential nectar sources (i.e., Flowers).
% Similarly; Flowers_i=[NectarSite_j] such that j={1,2,3,...,D}, containing D NectarSites.

% In this regard, a FlowerField matrix consists of Flowers row-vectors carrying possible nectar sources randomly scattered over a hypothetical Meadow.
% In the iABC algorithm, the concept of Meadow analogically corresponds to the 'search space' of the problem at hand.
FlowerField = rand(N, D) .* (ub-lb) +lb; 

% In iABC, the productivity level of each Flowers vector is measured using the objective function (ObjFun) at hand.
% Here, the 'userdata' variable is used to carry additional values that need to be provided to the objective function.
fitFlowers = feval(ObjFnc,FlowerField,userdata);

% Initialization of Global Solution:
% The most productive Flowers vector is referred to as 'the primary nectar source' (i.e., TheBestFlowers).
% The quality level for 'TheBestFlowers' is equal to the value 'fitTheBestFlowers'.
[fitTheBestFlowers, idx] = min(fitFlowers);

% In the iABC algorithm, the Flowers vector with the best quality level currently held by FlowerField is the best nectar source.
% In iABC, the best nectar source corresponds to the global best solution vector (i.e., gBest).
% The quality level of gBest corresponds to the global best solution value (i.e., gBestVal).
gBest = FlowerField(idx, :);
gBestVal = fitTheBestFlowers;

% The iABC algorithm defines the potential nectar sources over which Scout-Bees can fly across the Meadow
% to discover new nectar sources, using Flowers vectors called 'wildflowers'.
% The primary role of 'wildflowers' vectors is to preserve numerical diversity in the FlowerField matrix.
% The initial forms of 'wildflowers' vectors must be generated before the iterative search process is initiated.
wildflowers = lb + (ub -lb) .* rand(N, D); 

%% Main loop: Defines the iterative search process of iABC. This process
% consists of three sub-processes; (1) Employed-Bees, (2) Onlooker-Bees, and (3)
% Scout-Bees based search processes.
for iteration = 1:MaxCycle
    
    % The first step of the iterative search process is to store the information belonging to the current FlowerField matrix for later sharing with Onlooker-Bees.
    Flowers0=FlowerField;
    fitFlowers0=fitFlowers;


    %% Employed-Bee Phase: Analogically corresponds to bijective-bio interaction. This phase analogically models the nectar search of hypothetical bees called EmployedBees, flying over the Meadow toward a Flowers source that is necessarily different from their current Flowers vector.
    % EmployedBees search for potential nectar sources by flying over the Meadow along the direction vectors (i.e., evolutionary directions) defined by dx, starting from their current Flowers.
    % During this new nectar source search process, EmployedBees make equal-amplitude flights from the current Flowers vector toward nectar sources lying on the dx direction vectors, without favoring any particular nectar source.
    % The equal-amplitude flight process is controlled by the variable T (T=0.10).
    while 1, j=randperm(N); if sum(j==1:N)==0, break; end, end
    T=0.10; % This facilitates the solution of separable problems where the relationships between variables are loose.
    % In iABC, which nectar sources will be visited is controlled by the binary-valued M matrix. The amplitude of the corresponding flight (i.e., evolutionary step size) is controlled by the 'scale' variable.
    [M,scale] = generateMS(N, D,T);
    dx = FlowerField(j,:) - FlowerField; % Direction vector for the Employed-Bee phase
    EmployedBee = FlowerField + M .* scale .* dx;   % Morphogenesis process for the Employed-Bee phase
    % Update #1: The information provided by the Flowers sources (i.e., EmployedBee) obtained in the Employed-Bee phase is injected into FlowerField by applying greedy selection rules based on the quality level of the corresponding Flowers sources.
    [FlowerField, fitFlowers] = Update(FlowerField, fitFlowers, EmployedBee, lb, ub, ObjFnc, userdata);
    
    %% Onlooker-Bee Phase: This is a "top-best solutions"-based global search process that favors some nectar sources over others.
    % The purpose of the Onlooker-Bee process is to enable a hypothetical bee to leave its current Flowers and search for nectar sources in Flowers locations containing relatively more productive nectar sources across the Meadow.
    % The Onlooker-Bee phase is a partially elitist search process.
    while 1, j0 = randperm(N); if ( sum(j0 == 1 : N ) == 0 ), break; end, end
    while 1, j1 = randperm(N); if ( sum(j1 == j0 ) == 0 ), break; end, end 
    
    % The T=0.90 setting facilitates the solution of complex problems with hybrid structures among variables.
    T=0.90; % In the Onlooker-Bee phase, some nectar sources can be favored more than others. This process is controlled by the value T=0.90. The T value manages the process of differentiating the amplitudes of flights toward nectar sources.
    [M,scale] = generateMS(N, D, T);    
    
    % In this phase, the index numbers (i.e., ind) of the relatively more productive Flowers vectors are generated using the 'SelectByProbability' function.
    ind=SelectByProbability(fitFlowers0,ceil(N/2));
    
    % To obtain elitist direction vectors, N random selections are made from among the relatively most productive sources, and the indices of the Flowers vectors to fly toward are determined (i.e., indFlowers).
    indFlowers=ind(randi(numel(ind),1,N)); 

    % Whether the Onlooker-Bee phase will behave in an elitist manner in the current iteration is decided using a random mechanism.
    if rand<rand, j1=indFlowers; end
    
    % The dx direction vectors used in the Onlooker-Bee Phase are randomly formed in either a bijective (one-to-one projection) or surjective (onto projection) structure.
    % In elitist processes, the Flowers interactions required for dx occur in surjective form. This provides the opportunity to benefit more from a relatively more productive source.
    dx = Flowers0(j1,:) - FlowerField(j0,:);  % Direction vector for the Onlooker-Bee phase
    OnlookerFlowers = FlowerField + M .* scale .* dx;  % Morphogenesis process for the Onlooker-Bee phase
    
    % Update #2: The information provided by the Flowers vectors (i.e., OnlookerBee) obtained in the Onlooker-Bee phase is injected into FlowerField in two steps
    % by applying greedy selection rules based on the quality level of the corresponding Flowers sources.
    % Since this process enables the FlowerField vectors (i.e., Flowers vectors) stored at the beginning of the iteration to be evolved further in a different way, it contributes to the exploitation (local search) capability of iABC.
    [nectar1, fitnectar1] = Update(Flowers0, fitFlowers0, OnlookerFlowers, lb, ub, ObjFnc, userdata);
    
    % Update #3: In the final step of the Onlooker-Bee phase, the FlowerField information (i.e., the corresponding Flowers vectors and their quality values, fitFlowers) is updated.
    j=fitnectar1<fitFlowers;
    FlowerField(j,:)=nectar1(j,:);
    fitFlowers(j)=fitnectar1(j);
     

    %% Scout-Bee Phase: This is the process that enables iABC to avoid the problem of numerical diversity loss in the FlowerField matrix.
    % Analogically, it implies that the bees called ScoutBees are aware of the existence of wildflowers sources they have not directly reached.
    % ScoutBees determine the direction in which they will fly to obtain nectar by utilizing the wildflowers locations. The wildflowers locations are updated every 100 iteration steps. This enables ScoutBee vectors
    % to generate direction vectors without encountering the numerical diversity problem.
    
    if rand < 0.10        
        [M,scale] = generateMS(N, D,0.90);  % In the Scout-Bee phase, nectar sources are favored at different levels with randomly varying amplitudes.
        % This facilitates the solution of strongly-related problems where there are complex relationships among variables.
        dx = wildflowers - FlowerField; % Direction vector for ScoutBee
        ScoutBee = FlowerField + M .* scale .* dx;  % Morphogenesis process for the Scout-Bee phase
        % Update #4: In the final step of the Scout-Bee phase, the FlowerField information (i.e., the corresponding Flowers vectors and their quality values, fitFlowers) is updated.
        [FlowerField, fitFlowers] = Update(FlowerField, fitFlowers, ScoutBee, lb, ub, ObjFnc, userdata);
    end
    
    % Update #5: This is the phase for updating the global solutions.
    % In the iABC algorithm, the Flowers vector with the best quality level currently held by FlowerField is the best nectar source.
    % In iABC, the best nectar source corresponds to the global best solution vector (i.e., gBest).
    % The quality level of gBest corresponds to the global best solution value (i.e., gBestVal).
    [fBest, idx] = min(fitFlowers);
    if fBest < gBestVal, gBestVal = fBest; gBest = FlowerField(idx,:); end
    
    % Reporting of global solutions to the screen and the MATLAB Workspace
    out.gbest = gBest;   % Global best Flowers (nectar source)
    out.gval = gBestVal; % Quality of the global best Flowers vector
    assignin('base','iABC_out',out); % workspace reporter

    fprintf('Iter=%d -->  fMin=%5.16f\n', iteration, gBestVal);

    if mod(iteration,100)==0, wildflowers = lb + (ub -lb) .* rand(N, D); end % saves diversity of evolutionary-direction
end % end of iteration

end % end of function



%%                           SUB-FUNCTIONS

function [X, fitX] = Update(X, fitX, Y, lb, ub, fObj, userdata)
Y = max(lb, min(ub, Y));
fitY = feval(fObj, Y, userdata);
ind = fitY < fitX;
X(ind,:) = Y(ind,:);
fitX(ind) = fitY(ind);
end

function [M, scale] = generateMS(nnectar, nDim, T)
M = generateMask(nnectar, nDim);
q=[-5:1:-1 1:5]';
alpha=q(randi(10,nnectar,1));
if rand <= T  % The T value determines whether the bee collects nectar from all flowers in a local habitat (c=nDim) or from only one flower (c=1)
    c = 1;
else
    c = nDim;
end
while 1
    % scale: evolutionary step-size. Scale controls the length of the path 
    % bees must travel to reach the nectar source.
    scale = randn(nnectar, c) .^ alpha;
    scale(~isfinite(scale)) = randn(nnz(~isfinite(scale)),1);
    scale(~isfinite(scale)) = randn(nnz(~isfinite(scale)),1);
    if ~any(scale(:)==0) && ~any(scale(:)==1), break; end
end
end


function M = generateMask(A, B)
% M: Controls which flowers at the visited nectar source will have nectar collected from them
M=zeros(A,B);
k=abs(randi([0 1],1) - rand.^randi([2 10]) );
for i=1:A
    v=randperm(B);    
    h=ceil(k*B);
    M(i,v(1:h))=1;
end
end


function ind=SelectByProbability(x,N)
y=abs(x+min(x)+eps);
f=y.^5;
p=f./sum(f);
ind=nan(N,1);
for i=1:N
    [~,ind(i)]=min( abs( p - rand^randi([5 10]) ) );
end
end