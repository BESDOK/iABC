function out = iABC(fObj, userdata, nFood, nDim, lb, ub, MaxCycle)
% iABC - Improved Artificial Bee Colony Algorithm 
% Besdok, E. 2026 - Vectorized & Optimized Version

%% Başlangıç popülasyonu
Food = rand(nFood, nDim) .* (ub-lb) +lb;
fitFood = feval(fObj,Food,userdata);

% İlk global best
[fBest, idx] = min(fitFood);
gMin = fBest;
gBest = Food(idx, :);


%% Ana döngü
for t = 1:MaxCycle

    %% Colony
    Food0=Food;
    fitFood0=fitFood;


    %% Employed Bees  (Komşulukta yeni bir çözüm üretirler)
    % phase #1 : Employed arılar için multiple-besin kaynağından faydalanma süreci
    while 1, j=randperm(nFood); if sum(j==1:nFood)==0, break; end, end
    [M,scale] = generateMS(nFood, nDim);

    Employed1 = Food + M .* scale .* (Food(j,:) - Food); 
    [Food, fitFood] = Update(Food, fitFood, Employed1, lb, ub, fObj, userdata);
    
    %% Onlooker Bees (İyi çözümlere olasılıksal olarak yönelirler)
    % phase #2 : Employed arılar tek-besin kaynağından faydalanma süreci
    while 1, j0 = randperm(nFood); if ( sum(j0 == 1 : nFood ) == 0 ), break; end, end
    while 1, j1 = randperm(nFood); if ( sum(j1 == j0 ) == 0 ), break; end, end
    if rand<rand, j0=1:nFood; end


    [~,scale] = generateMS(nFood, nDim);
    if rand<rand, M = generateM1(nFood, nDim); end               
    
    % İyi çözümlere olasılıksal olarak yönelim süreci
    k=SelectByProbability(fitFood0,ceil(nFood/2));
    j3=k(randi(numel(k),1,nFood)); % stochastic phase : goto one of the top-best nectar source
    if rand<rand, j1=j3; end
 
    Employed2 = Food0 + M .* scale .* (Food0(j1,:) - Food0(j0,:) );
    [Food0, f1] = Update(Food0, fitFood0, Employed2, lb, ub, fObj, userdata);

    % Update #1 : Beslenme kaynağı kalite değerlerini koloniye taşıma süreci
    j=f1<fitFood;
    Food(j,:)=Food0(j,:);
    fitFood(j)=f1(j);

    % Update #2 : Küresel çözümleri güncellem süreci
    [fBest, idx] = min(fitFood);
    if fBest < gMin, gMin = fBest; gBest = Food(idx,:); end


    %% Scout Bees (Rastgele yeni çözüm üretirler)
    if rand < 0.10
        newFood = lb + (ub -lb) .* rand(nFood, nDim);
        [~,scale] = generateMS(nFood, nDim);

        % if rand<rand, M = generateM1(nFood, nDim);  end
        M = generateM1(nFood, nDim);

        Scout = Food + M .* scale .* (newFood - Food);
        [Food, fitFood] = Update(Food, fitFood, Scout, lb, ub, fObj, userdata);

        % Global best güncelle
        [fBest, idx] = min(fitFood);
        if fBest < gMin, gMin = fBest; gBest = Food(idx,:); end
    end

    %% Report to screen and workspace
    out.gmin = gMin;
    out.gbest = gBest;

    assignin('base', 'iABCout', out);

    % if mod(t, 100) == 0
    fprintf('Iter=%d -->  fMin=%5.16f\n', t, gMin);
    % end
end
end

%% ========================================================================
%                           ALT FONKSİYONLAR
%% ========================================================================

function [X, fitX] = Update(X, fitX, Y, lb, ub, fObj, userdata)
Y = max(lb, min(ub, Y));
fitY = feval(fObj, Y, userdata);
ind = fitY < fitX;
X(ind,:) = Y(ind,:);
fitX(ind) = fitY(ind);
end




function [M, scale] = generateMS(nFood, nDim)
M = generateMaskFast(nFood, nDim);
while 1
    alpha = randi([-5 5], nFood, 1);
    if alpha~=0, break, end
end
if rand < 0.10
    c = 1;
else
    c = nDim;
end
scale = randn(nFood, c) .^ alpha;
if rand<rand, scale=1./scale; end
end




function M = generateMaskFast(A, B)
M=zeros(A,B);
for i=1:A
    v=randperm(B);
    h=ceil(rand*B);
    M(i,v(1:h))=1;
end
end




function M = generateM1(nFood, nDim)
M = zeros(nFood, nDim);
cols = randi(nDim,nFood,1);                 % her satır için rastgele sütun
M((1:nFood)' + (cols-1)*nFood) = 1;         % sub2ind'in açık hali (garantili)
M=M==1;
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



