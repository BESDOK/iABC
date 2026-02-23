%{

Improved Artificial Bee Colony Algorithm (iABC)
Besdok, E. 2026 , Kayseri, TÜRKİYE


iABC algoritması reel değerli sayısal optimizasyon problemlerini çözmek üzere geliştirilmiş bir ABC algoritması varyantıdır. 
iABC algoritması single-objective, bounded, swarm-based, itarative population evolution tabanlı evrimsel bir global minimizer algorithm olarak tasarlanmıştır. 
iABC, eldeki bir probleme ait küresel minimumu ararken üç farklı aşamadan yararlanır: Employed-Bee search, Onlooker-Bee search, ve Scout-Bee search aşamaları.

%}


function out = iABC_Turkce(ObjFnc, userdata, N, D, lb, ub, MaxCycle)

% Setting of Initial Population:
% Bu süreç bal arılarının nektar toplayabileceği Flowers kümesini (i.e., FlowerField) initialize eder. Teknik olarak FlowerField [N D] boyutlu bir
% matristir. Bu matrisin her bir satırı D-boyutlu bir Flowers vektörüne
% karşılık gelir. Yani; FlowerField=[Flowers_i] öyle ki i={1,2,3,...,N} olmak üzere N adet potansiyel nektar kaynağı (i.e., Flowers) içerir. 
% Benzer şekilde; Flowers_i=[NectarSite_j ] öyle ki j={1,2,3,...,D} olmak üzere D sayıda NectarSite içermektedir.

% Bu açıdan bir FlowersField matrisi, hipotetik bir Meadow üzerine rasgele serpiştirilmiş halde bulunan olası nektar kaynaklarını taşıyan Flowers satır-vektörlerinden oluşmaktadır.
% iABC algoritmasında Meadow kavramı anolojik olarak eldeki probleme ait 'search space' kavramına karşılık gelmektedir.
FlowerField = rand(N, D) .* (ub-lb) +lb; 

% iABC'de, her bir Flowers vektörlerinin sahip olduğu verimlilik seviyesi eldeki objective function (ObjFun) fonksiyonu kullanılarak ölçülür. 
% Burada, 'userdata' değişkeni ilgili objective function'a ilave olarak sağlanması gereken değerleri taşımak için kullanılmaktadır.
fitFlowers = feval(ObjFnc,FlowerField,userdata);

% Initialization of Global Solution:
% En verimli olan Flowers vektörü 'the primary nectar source' (i.e., TheBestFlowers) olarak adlandırılır. 
% 'TheBestFlowers' için kalite seviyesi 'fitTheBestFlowers' değerine eşittir.
[fitTheBestFlowers, idx] = min(fitFlowers);

% iABC algoritmasında, FlowerField'ın bulunulan anda sahip olduğu en iyi kalite seviyesine sahip Flowers vektörü en-iyi nektar kaynağıdır. 
% iABC'de, en-iyi nektar kaynağı küresel en-iyi çözüm vektörüne (i.e., gBest) karşılık gelir. 
% gBest'in  kalite seviyesi ise küresel en-iyi çözüm değerine (i.e., gBestVal) karşılık gelir.
gBest = FlowerField(idx, :);
gBestVal = fitTheBestFlowers;

% iABC algoritması Scout-Bee'lerin yeni nektar kaynakları keşfetmek üzere
% Meadow üzerinde uçabilecekleri potansiyel nektar kaynaklarını
% 'wildflowers' olarak adlandırılan Flowers vektörler kullanarak tanımlamaktadır.
% 'wildflowers' vektörlerinin temel görevi FlowerField matrisindeki sayısal çeşitliliği (diversity) korumaktır. 
% 'wildflowers' vektörlerinin initil formları iteratif arama süreci başlatılmadan üretilmelidir.
wildflowers = lb + (ub -lb) .* rand(N, D); 

%% Main loop: iABC'nin iteratif arama sürecini tanımlar. Bu süreç, üç
% alt-süreçten oluşur; (1) Employed-Bees, (2) Onlooker-Bees, ve (3)
% Scout-Bees tabanlı arama sürecleri.
for iteration = 1:MaxCycle
    
    % İteratif arama sürecinin ilk adımı eldeki FlowerField matrisine ait bilgiyi daha sonra Outlooker-Bee'lerle paylaşmak üzere saklamaktır.
    Flowers0=FlowerField;
    fitFlowers0=fitFlowers;


    %% Employed-Bee Aşaması: Anolojik olarak bijective-bio interaction'a karşılık gelir. Bu aşama, EmployedBee olarak adlandırılan hypotetik arıların o an bulunduğu Flowers vektöründen mutlaka farklı bir
    % Flowers kaynağına doğru Meadow üzerinde uçarak nektar aramasını anolojik olarak modellemektedir. EmployedBee'ler, o an bulundukları Flowers'dan başlayarak dx tarafından tanımlanan doğrultu-vektörleri (i.e., evolutionary directions) boyunca Meadow üzerinde uçarak potansiyel
    % nektar kaynakları ararlar. EmployedBee'ler, bu yeni nektar kaynağı arama süreci boyunca eldeki Flowers vektörününden dx doğrultu vektörleri üzerinde bulunan nektar kaynaklarına doğru eşit genlikte uçuşlar yapar ve herhangi bir nektar kaynağını kayırmazlar. 
    % İlgili eşit genlikte uçuş yapması süreci T değişkeni tarafından kontrol edilmektedir (T=0.10).
    while 1, j=randperm(N); if sum(j==1:N)==0, break; end, end
    T=0.10; % % Bu durum değişkenler arasındaki ilişkilerin gevşek olduğu sperable problemlerin çözümünü kolaylaştırmaktadır.
    % iABC'de hangi nektar kaynaklarına uçuş yapılacağı binary-valued M matrisi tarafından kontrol edilmektedir. İlgili uçuşun genliği  (i.e, evolutionary step size) 'scale' değişkeni tarafından kontrol edilmektedir.
    [M,scale] = generateMS(N, D,T);
    dx = FlowerField(j,:) - FlowerField; % Employed-Bee aşaması için doğrultu vektörü
    EmployedBee = FlowerField + M .* scale .* dx;   % Employed-Bee aşaması için morphogenesis süreci
    % Update #1: Employed-Bee aşamasında elde edilen Flowers kaynaklarının (i.e., EmployedBee) sağladığı bilgi ilgili Flowers kaynaklarının sağladığı kalite seviyesine göre aç-gözlü seçim kuralları işletilerek FlowerField içine enjekte edilir.
    [FlowerField, fitFlowers] = Update(FlowerField, fitFlowers, EmployedBee, lb, ub, ObjFnc, userdata);
    
    %% Onlooker-Bee Aşaması: Bu aşama nektar kaynaklarından bazılarını diğerlerine göre daha çok kayıran "top-best solutions" tabanlı küresel bir arama sürecidir. 
    % Onlooker-Bee sürecinin amacı; bir hipotetik arının o an bulunduğu Flowers'dan ayrılarak görece daha verimli nektar kaynakları içeren Flowers lokasyonlarına doğru Meadow üzerinde nektar kaynağı aramasını sağlamaktır. 
    % Onlooker-Bee aşaması kısmen elitist bir arama sürecidir. 
    while 1, j0 = randperm(N); if ( sum(j0 == 1 : N ) == 0 ), break; end, end
    while 1, j1 = randperm(N); if ( sum(j1 == j0 ) == 0 ), break; end, end 
    
    % T=90 durumu değişkenler arasında hybrid yapılı kompleks problemlerin çözümünü kolaylaştırmaktadır.
    T=0.90; % Onlooker-Bee aşamasında nektar kaynaklarından bazıları diğerlerinden daha fazla kayrılabilir. Bu süreç T=0.90 değeri ile kontrol edilmektedir. T değeri nektar kaynaklarına doğru yapılan uçuşların genliklerini farklılaştırma sürecini yönetmektedir.
    [M,scale] = generateMS(N, D, T);    
    
    % Bu aşamada görece daha verimli olan Flowers vektörlerine ait indis numaraları (i.e., ind) 'SelectByProbability' fonksiyonu kullanılarak üretilmektedir.
    ind=SelectByProbability(fitFlowers0,ceil(N/2));
    
    % Elitist doğrultu vektörleri elde etmek amacıyla görece en verimli kaynaklar arasından N defa rasgele seçim yapılır ve kendisine doğru uçulacak Flowers vektörlerine ait indisler belirlenir (i.e., indFlowers).
    indFlowers=ind(randi(numel(ind),1,N)); 

    % Onlooker-Bee aşamasının bulunulan iterasyonda elitist davranıp davranmayacağına rasgele bir mekanizama kullanılarak karar verilir
    if rand<rand, j1=indFlowers; end
    
    % Onlooker-Bee Aşamasında kullanılan dx doğrultu vektörleri rasgele olarak bijektif (one-to-one projection) yapıda veya surjective (onto projection) yapı oluşur. 
    % Elitist süreçlerde dx için gerekli Flowers interasksiyonları surjective formda oluşur. Bu ise görece daha verimli bir kaynaktan, görece daha çok yararlanma olanağı sağlar.
    dx = Flowers0(j1,:) - FlowerField(j0,:);  % Onlooker-Bee aşaması için doğrultu vektörü
    OnlookerFlowers = FlowerField + M .* scale .* dx;  % Onlooker-Bee aşaması için morphogenesis süreci
    
    % Update #2: Onlooker-Bee aşamasında elde edilen Flowers vektörlerinin (i.e., OnlookerBee) sağladığı bilgi ilgili Flowers kaynaklarının 
    % sağladığı kalite seviyesine göre aç-gözlü seçim kuralları işletilerek iki aşamada FlowerField içine enjekte edilir.
    % Bu süreç iterasyon başlangıcında saklanan FlowerField vektörlerini (i.e., Flowers vektörleri) farklı bir yolla daha evolve etmeyi sağladığından iABC'nin exploitation (local search) yeteneğine katkıda bulunmaktadır
    [nectar1, fitnectar1] = Update(Flowers0, fitFlowers0, OnlookerFlowers, lb, ub, ObjFnc, userdata);
    
    % Update #3: Onlooker-Bee aşamasının son adımında FlowerField bilgileri (i.e., ilgili Flowers vektörleri ve karşılık gelen kalite değerleri, fitFlowers) güncellenir.
    j=fitnectar1<fitFlowers;
    FlowerField(j,:)=nectar1(j,:);
    fitFlowers(j)=fitnectar1(j);
     

    %% Scout-Bee Aşaması: iABC'nin FlowerField matrisinde sayısal çeşitlilik yoksunluğu probleminden kaçınmasını sağlaya süreçtir. 
    % Anolojik olarak ScoutBee olarak adlandırılan arıların doğrudan ulaşmadıkları wildflowers kaynaklarının varlıklarının farkında olduklarını ifade eder. 
    % ScoutBee nektar elde etmek için uçacağı yönü wildflowers konumlarından yararlanarak belirler.  wildflowers konumları her 100 iterayon adımında bir değiştirilmektedir. Bu durum ScoutBee vektörlerinin
    % numrical diversity-problem ile karşılaşmaksızın doğrultu vektörleri üretebilmelerini sağlamaktadır. 
    
    if rand < 0.10        
        [M,scale] = generateMS(N, D,0.90);  % Scout-Bee aşamasında nektar kaynakları rasgele genliklerle değişen farklı seviyelerde kayrılırlar. 
        % Bu durum değişkenler arasında kompleks ilişkiler olan strongly-related problemlerin çözümünü kolaylaştırmaktadır.
        dx = wildflowers - FlowerField; % ScoutBee için doğrultu vektörü
        ScoutBee = FlowerField + M .* scale .* dx;  % Scout-Bee aşaması için morphogenesis süreci
        % Update #4: Scout-Bee aşamasının son adımında FlowerField bilgileri (i.e., ilgili Flowers vektörleri ve karşılık gelen kalite değerleri, fitFlowers) güncellenir.
        [FlowerField, fitFlowers] = Update(FlowerField, fitFlowers, ScoutBee, lb, ub, ObjFnc, userdata);
    end
    
    % Update #5: Küresel çözümlerin güncellenemsi aşamasıdır
    % iABC algoritmasında, FlowerField'ın bulunulan anda sahip olduğu en iyi kalite seviyesine sahip Flowers vektörü en-iyi nektar kaynağıdır. 
    % iABC'de, en-iyi nektar kaynağı küresel en-iyi çözüm vektörüne (i.e., gBest) karşılık gelir. 
    % gBest'in  kalite seviyesi ise küresel en-iyi çözüm değerine (i.e., gBestVal) karşılık gelir.
    [fBest, idx] = min(fitFlowers);
    if fBest < gBestVal, gBestVal = fBest; gBest = FlowerField(idx,:); end
    
    % Küresel çözümlerin ekrana ve Matlab Workspace ortamına raporlanması
    out.gbest = gBest;   % Küresel en-iyi Flowers (nektar kaynağı)
    out.gval = gBestVal; % Küresel en-iyi Flowers vektörü kalitesi
           

    fprintf('Iter=%d -->  fMin=%5.16f\n', iteration, gBestVal);

    if mod(iteration,100)==0, wildflowers = lb + (ub -lb) .* rand(N, D); end % saves diversity of evolutionary-direction
end % iterasyon sonu

end % fonksiyon sonu



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