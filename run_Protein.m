% A = SDU + N/e: nxd
% S: nxm, S_ij ~ N(0,1)
% D: mxm, D_ii = 1-(i-1)/m
% U: mxd, random and UtU = Id
% N: nxd, N_ij ~ N(0,1) 
% Set d = 1000; k=m;
% n = 10,000, ..., 100,000;
% m = 10, 20, 50;
% e = 5, 10, 15, 20
% Methods: 
% Norm Sampling
% DCT
% SpEmb
% FD (FreqDir)
% SpFD (RandFreqDir)

clear; clc;
repj = 5;
rk = 50;
load('protein.mat');
n = size(A,1);
filename = sprintf('%s.mat', 'final_protein_result');
save(filename, 'rk');
sketch = rk:10:170;
len_sk1 = length(sketch);

Spemb_Result.relErr2 = zeros(repj, len_sk1);
Spemb_Result.relErrF = zeros(repj, len_sk1);
Spemb_Result.norm2 = zeros(repj, len_sk1);
Spemb_Result.normF = zeros(repj, len_sk1);
Spemb_Result.timeV = zeros(repj, len_sk1);
Spemb_Result.timeform = zeros(repj, len_sk1);
Spemb_Result.timetot = zeros(repj, len_sk1);

normSamp_Result = Spemb_Result;
DCT_Result = Spemb_Result;

% alp = 1 -> FD, alp = n/l -> spemb
% n = 10000, l = 10:200, alp = 1: 1000, 1:50, k = 10, n/l*
% pick three alp = n/l/5, n/l/10, n/l/50
% This implies cut n into 5, 10, 50 parts and map them to l. 
alp = [5, 10, 50];
RandFreqDir5_Result = Spemb_Result;
RandFreqDir10_Result = Spemb_Result;
RandFreqDir50_Result = Spemb_Result;
% FreqDir no need to repeat its procedure as it is not random, but as time
% is not always accurate, we repeat it.
FreqDir_Result = Spemb_Result;
MSFD_Result = Spemb_Result;

% In this, we set l = r; this is for the sake of having the sketch of the
% same dimension, so use freq(A, tarNumRow(i)), this way, FD & SpFD will be
% even more accurate, but sacrifice the running time.
% If we set 2l = r, then use freq(A, 1/2*tarNumRow(i))
% This way, FD & SpFD will be faster but less accurate. But this way can be
% explained as same memory space. 

% Pick V(:,1:l) this is to ensure pick from range space with the same
% dimension, very bookish
% in real application, choose 2l which could be implemented in when 2l = r.

% Exact method

Exact.time = zeros(3, 1);
for i = 1:3
    fprintf('SVD starts...\n')
    tic
    [u, s, v] = svd(A, 'econ');
    Ak = u(:, 1:rk) * s(1:rk,1:rk) * v(:,1:rk)';
    Exact.time(i) = toc;
    Exact.norm2 = s(rk+1,rk+1);
    Exact.normF = sqrt(sum(diag(s(rk+1:end, rk+1:end)^2)));
    save(filename, 'Exact', '-append'); 
end

for j = 1:repj
    
    % note for freq direction method, they keep the top l columns of V,
    % even for the last step, V does have 2l columns
    for i = 1:len_sk1
        l = sketch(i);
        
        fprintf('%d %d normSamp starts ...\n', i, j)
        [A_rand, time_rand] = normSamp(A, l);
        [Ak_rand, time1, timet] = randlowrank(A, A_rand, rk);
        normSamp_Result.norm2(j, i) = norm(A-Ak_rand, 2);
        normSamp_Result.normF(j, i) = norm(A-Ak_rand,'fro');
        normSamp_Result.relErr2(j, i) = normSamp_Result.norm2(j, i)/Exact.norm2;
        normSamp_Result.relErrF(j, i) = normSamp_Result.normF(j, i)/Exact.normF;
        normSamp_Result.timeform(j, i) = time_rand;        
        normSamp_Result.timeV(j, i) = time_rand+time1;        
        normSamp_Result.timetot(j, i) = time_rand+timet;
        save(filename, 'normSamp_Result', '-append');
        

        fprintf('%d %d DCT starts ...\n', i, j)
        [A_rand, time_rand] = DisCT(A, l);
        [Ak_rand, time1, timet] = randlowrank(A, A_rand, rk);
        DCT_Result.norm2(j, i) = norm(A-Ak_rand,2);
        DCT_Result.normF(j, i) = norm(A-Ak_rand,'fro');
        DCT_Result.relErr2(j, i) = DCT_Result.norm2(j, i)/Exact.norm2;
        DCT_Result.relErrF(j, i) = DCT_Result.normF(j, i)/Exact.normF;
        DCT_Result.timeform(j, i) = time_rand;           
        DCT_Result.timeV(j, i) = time_rand+time1;        
        DCT_Result.timetot(j, i) = time_rand+timet;
        save(filename, 'DCT_Result', '-append');

        fprintf('SpEmb starts...\n')
        [A_rand, time_rand] = Sparse(A, l);
        [Ak_rand, time1, timet] = randlowrank(A, A_rand, rk);
        Spemb_Result.norm2(j, i) = norm(A-Ak_rand,2);
        Spemb_Result.normF(j, i) = norm(A-Ak_rand,'fro');
        Spemb_Result.relErr2(j, i) = Spemb_Result.norm2(j, i)/Exact.norm2;
        Spemb_Result.relErrF(j, i) = Spemb_Result.normF(j, i)/Exact.normF;
        Spemb_Result.timeform(j, i) = time_rand; 
        Spemb_Result.timeV(j, i) = time_rand+time1;        
        Spemb_Result.timetot(j, i) = time_rand+timet;
        save(filename, 'Spemb_Result', '-append');
       
        fprintf('FreqDir starts...\n')
        [~, V_fr, time_fr] = freqDir(A, l);
        [Ak_freqV, time_freqV] = lowrankV(A, V_fr(:,1:l), rk);
        FreqDir_Result.norm2(j, i) = norm(A-Ak_freqV, 2);
        FreqDir_Result.normF(j, i) = norm(A-Ak_freqV, 'fro');
        FreqDir_Result.relErr2(j, i) = FreqDir_Result.norm2(j, i)/Exact.norm2;
        FreqDir_Result.relErrF(j, i) = FreqDir_Result.normF(j, i)/Exact.normF;
        FreqDir_Result.timeV(j, i) = time_fr;
        FreqDir_Result.timetot(j, i) = time_fr + time_freqV;
        save(filename, 'FreqDir_Result', '-append');
        
        fprintf('Mina SFD starts...\n')
        [~, V_mfr, time_msfr] = MSFD(A, l);
        [Ak_MSFDV, time_mV] = lowrankV(A, V_mfr(:, 1:l), rk); 
        MSFD_Result.norm2(j, i) = norm(A-Ak_MSFDV, 2);
        MSFD_Result.normF(j, i) = norm(A-Ak_MSFDV, 'fro');
        MSFD_Result.relErr2(j, i) = MSFD_Result.norm2(j, i)/Exact.norm2;
        MSFD_Result.relErrF(j, i) = MSFD_Result.normF(j, i)/Exact.normF;
        MSFD_Result.timeV(j, i) = time_msfr;
        MSFD_Result.timetot(j, i) = time_msfr + time_mV;
        save(filename, 'MSFD_Result', '-append');

        fprintf('RandFreqDir 5 cats starts ... \n')
        [~, V_rfr, time_rfr] = randFreqDirP(A, l, ceil(n/alp(1)));
        [Ak_rfreq, time_rfreq] = lowrankV(A, V_rfr(:,1:l), rk);
        RandFreqDir5_Result.norm2(j, i) = norm(A-Ak_rfreq, 2);
        RandFreqDir5_Result.normF(j, i) = norm(A-Ak_rfreq, 'fro');
        RandFreqDir5_Result.relErr2(j, i) = RandFreqDir5_Result.norm2(j, i)/Exact.norm2;
        RandFreqDir5_Result.relErrF(j, i) = RandFreqDir5_Result.normF(j, i)/Exact.normF;
        RandFreqDir5_Result.timeV(j, i) = time_rfr;
        RandFreqDir5_Result.timetot(j, i) = time_rfr + time_rfreq;
        save(filename, 'RandFreqDir5_Result', '-append');

        fprintf('RandFreqDir 10 cats starts ... \n')
        [~, V_rfr, time_rfr] = randFreqDirP(A, l, ceil(n/alp(2)));
        [Ak_rfreq, time_rfreq] = lowrankV(A, V_rfr(:,1:l), rk);
        RandFreqDir10_Result.norm2(j, i) = norm(A-Ak_rfreq, 2);
        RandFreqDir10_Result.normF(j, i) = norm(A-Ak_rfreq, 'fro');
        RandFreqDir10_Result.relErr2(j, i) = RandFreqDir10_Result.norm2(j, i)/Exact.norm2;
        RandFreqDir10_Result.relErrF(j, i) = RandFreqDir10_Result.normF(j, i)/Exact.normF;
        RandFreqDir10_Result.timeV(j, i) = time_rfr;
        RandFreqDir10_Result.timetot(j, i) = time_rfr + time_rfreq;
        save(filename, 'RandFreqDir10_Result', '-append');

        fprintf('RandFreqDir 50 cats starts ... \n')
        [~, V_rfr, time_rfr] = randFreqDirP(A, l, ceil(n/alp(3)));
        [Ak_rfreq, time_rfreq] = lowrankV(A, V_rfr(:,1:l), rk);
        RandFreqDir50_Result.norm2(j, i) = norm(A-Ak_rfreq, 2);
        RandFreqDir50_Result.normF(j, i) = norm(A-Ak_rfreq, 'fro');
        RandFreqDir50_Result.relErr2(j, i) = RandFreqDir50_Result.norm2(j, i)/Exact.norm2;
        RandFreqDir50_Result.relErrF(j, i) = RandFreqDir50_Result.normF(j, i)/Exact.normF;
        RandFreqDir50_Result.timeV(j, i) = time_rfr;
        RandFreqDir50_Result.timetot(j, i) = time_rfr + time_rfreq;
        save(filename, 'RandFreqDir50_Result', '-append');
    end
    


end