clear all; clc; close all;

k = 50;
r = 50:10:170;

load('final_protein_result.mat');

norm_rel2 = median(normSamp_Result.relErr2);
norm_relF = median(normSamp_Result.relErrF);
norm_time = median(normSamp_Result.timetot);

dct_rel2 = median(DCT_Result.relErr2);
dct_relF = median(DCT_Result.relErrF);
dct_time = median(DCT_Result.timetot);

spemb_rel2 = median(Spemb_Result.relErr2);
spemb_relF = median(Spemb_Result.relErrF);
spemb_time = median(Spemb_Result.timetot);

freq_rel2 = median(FreqDir_Result.relErr2);
freq_relF = median(FreqDir_Result.relErrF);
freq_time = median(FreqDir_Result.timetot);

spfd5_rel2 = median(RandFreqDir5_Result.relErr2);
spfd5_relF = median(RandFreqDir5_Result.relErrF);
spfd5_time = median(RandFreqDir5_Result.timetot);

spfd10_rel2 = median(RandFreqDir10_Result.relErr2);
spfd10_relF = median(RandFreqDir10_Result.relErrF);
spfd10_time = median(RandFreqDir10_Result.timetot);

spfd50_rel2 = median(RandFreqDir50_Result.relErr2);
spfd50_relF = median(RandFreqDir50_Result.relErrF);
spfd50_time = median(RandFreqDir50_Result.timetot);



figure(1)   % relative error 2 norm
plot(r, norm_rel2, 'bo-', r,dct_rel2,'c*-', r,spemb_rel2,'rx-', r, freq_rel2, 'ks-', r,spfd5_rel2,'gd-',...
    r,spfd10_rel2,'gd--', r,spfd50_rel2,'gd:','LineWidth',1.5);
set(gca, 'fontsize', 18);
ylabel('Relative error (2 norm)', 'interpreter','LaTex', 'fontsize', 18);
xlabel('sketch size $\ell$','interpreter','LaTex', 'fontsize', 18);
l = legend('NormSamp', 'DCT', 'SpEmb', 'FD', 'SpFD5', 'SpFD10', 'SpFD50');
set(l, 'fontsize', 14)

figure(2)   % relative error F norm
plot(r, norm_relF, 'bo-', r,dct_relF,'c*-', r,spemb_relF,'rx-', r, freq_relF, 'ks-', r,spfd5_relF,'gd-',...
    r,spfd10_relF,'gd--', r,spfd50_relF,'gd:','LineWidth',1.5);
set(gca, 'fontsize', 18);
ylabel('Relative error (F norm)', 'interpreter','LaTex', 'fontsize', 18);
xlabel('sketch size $\ell$','interpreter','LaTex', 'fontsize', 18);
l = legend('NormSamp', 'DCT', 'SpEmb', 'FD', 'SpFD5', 'SpFD10', 'SpFD50');
set(l, 'fontsize', 14);

figure(3)   % total time
plot(r, norm_time, 'bo-', r,dct_time,'c*-', r,spemb_time,'rx-', r, freq_time, 'ks-', r,spfd5_time,'gd-',...
    r,spfd10_time,'gd--', r,spfd50_time,'gd:','LineWidth',1.5);
set(gca, 'fontsize', 18);
ylabel('Running Time (secs)', 'interpreter','LaTex', 'fontsize', 18);
xlabel('sketch size $\ell$','interpreter','LaTex', 'fontsize', 18);
l = legend('NormSamp', 'DCT', 'SpEmb', 'FD', 'SpFD5', 'SpFD10', 'SpFD50');
set(l, 'fontsize', 14);


