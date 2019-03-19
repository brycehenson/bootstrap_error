%ok trying to figure out how to estimate A(x) (opperated elementwise) from the moments of the 
%result of operation A(S) operated on many subsets


data=rand([1e4,1]);
anal_opp=@(x) mean(x);
%the unorrected standard deviation is a good example
%biased sample variance of the subset
var_biased=moment(data,2);
data_size=size(data,1);
unbias_samp_var=(data_size./(data_size-1)).* var_biased;
unbias_samp_std=sqrt(unbias_samp_var)

%%
%now we try to replicate this but only sampling larger subsets of the data
samples=1000;
sample_frac=0.1;
sample_number=round(data_size*sample_frac);
do_replace=true;
data_smpl=[];
for ii=1:samples
opp_smpl(ii)=anal_opp(randsample(data,sample_number,do_replace));
end

%%
var_biased_smpl=moment(opp_smpl,2);
%correct for the bias
unbias_samp_var_smpl=(samples./(samples-1)).* var_biased_smpl;
%correct for the average over everything in the subset
unbias_samp_var_smpl=unbias_samp_var_smpl*sample_number;
unbias_samp_std_smpl=sqrt(unbias_samp_var_smpl)


%% ok so does the same approach work for higher order moments?? 
% wa can correct for the finite sample size with http://mathworld.wolfram.com/SampleCentralMoment.html

data=(rand([1e6,1])-0.5)*1;
%data=(normrnd(0,1,[1e6,1])-0)*10;
anal_opp=@(x) mean(x);
%the unorrected standard deviation is a good example
%biased sample variance of the subset
mean_data=mean(data);
fourth_biased=moment(data-mean_data,4);
data_size=size(data,1);
fourth_unbiased=fourth_biased*(data_size)^2/((data_size-1)*(data_size-2));
fourth_unbiased

%%
%now we try to replicate this but only sampling larger subsets of the data
samples=5000;
sample_frac=0.1;
sample_number=round(data_size*sample_frac);
do_replace=true;
opp_smpl=[];
for ii=1:samples
opp_smpl(ii)=anal_opp(randsample(data,sample_number,do_replace));
end
mean_opp=mean(opp_smpl);
hist(opp_smpl,1e2)
fourth_biased_sub=moment(opp_smpl-mean_opp,4);
%correct for the bias
fourth_unbiased_sub=fourth_biased_sub*(samples^2)/((samples-1)*(samples-2));
%correct for the average over everything in the subset
fourth_unbiased_est=fourth_unbiased_sub*samples
%%


