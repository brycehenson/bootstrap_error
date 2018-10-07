%stat_bootstrap_tests
%develop statistical bootstraping to deterine the standard error of complicated functions
%and check it converges to the correct values for known distributions

%there are two ways to bootstrap
% 1. with replacements
%   * estimate the population sd by std(anal(data_sample_with_rep))*sqrt(n_sample)
% 2. without replacements
%  * estimate the population sd by
%  std(anal(data_sample_no_rep))*sqrt(n_sample)/finite_pop_correction
%calculates finite population correction
%https://www.jstor.org/stable/2340569?origin=crossref&seq=2#metadata_info_tab_contents

% Known BUGS/ Possible Improvements
%   -make a function that integerates all this into something that can be easily wrapped arround a
%   dataset
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-10-04

%repeatability
%rng(round(pi*exp(1)*1e3));

% %for simple things tell us what the actual stadnard error in the mean is

%normal data operation=mean
data= normrnd(0,1,[1e3,1]);
anal_opp=@(x) mean(x);
real_ste=1/sqrt(numel(data));

%unit interval uniform
%data=rand([1e2,1]);
%anal_opp=@(x) mean(x);
%real_ste=sqrt((1/12))/sqrt(numel(data));


sample_ste=std(data,1)/sqrt(numel(data));
sample_frac_vec=linspace(1e-2,1,60);
repeat_samp_prefactor=1e2;
unc_frac=NaN(numel(sample_frac_vec),2);
unc_pop=NaN(size(sample_frac_vec));
n_total=numel(data);

fprintf('%02u',0)
for ii=1:numel(sample_frac_vec)
    sample_frac=sample_frac_vec(ii);
    n_sample=floor(sample_frac*n_total);
    if n_sample>3
        %the scaling of std(std(x)) is 1/n so lets take more data at smaller sample_frac
        repeat_samp=round(repeat_samp_prefactor*1/sample_frac);
        anal_sample_no_rep=NaN(repeat_samp,1);
        anal_sample_with_rep=anal_sample_no_rep;
        for jj=1:repeat_samp
            anal_sample_with_rep(jj)=anal_opp(randsample(data,n_sample,true));
            anal_sample_no_rep(jj)=anal_opp(randsample(data,n_sample));
        end
        %estimape the pop std for with replacements
        unc_frac(ii,1)=std(anal_sample_with_rep)*sqrt(n_sample);
        finite_pop_corr=sqrt((n_total-n_sample)/(n_total-1));
        %this does not work well at subset fractions aproaching 1
        if finite_pop_corr>1e-4 
            unc_frac(ii,2)=std(anal_sample_no_rep)*sqrt(n_sample)/finite_pop_corr;
        end
    end
    fprintf('\b\b%02u',ii);
end
fprintf('..Done\n')
est_anal_unc=[];
mean_est_anal_unc=[];
est_anal_unc(:,1)=unc_frac(:,1)/sqrt(numel(data));
mean_est_anal_unc(1)=nanmean(est_anal_unc(:,1));
est_anal_unc(:,2)=unc_frac(:,2)/sqrt(numel(data));
mean_est_anal_unc(2)=nanmean(est_anal_unc(:,2));

%%



figure(1);
clf

plot(sample_frac_vec,est_anal_unc(:,1),'r');
hold on
plot(sample_frac_vec,est_anal_unc(:,2),'b');
xl=xlim(gca);
line(xl,[1,1]*nanmean(est_anal_unc(:,1)),'Color','r','LineWidth',2)
line(xl,[1,1]*nanmean(est_anal_unc(:,2)),'Color','b','LineWidth',2)
line(xl,[1,1]*real_ste,'Color','k','LineWidth',2)
line(xl,[1,1]*sample_ste,'Color','m','LineWidth',2)

hold off
legend('with replacement','without replacement','rep avg','no rep avg','dist SE','data sample SE')
xlabel('frac data')
ylabel('std subset * sqrt(n)')

%% Making a function out of it

fignum=10;
data=normrnd(0,1,[1e4,1]);
anal_opp=@(x) mean(x);
real_dist_ste=1/sqrt(numel(data));
real_samp_se=std(data,1)/sqrt(numel(data));

boot=boostrap_se(anal_opp,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.8],...
    'num_samp_frac',1e3,...
    'num_samp_rep',1e1,...
    'plot_fig_num',fignum,...
    'true_dist_se',real_dist_ste,...
    'true_samp_se',real_samp_se)


(real_dist_ste-boot.se_opp)/real_dist_ste
(real_dist_ste-boot.se_opp)/boot.se_se_opp
