%boostrap_test_1
% testing the standard error values produced from the bootstraping 
%test to check if the estimated SE from the bootstrap is a reasonable
%estimation of the error



%% NORMAL DISTRIBUTION

%define the distribution generator with a mock function
gen_data=@(x) normrnd(0,1,[1e2,1]); %the input does nothing
data=gen_data(1);
%define the analysis operation as the mean
anal_opp=@(x) mean(x);
%for this distribution the error in the meaan is simple
analytic_dist_ste=1/sqrt(numel(data));
%because we are only calculating the mean we can get what the statisitcal answer is for this
%dataset
analytic_samp_se=std(data,1)/sqrt(numel(data));

%% Check that analytic SE is correct
% just to doubly make sure lets check that the analysic prediction for the SE is correct
% to do so we draw datasets (Monte Carlo) from the distirbution and calculate the mean
% then the error is the std dev of these means
dist_sample_res=[];
iimax=1e3;
for ii=1:iimax
data=gen_data(1);
dist_sample_res(ii)=anal_opp(data);
end
mc_samp_se=std(dist_sample_res);

logical_str = {'FAIL','pass'};
test_tol=0.1;
is_analysic_se_right= (abs(diff([mc_samp_se,analytic_dist_ste]))/...
    mean([mc_samp_se,analytic_samp_se]))<test_tol;
fprintf('INFO: analytic SE = %.3f\n',analytic_dist_ste)
fprintf('INFO: MC SE       = %.3f\n',mc_samp_se)
fprintf('TEST: analytic = MC standard error   ±%.0f%% :%s \n',test_tol*1e2,logical_str{is_analysic_se_right+1})

%%



boot=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',false,...
    'samp_frac_lims',[0.005,0.9],...
    'num_samp_frac',1e2,...
    'num_samp_rep',1e2,...
    'true_dist_se',analytic_dist_ste,...
    'true_samp_se',analytic_samp_se,...
    'mean_se_for_se_se',0);

test_tol=0.1;
is_boot_right= frac_diff(boot.results.se_fun_whole,analytic_dist_ste)<test_tol;
fprintf('INFO: boot SE     = %.3f\n',boot.results.se_fun_whole)
fprintf('TEST: boot = analytic standard error ±%.0f%% :%s \n',test_tol*1e2,logical_str{is_boot_right+1})

%
%see if bootstrap SE is within a few of the estimated SE in the SE from the
%analytic value for the SE
test_tol=3;

fprintf('INFO: boot (sample frac) SE in SE = %.1e\n', boot.results.se_se_fun_whole_unweighted)
fprintf('INFO: boot (estimated) SE in SE   = %.1e\n', boot.results.se_se_fun_whole_weighted_norm)
fprintf('INFO: boot- analytic SE         = %.1e\n', diff([boot.results.se_fun_whole,analytic_dist_ste]))
fprintf('INFO: (boot- analytic SE)/(sample frac) std SE    = %.1g\n', diff([boot.results.se_fun_whole,analytic_dist_ste])/boot.results.se_se_fun_whole_unweighted)
is_boot_within_se= abs(diff([boot.results.se_fun_whole,analytic_dist_ste]))<boot.results.se_se_fun_whole*test_tol;
fprintf('TEST: boot = analytic standard error w/in %.0f (est)std :%s \n',test_tol,logical_str{is_boot_within_se+1})
fprintf('INFO: (boot- analytic SE)/estimated) std SE       = %.1g\n', diff([boot.results.se_fun_whole,analytic_dist_ste])/boot.results.se_se_fun_whole_weighted_norm)
is_boot_within_se= abs(diff([boot.results.se_fun_whole,analytic_dist_ste]))<boot.results.se_se_fun_whole_weighted_norm*test_tol;
fprintf('TEST: boot = analytic standard error w/in %.0f (est)std :%s \n',test_tol,logical_str{is_boot_within_se+1})

%%
% for each sample faction find the number of se that the predicted se in the total was different from the mean across
% sample fractions, normalized by the estimated se in the se from that sample fraction
%histogram the residuals to see if the moemnt based error is roughly correct
%
stfig('error histogram');
subplot(2,1,1)
hist((boot.sampling.projected_whole_se_norm_unbias-mean(boot.sampling.projected_whole_se_norm_unbias))./(boot.sampling.std_projected_whole_se_norm),1e2)
xlabel('number std (projected se for sample fraction - mean se)')
ylabel('counts')
title('normal distribution correction')
subplot(2,1,2)
hist((boot.sampling.std_projected_whole_se_norm-mean(boot.sampling.std_projected_whole_se_norm))./(boot.sampling.std_projected_whole_se_norm),1e2)
hist(abs((boot.se_opp-boot.opp_frac_est_se(:,2))./boot.std_se_opp_unweighted),1e2)
xlabel('number std (measured from sample frac)')
ylabel('counts')

%compare the error in the SE found using the moments to that computed using the spread in data 
(boot.se_se_opp)/boot.se_se_opp_unweighted
%find the fraction error in the error estimation
(analytic_dist_ste-boot.se_opp)/analytic_dist_ste
%find the number of estimated SD the real value is away
(analytic_dist_ste-boot.se_opp)/boot.se_se_opp

%%

