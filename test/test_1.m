%boostrap_test_1
%test to check if the estimated SE from the bootstrap is a reasonable
%estimation of the error

%define the distribution generator with a mock function
gen_data=@(x) normrnd(0,1,[1e2,1]);
data=gen_data(1);
%define the operation
anal_opp=@(x) mean(x);
analytic_dist_ste=1/sqrt(numel(data));
analytic_samp_se=std(data,1)/sqrt(numel(data));

%% Check that analytic SE is correct
%use a simple mutliple sampling of the distribution to give a Monte Carlo
%answer
dist_sample_res=[];
iimax=1e3;
for ii=1:iimax
data=gen_data(1);
dist_sample_res(ii)=anal_opp(data);
end
mc_samp_se=std(dist_sample_res);

LogicalStr = {'FAIL','pass'};
test_tol=0.1;
is_analysic_se_right= (abs(diff([mc_samp_se,analytic_dist_ste]))/...
    mean([mc_samp_se,analytic_samp_se]))<test_tol;
fprintf('analytic and MC standard error =±%.0f%% :%s \n',test_tol*1e2,LogicalStr{is_analysic_se_right+1})

%%


fignum=10;
boot=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',false,...
    'samp_frac_lims',[0.005,0.9],...
    'num_samp_frac',1e2,...
    'num_samp_rep',1e1,...
    'plot_fig_num',fignum,...
    'true_dist_se',analytic_dist_ste,...
    'true_samp_se',analytic_samp_se,...
    'mean_se_for_se_se',0);

test_tol=0.1;
is_boot_right= (abs(diff([boot.se_opp,analytic_dist_ste]))/...
    mean([boot.se_opp,analytic_dist_ste]))<test_tol;
fprintf('boot and analytic standard error =±%.0f%% :%s \n',test_tol*1e2,LogicalStr{is_analysic_se_right+1})


%%
%see if bootstrap SE is within a few of the estimated SE in the SE from the
%analytic value for the SE
test_tol=2;
is_boot_within_se= abs(diff([boot.se_opp,analytic_dist_ste]))<boot.se_se_opp*test_tol;
%? does not always render right
fprintf('boot and analytic standard error within %.0f sigma :%s \n',test_tol,LogicalStr{is_analysic_se_right+1})

%%

%histogram the residuals to see if the moemnt based error is roughly correct
sfigure(2);
subplot(2,1,1)
hist(abs((boot.se_opp-boot.opp_frac_est_se(:,2))./boot.opp_frac_est_se(:,3)),1e2)
subplot(2,1,2)
hist(abs((boot.se_opp-boot.opp_frac_est_se(:,2))./boot.std_se_opp_unweighted),1e2)


%compare the error in the SE found using the moments to that computed using the spread in data 
(boot.se_se_opp)/boot.se_se_opp_unweighted
%find the fraction error in the error estimation
(analytic_dist_ste-boot.se_opp)/analytic_dist_ste
%find the number of estimated SD the real value is away
(analytic_dist_ste-boot.se_opp)/boot.se_se_opp

%%

