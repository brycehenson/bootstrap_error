
data=rand([1e3,1]);
anal_opp=@(x) mean(x);
real_samp_se=std(data,1)/sqrt(numel(data));
real_dist_ste=sqrt((1/12))/sqrt(numel(data));
%unifrm distributions give an oversetimated error in the SE


[boot,boot_detailed]=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.1],...
    'num_samp_frac',1e2,...
    'num_samp_rep',1e2,...
    'true_dist_se',real_dist_ste,...
    'true_samp_se',real_samp_se,...
    'verbose',10,...
    'save_multi_out',1)

%find the fraction error in the error estimation
(real_dist_ste-boot.results.se_fun_whole)/real_dist_ste
%find the number of estimated SD the real value is away
(real_dist_ste-boot.results.se_fun_whole)/boot.results.se_se_fun_whole


%%
data=rand([1e3,1]);
anal_opp=@(x) mean(x);
real_samp_se=std(data,1)/sqrt(numel(data));
real_dist_ste=sqrt((1/12))/sqrt(numel(data));
%unifrm distributions give an oversetimated error in the SE


boot=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.1],...
    'num_samp_frac',1e2,...
    'num_samp_rep',10,...
    'verbose',10)

%find the fraction error in the error estimation
(real_dist_ste-boot.results.se_fun_whole)/real_dist_ste
%find the number of estimated SD the real value is away
(real_dist_ste-boot.results.se_fun_whole)/boot.results.se_se_fun_whole

%% try with normaly distributed data
data=randn([1e3,1]);
anal_opp=@(x) mean(x);
real_samp_se=std(data,1)/sqrt(numel(data));
real_dist_ste=1/sqrt(numel(data));
%unifrm distributions give an oversetimated error in the SE


boot=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.1],...
    'num_samp_frac',30,...
    'num_samp_rep',30,...
    'true_dist_se',real_dist_ste,...
    'true_samp_se',real_samp_se,...
    'verbose',10)

%find the fraction error in the error estimation
(real_samp_se-boot.results.se_fun_whole)/real_dist_ste
%find the number of estimated SD the real value is away
(real_samp_se-boot.results.se_fun_whole)/boot.results.se_se_fun_whole

%% Test that it can bootstrap with a single sample size
data=randn([1e3,1]);
anal_opp=@(x) mean(x);
real_samp_se=std(data,1)/sqrt(numel(data));
real_dist_ste=1/sqrt(numel(data));
%unifrm distributions give an oversetimated error in the SE


boot=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.1],...
    'num_samp_frac',1,...
    'num_samp_rep',30,...
    'true_dist_se',real_dist_ste,...
    'true_samp_se',real_samp_se,...
    'verbose',10)

%find the fraction error in the error estimation
(real_samp_se-boot.results.se_fun_whole)/real_dist_ste
%find the number of estimated SD the real value is away
(real_samp_se-boot.results.se_fun_whole)/boot.results.se_se_fun_whole

%% Repeat to see what the distibution in the output values are
est_se=[];
est_se_se=[];
sigma_err=[];
fprintf('%04u',0)
for ii=1:100
data=rand([1e3,1]);
anal_opp=@(x) mean(x);
real_samp_se=std(data,1)/sqrt(numel(data));
real_dist_ste=sqrt((1/12))/sqrt(numel(data));
%unifrm distributions give an oversetimated error in the SE

boot=bootstrap_se(anal_opp,data,...
    'plots',false,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.1],...
    'num_samp_frac',1e2,...
    'num_samp_rep',1e2,...
    'true_dist_se',real_dist_ste,...
    'true_samp_se',real_samp_se,...
    'verbose',0);

est_se(ii)=boot.results.se_fun_whole_weighted_arb;
est_se_se(ii)=boot.results.se_se_fun_whole_weighted_arb;
sigma_err(ii)=(real_dist_ste-boot.results.se_fun_whole_weighted_arb)/boot.results.se_se_fun_whole_weighted_arb;
fprintf('\b\b\b\b\b%04u',ii)
end
fprintf('\n')
%%
stfig('results distribution')
histogram(col_vec(sigma_error),round(numel(sigma_error)/5))
xlabel('number of standard deviations')
fprintf('std of est se values %f, mean se se %f \n',std(est_se),mean(est_se_se))
fprintf('sigma error from mean %f\n',std(sigma_err))
fprintf('val sd / mean err  %f\n',std(est_se)/mean(est_se_se))



%%

data=num2cell(rand([1e2,2]),2);
d_flat=cell2mat(data);
d_flat=d_flat(:);
real_samp_se=std(d_flat,1)/sqrt(numel(d_flat));
real_dist_ste=sqrt((1/12))/sqrt(numel(d_flat));
%unifrm distributions give an oversetimated error in the SE

fignum=10;
boot=bootstrap_se(@(x) mean(mean(cell2mat(x))),data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.9],...
    'num_samp_frac',1e1,...
    'num_samp_rep',1e2,...
    'true_dist_se',real_dist_ste,...
    'true_samp_se',real_samp_se,...
    'verbose',10)

%histogram the residuals to see if the moemnt based error is roughly correct
%compare the error in the SE found using the moments to that computed using the spread in data 
(boot.results.se_se_fun_whole_weighted_arb-boot.results.se_se_fun_whole_unweighted)/...
    boot.results.se_se_fun_whole_unweighted
%find the fraction error in the error estimation
(real_dist_ste- boot.results.se_fun_whole)/real_dist_ste
%find the number of estimated SD the real value is away
(real_dist_ste- boot.results.se_fun_whole)/boot.results.se_se_fun_whole_weighted_arb


%%


data=rand([1e3,1])-0.5;
gen_norm_measure=@(x,n) norm(x,n)*(numel(x)^(-1/n));
anal_opp=@(x) gen_norm_measure(x,6); %norm(x,2);
%unifrm distributions give an oversetimated error in the SE


boot=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.5,1.5],...
    'num_samp_frac',1e2,...
    'num_samp_rep',1e2,...
    'verbose',10)

%compare the error in the SE found using the moments to that computed using the spread in data 
(boot.results.se_fun_whole_unweighted-boot.results.se_fun_whole_weighted_arb)/boot.results.se_fun_whole_weighted_arb

%%

data=num2cell(rand([1e2,2]),2);
d_flat=cell2mat(data);
d_flat=d_flat(:,1);
real_samp_se=std(d_flat,1)/sqrt(numel(d_flat));
real_dist_ste=sqrt((1/12))/sqrt(numel(d_flat));
%unifrm distributions give an oversetimated error in the SE

fignum=10;
boot=bootstrap_se(@first_elm_avg,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.9],...
    'num_samp_frac',1e2,...
    'num_samp_rep',1e2,...
    'true_dist_se',real_dist_ste,...
    'true_samp_se',real_samp_se,...
    'save_multi_out',1,...
    'verbose',10)

%compare the error in the SE found using the moments to that computed using the spread in data 
(boot.results.se_fun_whole_unweighted-boot.results.se_fun_whole_weighted_arb)/boot.results.se_fun_whole_weighted_arb
(real_dist_ste- boot.results.se_fun_whole)/real_dist_ste
% number of sigma away from true value
(real_dist_ste- boot.results.se_fun_whole)/boot.results.se_se_fun_whole


function [out1,out2]=first_elm_avg(in) 
    d_all=cell2mat(in);
    out1=mean(d_all(:,1));
    out2=mean(d_all(:,2));
end

function out=both_elm_avg(in) 
    d_flat=cell2mat(in);
    d_flat=d_flat(:);
    out=mean(d_flat);
end
