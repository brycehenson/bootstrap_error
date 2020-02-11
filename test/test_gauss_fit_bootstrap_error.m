%test_gauss_fit_bootstrap_error

data=randn([1e3,1]);



anal_opp=@(x) test_gauss_fit_bootstrap_error_fit_wrapper(x)

real_samp_se=std(data,1)/sqrt(numel(data));
real_dist_ste=sqrt((1/12))/sqrt(numel(data));
%unifrm distributions give an oversetimated error in the SE

boot=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.005,0.1],...
    'num_samp_frac',1e2,...
    'num_samp_rep',1e2,...
    'verbose',10)

%find the fraction error in the error estimation
(real_dist_ste-boot.results.se_fun_whole)/real_dist_ste
%find the number of estimated SD the real value is away
(real_dist_ste-boot.results.se_fun_whole)/boot.results.se_se_fun_whole
