%test_gauss_fit_bootstrap_error

paren_fun = @(x, varargin) x(varargin{:});
curly_fun = @(x, varargin) x{varargin{:}};

data=[randn([1e3,1]);rand([1e3,1])];

%test_gauss_fit_bootstrap_error_fit_wrapper(data)
%
anal_opp=@(x) paren_fun(test_gauss_fit_bootstrap_error_fit_wrapper(x),2);
%
real_samp_se=std(data,1)/sqrt(numel(data));
real_dist_ste=sqrt((1/12))/sqrt(numel(data));
%unifrm distributions give an oversetimated error in the SE

boot=bootstrap_se(anal_opp,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.1,1],...
    'samp_frac_log',1,...
    'num_samp_frac',1e2,...
    'num_samp_rep',1e2,...
    'verbose',10)

%find the fraction error in the error estimation
(real_dist_ste-boot.results.se_fun_whole)/real_dist_ste
%find the number of estimated SD the real value is away
(real_dist_ste-boot.results.se_fun_whole)/boot.results.se_se_fun_whole
