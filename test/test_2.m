%test_2
%test the estimates of the error in the esimated SE by nesting the
%boostraper

%too do 
%   - look at derivation of error agian
%   - handle no fractions case

%define the distribution generator with a mock function
gen_data=@(x) normrnd(0,1,[1e3,1]);
data=gen_data(1);
%define the operation
anal_opp=@(x) mean(x);


fignum=50;
[boot,nested_boot_record]=bootstrap_se(@wrapped_boot,data,...
    'plots',true,...
    'replace',true,...
    'samp_frac_lims',[0.05,0.9],...
    'num_samp_frac',10,...
    'num_samp_rep',1e2,...
    'plot_fig_num',fignum,...
    'save_multi_out',true,...
    'mean_se_for_se_se',0,...
    'opp_arguments',{anal_opp});


%now i try and see which metric was the best at identifying the se_se


%%
a=nested_boot_record.out_cell(4,2);
a{1}
%%
iimax=size(nested_boot_record.out_cell,1);
moment_est_se_se=nan(iimax,1);
for ii=1:iimax
    cstat=nested_boot_record.out_cell(ii,:);
    stat_tmp=arrayfun(@(x) cstat{x}.se_se_opp,1:numel(cstat));
    moment_est_se_se(ii,1)=mean(stat_tmp);
    moment_est_se_se(ii,2)=std(stat_tmp);
end
figure(51) 
pred_ratio=moment_est_se_se./bootstrap_est_se_se;
pred_ratio=pred_ratio.*nested_boot_record.sample_num_vec/numel(data);
plot(nested_boot_record.sample_frac_vec,pred_ratio(:,1),'k')
hold on
plot(nested_boot_record.sample_frac_vec,pred_ratio(:,1)+pred_ratio(:,2),'r')
plot(nested_boot_record.sample_frac_vec,pred_ratio(:,1)-pred_ratio(:,2),'r')
hold off
ylabel('predicted se in the se vs the measured value')
xlabel('fraction of data')
bootstrap_est_se_se=boot.se_opp;
moment_est_se_se/bootstrap_est_se_se

%%
function [se_val,boot]=wrapped_boot(data,anal_opp)

fignum=9;
analytic_dist_ste=1/sqrt(numel(data));
boot=bootstrap_se(anal_opp,data,...
    'plots',false,...
    'replace',true,...
    'samp_frac_lims',0.4,...%[0.005,0.9]
    'num_samp_frac',1e1,...
    'num_samp_rep',1e1,...
    'plot_fig_num',fignum,...
    'true_dist_se',analytic_dist_ste,...
    'save_multi_out',true,...
    'verbose',false,...
    'mean_se_for_se_se',0);

se_val=boot.se_opp;
end