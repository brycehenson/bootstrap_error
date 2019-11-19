function [out,detailed_data_out]=bootstrap_se(est_fun,data,varargin)
%funciton calculates the standard error of performing est_fun(data)
%data is a vector of cells or scalars

%est_fun (estimator opperation) is a function that takes a vector or a cell matrix and retuns a scalar as the first output
%this code can store an arbitraty number of function outputs in a cell array and return them for use
% Known BUGS/ Possible Improvements
%   -handle arbitrary number of scalars (as a vector or matrix) for simulanious multi paramerter boostraping
%   - [ ] simplify codeflow with subfunctions
%   - [x] add in functionality to look at the mean with the sample fraction
%   - [x] add in biasing corrections
%     - [x] corrections for standard deviation assuming normality
%   - [x] error in the error
%     - [x] unbiased std of the estimated std of the distribution assuming normality
%     - [x] approx std of the estimated std for arb distribution
%     - [x] c4 unbiasing function
%   - [ ] multipe output from estimation function and simulanious boostrap
% Author: Bryce Henson
% email: Bryce.Henson[the a swirly]live.com you must have
% '[matlab][bootstrap_error]' in the subject line for me to respond
% Last revision:2019-06-02

min_sample_num=5; %minimum data size for resampling

p = inputParser;
is_lims=@(x) (isequal(size(x),[1,2]) && (isnumeric(x) && x(2)<10 && x(1)>0))...
        || (isscalar(x) && x<=10);
is_c_logical=@(in) isequal(in,true) || isequal(in,false); %can x be cast as a logical
addOptional(p,'replace',true,is_c_logical);
addOptional(p,'num_samp_frac',10,@(x) isnumeric(x) & x>=1);
addOptional(p,'num_samp_rep',10,@(x) isnumeric(x) & x>=1);
addOptional(p,'samp_frac_lims',0.1,is_lims);
addOptional(p,'plots',false,is_c_logical);
%optional arguments for diagnostic/test plots
addOptional(p,'true_dist_se',nan,@(x) isnumeric(x) & x>0);
addOptional(p,'true_samp_se',nan,@(x) isnumeric(x) & x>0);
addOptional(p,'mean_se_for_se_se',false,is_c_logical);
addOptional(p,'save_multi_out',false,is_c_logical);
addOptional(p,'save_input_data',false,is_c_logical);
addOptional(p,'verbose',1,@(x) isscalar(x) & x>=0);
addOptional(p,'opp_arguments',{},@(x) iscell(x) & size(x,1)==1);
addOptional(p,'plot_fig_name','',@ischar);
addOptional(p,'do_mean_fit',1,is_c_logical)
addOptional(p,'use_weighting',true,is_c_logical)
addOptional(p,'norm_weighting',false,is_c_logical)
parse(p,varargin{:});

do_mean_fit=p.Results.do_mean_fit;
do_plots=coerce_logical(p.Results.plots);
do_replace=coerce_logical(p.Results.replace);
use_mean_se_for_se_se=coerce_logical(p.Results.mean_se_for_se_se);
use_weighting=coerce_logical(p.Results.use_weighting);
norm_weighting=coerce_logical(p.Results.norm_weighting);

if use_mean_se_for_se_se
    warning('mean_se_for_se_se is depricated')
end

save_multi_out=coerce_logical(p.Results.save_multi_out);
save_input_data=coerce_logical(p.Results.save_input_data);
verbose=p.Results.verbose;
opp_arguments=p.Results.opp_arguments;
repeat_samp_prefactor=p.Results.num_samp_rep;
%input taken care of

%if the number of fractions to sample is one take the mean of the limits
if p.Results.num_samp_frac==1
     sample_frac_vec=mean(p.Results.samp_frac_lims);
    %sample over multiple fractions
elseif size(p.Results.samp_frac_lims,2)==2 
    sample_frac_vec=linspace(p.Results.samp_frac_lims(1),...
    p.Results.samp_frac_lims(2),p.Results.num_samp_frac)';
%if the lims are a scalar then just do one fraction

end
%overwrite plots if only samplign at one fraction of the data
%if size(sample_frac_vec,1)==1
%   do_plots=false;
%end

n_total=numel(data);
sample_num_vec=floor(sample_frac_vec*n_total);
sample_frac_vec=sample_num_vec/n_total;
%cull anything below min_num
mask=sample_num_vec>=min_sample_num;
sample_num_vec=sample_num_vec(mask);
sample_frac_vec=sample_frac_vec(mask);
iimax=numel(sample_frac_vec);

% the number of repeated sampling for a given sample size
% i have built the ability for this to vary but decided against it to prevent biasing
repeat_samp=repeat_samp_prefactor+0*sample_frac_vec; 
            
%catch when there are no valid sample sizes
if iimax==0
    error('no sample fracions')
end


%find the output size of the passed estimator function as in function [a,b,c]=function(inputs)
%should build in optional argument to specify this
est_fun_output_size=nargout(est_fun);
%this will only take the first output of an an anonymous function
if est_fun_output_size==-1
    est_fun_output_size=1;
end

out_cell=cell(iimax,repeat_samp_prefactor);
in_cell=out_cell;
out_cell_tmp=cell(1,est_fun_output_size);
if est_fun_output_size>1 && save_multi_out
    est_fun_multi_out=true;
else
    est_fun_multi_out=false;
end

%% find the size of the scalar output by calling the estimation function once
data_smpl=randsample(data,min_sample_num,do_replace);
if est_fun_multi_out
    [out_cell_tmp{:}]=est_fun(data_smpl,opp_arguments{:});
    out_val_tmp=out_cell_tmp{1};
else
    out_val_tmp=est_fun(data_smpl,opp_arguments{:});
end 
if ~isvector(out_val_tmp)
    erro('first output of function is not a scalar or vector')
end
out_val_tmp=col_vec(out_val_tmp);
output_val_size=numel(out_val_tmp);
%%
%prealocate the moments of the distribution
mean_sub=NaN(numel(sample_frac_vec),output_val_size);
moments_sub=NaN(numel(sample_frac_vec),3,output_val_size);


if verbose>0
    fprintf('Bootstrapping with different sample fractions %04u:%04u',0)
end
for ii=1:iimax
    n_sample=sample_num_vec(ii);
    %std means nothing for n<3
    %the finte sample correaction for the no replacements method breaks when n_sample=ntot
    if n_sample>3 && (n_sample<n_total || do_replace)
        est_fun_res_sub=NaN(repeat_samp(ii),output_val_size); %the results of the estimation function
        if do_replace
            finite_pop_corr=1;
        else
            finite_pop_corr=(n_total-n_sample)/(n_total-1);
        end 
        for jj=1:repeat_samp(ii)
            %calculate the analysis operation on the subset of dat
            data_smpl=randsample(data,n_sample,do_replace);
            %then we assign the output of the est_fun on the data_smp to est_fun_res_sub
            %matalb cant do [out_cell_tmp{:}]=scalar so a case statement
            %is needed
            if est_fun_multi_out
                [out_cell_tmp{:}]=est_fun(data_smpl,opp_arguments{:});
                out_val_tmp=col_vec(out_cell_tmp{1});
                if numel(out_val_tmp)~=output_val_size
                    error('output size wrong')
                end
                est_fun_res_sub(jj,:)=out_val_tmp;
                out_cell{ii,jj}=out_cell_tmp{2:end};
            else
                out_val_tmp=est_fun(data_smpl,opp_arguments{:});
                out_val_tmp=col_vec(out_val_tmp);
                if numel(out_val_tmp)~=output_val_size
                    error('output size wrong')
                end
                est_fun_res_sub(jj,:)=out_val_tmp;
            end 
            if save_input_data
                in_cell{ii,jj}=data_smpl;
            end
        end
        % calculate statistics on the results with a given sample size
        mean_sub(ii,:)=mean(est_fun_res_sub,1);
        % biased sample variance of the subset
        moments_sub(ii,1,:)=moment(est_fun_res_sub,2,1);
        moments_sub(ii,2,:)=moment(est_fun_res_sub,3,1);
        moments_sub(ii,3,:)=moment(est_fun_res_sub,4,1);

        %use finte population correction to estimate the population std using sampling without
        %replacements, if do_replace finite_pop_corr=1
        moments_sub(ii,1,:)=moments_sub(ii,1,:)/finite_pop_corr;
    end
    if verbose>0, fprintf('\b\b\b\b%04u',ii), end
end
 fprintf('\n')
 
out=[];
% apply the correction for the central moment
unbias_factor=(repeat_samp./(repeat_samp-1));
unbias_moments_sub= moments_sub.*repmat(unbias_factor,1,size(moments_sub,2));
%unbiased sample variance for the results of each bootstrap of a given size:
unbias_samp_var=unbias_moments_sub(:,1,:);
unbias_samp_var=permute(unbias_samp_var,[1,3,2]); % permute the dimensions to get the rid of the singleton dim
% biased sample standard deviation
std_est_subsamp=sqrt(unbias_samp_var); 

% baised sample standard error for each bootstrap size size()=[numel(sample_frac_vec), output_val_size]
ste_est_subsamp=std_est_subsamp./repmat(sqrt(repeat_samp),[1,output_val_size]);

%now calulate the standard error in the anal.
%operation on the whole datset assuming mean like scaling
% if this is flat with sample size then the estimator is mean-like (which is good)
mean_like_scaling_factor=sqrt(sample_num_vec)./sqrt(n_total);
mean_like_scaling_factor=repmat(mean_like_scaling_factor,[1,output_val_size]);
est_se_opp=std_est_subsamp.*mean_like_scaling_factor;

%% unbiasing normaly distributed data
% correct for the bias of the std https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
% assumes that the output is normaly distributed
unbias_normal_std_est_subsamp=std_est_subsamp./repmat(normal_correction_c4(repeat_samp),[1,output_val_size]);
ste_est_subsamp_normunbias=unbias_normal_std_est_subsamp./repmat(sqrt(repeat_samp),[1,output_val_size]);
est_se_opp_norm_unbias=unbias_normal_std_est_subsamp.*mean_like_scaling_factor;

%now we wish to give our best guess on what the standard error on the whole dataset will be
% as a first guess we could take the mean across sampling fraction
est_se_opp_mean_unweighted=mean(est_se_opp,1,'omitnan');
% and the asociated error
est_std_opp_se_unweighted= std(est_se_opp,[],1);
est_se_opp_se_unweighted=est_std_opp_se_unweighted./sqrt(size(est_se_opp,1));

%% finding the standard devitation for the estimated standard deviation in the whole dataset
% weighted SE(fun(data)),SE(SE(fun(data)))
% we could weight each estimated value, to do that we need some way of estimating what the expected standard deviation
% is for estimating the standard deviaiton in each subsample
% TODO come up with reasonable estimate of the se in the sample std
%https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance


% if we assume normality then things are pretty easy
%var_samp_var_norm=2*(std_est_subsamp.^4)/(repeat_samp-1);
% now we will do a crude aproximation and propagate this as var[s^2] to std(sigma)
%using a=b^(1/2) sigma(a)=(1/2)*sigma(b)*b^(-1/2)  sigma(a)=(1/2)*sigma(b)/a
%se_samp_std_norm_approx=(1/2)*sqrt(var_samp_var_norm)./std_est_subsamp;
%std_se_opp_norm_approx=se_samp_std_norm_approx*mean_like_scaling_factor;

% to do this properly we could do a change of variables on the chi-squared distribution
%https://en.wikipedia.org/wiki/Probability_density_function#Dependent_variables_and_change_of_variables
% see the mathematica notebook in the derivation folder
% this results in an expected standard deviation in the sample estimated
% standard deviation of sigma* Sqrt[1 - c4^2]


se_samp_std_norm=unbias_normal_std_est_subsamp.*repmat(sqrt(1-normal_correction_c4(repeat_samp).^2),[1,output_val_size]);
std_se_opp_norm=se_samp_std_norm.*mean_like_scaling_factor;

%TODO build argument to unc_wmean for dimension along which to operate so i dont have to have this ugly loop here
% est_se_opp_mean_weighted_norm=nan(output_val_size,1);
% est_se_opp_se_weighted_norm=nan(output_val_size,1);
% lets now find the weighted values of the predicted SE values
[est_se_opp_mean_weighted_norm,est_se_opp_se_weighted_norm]=...
    unc_wmean(unbias_normal_std_est_subsamp,std_se_opp_norm);

% Bryce finished vectorizing up to here will leave for kieran to finish/ merge changes

% if we do not assume normality
% we can use 4th centeral (unbiased) moment over the subsamples 
% tmp
unbias_moments_sub_3 = unbias_moments_sub(:,3,:);
unbias_moments_sub_3=permute(unbias_moments_sub_3,[1,3,2]); % permute the dimensions to get the rid of the singleton dim
var_samp_var_arb=(1./repeat_samp).*(unbias_moments_sub_3-(std_est_subsamp.^4).*((repeat_samp-3)./(repeat_samp-1)));
se_samp_std_arb=(1/2).*sqrt(var_samp_var_arb)./std_est_subsamp;
std_se_opp_arb=se_samp_std_arb.*mean_like_scaling_factor;
% this seems to do very well in my tests

% and then compute the weighted mean
[est_se_opp_mean_weighted_arb,est_se_opp_se_weighted_arb]=unc_wmean(est_se_opp,std_se_opp_arb);



%% save some outputs
% two outputs one out.sampling that gives the details of the sampling process
% the second out.result will contain the useful results est se in full anal ect
% the main output will be 

out.sampling.moments_sub=moments_sub;
out.sampling.mean=mean_sub;
out.sampling.std=std_est_subsamp;
out.sampling.ste=ste_est_subsamp;

out.sampling.std_norm_unbias=unbias_normal_std_est_subsamp;
out.sampling.ste_norm_unbias=ste_est_subsamp_normunbias;

out.sampling.se_std_norm=se_samp_std_norm;
out.sampling.se_std_arb=se_samp_std_arb;
out.sampling.projected_whole_se=est_se_opp;
out.sampling.projected_whole_se_norm_unbias=est_se_opp_norm_unbias;

out.sampling.std_projected_whole_se_norm=std_se_opp_norm;
out.sampling.std_projected_whole_se_arb=std_se_opp_arb;
out.sampling.sample_repeats=repeat_samp;
out.sampling.sample_size=sample_num_vec;

out.results.se_fun_whole_unweighted=est_se_opp_mean_unweighted;
out.results.se_se_fun_whole_unweighted=est_se_opp_se_unweighted;

out.results.se_fun_whole_weighted_arb=est_se_opp_mean_weighted_arb;
out.results.se_se_fun_whole_weighted_arb=est_se_opp_se_weighted_arb;

out.results.se_fun_whole_weighted_norm=est_se_opp_mean_weighted_norm;
out.results.se_se_fun_whole_weighted_norm=est_se_opp_se_weighted_norm;

% set the outputs based on what input options were chosen
if use_weighting
    if norm_weighting
        out.results.se_fun_whole=out.results.se_fun_whole_weighted_norm;
        out.results.se_se_fun_whole=out.results.se_se_fun_whole_weighted_norm;
    else
        out.results.se_fun_whole=out.results.se_fun_whole_weighted_arb;
        out.results.se_se_fun_whole=out.results.se_se_fun_whole_weighted_arb;
    end
else
    out.results.se_fun_whole=est_se_opp_mean_unweighted;
    out.results.se_se_fun_whole=est_se_opp_se_unweighted;
end



%out.opp_frac_est_se(:,1)=sample_frac_vec;
%out.opp_frac_est_se(:,2)=est_se_opp;

%the mean of the estimation function as a function of data subsample
%size tells us about the bais of the estimation function


%% try to fit the dependence of mean est fun (subset) so that we can estimate the bias with sample size
if do_mean_fit && numel(sample_frac_vec)>1 && size(unbias_samp_var,2)<2
    %TODO: a model that is asmyptotic to a value
    modelfun=@(b,x) b(1)+b(2).*x;
    weights=1./(out.sampling.ste.^2);
    not_nan_mask=~isnan(out.sampling.sample_size) & ~isnan(out.sampling.mean) & ~isnan(weights);
    predictor=out.sampling.sample_size(not_nan_mask);
    response=out.sampling.mean(not_nan_mask);
    weights=weights(not_nan_mask);
    weights=weights./nansum(weights);
    beta0=[nanmean(response),0];
    cof_names={'offset','grad'};
    opt = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
    fitobject=fitnlm(predictor,response,modelfun,beta0,...
        'Weights',weights,'options',opt,...
        'CoefficientNames',cof_names);
    out.est_mean_dep_fit=fitobject;
    %%itparam=fitobject.Coefficients;
    %osc_fit.model_coefs(ii,:,:)=[fitparam.Estimate,fitparam.SE];
    sigma_threshold=3;%number of standard deviations away from zero to be signfigant
    is_grad_sig=abs(fitobject.Coefficients.Estimate(2))>fitobject.Coefficients.SE(2)*sigma_threshold;
    if is_grad_sig && verbose>0
        warning(['%s: fit to mean result of est fun on data subset \n'...
                'shows that the gradient with data size is not within %.0f sd of zero \n',...
                'you may have a biased estimator\n'],mfilename,sigma_threshold)
    end
end

if do_plots && numel(sample_frac_vec)>1 && size(unbias_samp_var,2)<2
    if isempty(p.Results.plot_fig_name) || isequal(p.Results.plot_fig_name,'')
        stfig('bootstrap results','add_stack',1);
    else
        stfig(p.Results.plot_fig_name,'add_stack',1);
    end
    clf
    subplot(2,1,1)
    hold on
    legends={};
    
    errorbar(out.sampling.sample_size,out.sampling.projected_whole_se,...
    out.sampling.std_projected_whole_se_arb,'ko',...
    'CapSize',3,...
    'MarkerSize',6,...
    'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1,1,1]*0.3);
    
    xl=xlim(gca);
     legends{end+1}='Est SE';
    
    errorbar(out.sampling.sample_size,out.sampling.projected_whole_se_norm_unbias,...
    out.sampling.std_projected_whole_se_norm,'ro',...
    'CapSize',3,...
    'MarkerSize',6,...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[1,0,0]*0.3);
     legends{end+1}='Est SE normality';
    %
    
    line(xl,[1,1]*est_se_opp_mean_unweighted,'Color','k','LineWidth',2)
     legends{end+1}='mean Est SE';
    line(xl,[1,1]*(est_se_opp_mean_unweighted-est_std_opp_se_unweighted),'Color','m','LineWidth',2)
     legends{end+1}='+std Est SE';
    line(xl,[1,1]*(est_se_opp_mean_unweighted+est_std_opp_se_unweighted),'Color','m','LineWidth',2)
     legends{end+1}='-std Est SE';
    if ~isnan(p.Results.true_dist_se)
        legends=[legends,'true dist SE'];
        line(xl,[1,1]*p.Results.true_dist_se,'Color','r','LineWidth',2)
    end
    if ~isnan(p.Results.true_samp_se)
        legends=[legends,'true Samp SE'];
        line(xl,[1,1]*p.Results.true_samp_se,'Color','b','LineWidth',2)
    end
    if max(xl)>=n_total
        yl=ylim(gca);
        line([1,1]*n_total,yl,'Color',[1,1,1]*0.7,'LineWidth',2);
    %bring the point which has all the data on top of the line so that the error bar can be seen
        chi=get(gca, 'Children');
        set(gca, 'Children',flipud(chi))
        legends{end+1}='total data size';
        legends=fliplr(legends);
    end
    
    legend(legends)
    hold off
    xlabel(sprintf('subsample size (whole data set =%u, vert line)',n_total))
    ylabel('est. SE in estimator function on whole data set')


    subplot(2,1,2)
    hold on
    legends={};
    legends{end+1}='Est mean';
    errorbar(out.sampling.sample_size,out.sampling.mean,out.sampling.ste,'ko',...
    'CapSize',3,...
    'MarkerSize',6,...
    'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1,1,1]*0.3)


    errorbar(out.sampling.sample_size,out.sampling.mean,out.sampling.ste_norm_unbias,'ro',...
    'CapSize',3,...
    'MarkerSize',6,...
    'LineWidth',1.5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[1,0,0]*0.3)
    legends{end+1}='Est mean normality';
    xlabel(sprintf('subsample size (whole data set =%u, vert line)',n_total))
    ylabel('mean est fun of subsample')

    if do_mean_fit
        x_plot_fit=col_vec(linspace(min(out.sampling.sample_size),max(out.sampling.sample_size),1e4));
        [y_plot_fit_val,y_plot_fit_ci]=predict(fitobject,x_plot_fit,'Prediction','curve','Alpha',1-erf(1/sqrt(2)));
        plot(x_plot_fit,y_plot_fit_val,'r')
         legends{end+1}='fit';
        plot(x_plot_fit,y_plot_fit_ci(:,1),'g')
         legends{end+1}='fit+se';
        plot(x_plot_fit,y_plot_fit_ci(:,2),'g')
         legends{end+1}='fit-se';
    end
    
    
    if max(xl)>=n_total
        yl=ylim(gca);
        line([1,1]*n_total,yl,'Color',[1,1,1]*0.7,'LineWidth',2);
    %bring the point which has all the data on top of the line so that the error bar can be seen
        chi=get(gca, 'Children');
        set(gca, 'Children',flipud(chi))
        legends=[legends,'total data size'];
    end
    hold off
    drawnow
    legend(legends)
    %check that the estimated errors look about right
    %hist(abs((boot.se_opp-boot.opp_frac_est_se(:,2))./boot.opp_frac_est_se(:,3)),1e2)
    
    if verbose>4   
        stfig('bootstrap diagnostics','add_stack',1);
        
        error_std_whole_to_mean=out.sampling.projected_whole_se-est_se_opp_mean_unweighted;
        sigma_diff_arb=error_std_whole_to_mean./std_se_opp_arb;
        sigma_diff_norm=error_std_whole_to_mean./std_se_opp_norm;
        fprintf('std of error/est error for arb dist %f , norm dist %f \n',std(sigma_diff_arb),std(sigma_diff_norm))
        subplot(2,2,1)
        plot(out.sampling.sample_size,sigma_diff_arb,'r')
        hold on
        plot(out.sampling.sample_size,sigma_diff_norm,'b')
        hold off
        xlabel(sprintf('subsample size (whole data set =%u)',n_total))
        ylabel('error in est SE/predicted error')
        subplot(2,2,2)
        histogram(sigma_diff_arb,round(numel(sigma_diff_arb)/10),'FaceColor','r')
        xlabel('error in est SE/predicted error')
        hold on
        histogram(sigma_diff_norm,round(numel(sigma_diff_norm)/10),'FaceColor','b')
        hold off
        subplot(2,2,3)
        mean_val=nanmean(out.sampling.mean);
        error_mean=out.sampling.mean-mean_val;
        sigma_diff_val=error_mean./out.sampling.ste;
        plot(out.sampling.sample_size,sigma_diff_val,'r')
        xlabel(sprintf('subsample size (whole data set =%u)',n_total))
        ylabel('error in est val/predicted error')
        subplot(2,2,4)
        histogram(sigma_diff_val,round(numel(sigma_diff_val)/10))
        xlabel('error in est val/predicted error')
        fprintf('std of val error %f \n',std(sigma_diff_val))
    end
     
end
detailed_data_out=[];
if save_multi_out, detailed_data_out.out_cell=out_cell; end
if save_input_data, detailed_data_out.in_cell=in_cell; end
detailed_data_out.sample_num_vec=sample_num_vec;
detailed_data_out.sample_frac_vec=sample_frac_vec;

if verbose>0, fprintf('..Done\n'), end

%
end



function in=coerce_logical(in)
if ~islogical(in)
    if ischar(in)
        in=isequal(in,'true');
    elseif isnumeric(in)
        in=logical(in);
    end
end
end