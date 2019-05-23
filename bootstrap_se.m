function [out,detailed_data_out]=boostrap_se(est_fun,data,varargin)
%funciton calculates the standard error of performing est_fun(data)
%data is a vector of cells or scalars

%est_fun (estimator opperation) is a function that takes a vector or a cell matrix and retuns a scalar as the first output
%this code can store an arbitraty number of function outputs in a cell array and return them for use
% Known BUGS/ Possible Improvements
%   -handle arbitrary number of scalars (as a vector or matrix) for simulanious multi paramerter boostraping
%   - [ ] simplify codeflow with subfunctions
%   - [ ] add in functionality to look at the mean with the sample fraction
%   - [ ] design a way to schedule in the number of evaluations for each sample fraction
% Author: Bryce Henson
% email: Bryce.Henson[the a swirly]live.com you must have
% '[matlab][bootstrap_error]' in the subject line for me to respond
% Last revision:2018-10-08

min_frac=1e-4;
min_num=5;

p = inputParser;
is_lims=@(x) (isequal(size(x),[1,2]) && (isnumeric(x) && x(2)<10 && x(1)>1e-4))...
        || (isscalar(x) && x<=1 && x>min_frac);
is_c_logical=@(in) isequal(in,true) || isequal(in,false); %can x be cast as a logical
addOptional(p,'replace',false,is_c_logical);
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
parse(p,varargin{:});

do_mean_fit=p.Results.do_mean_fit;
do_plots=coerce_logical(p.Results.plots);
do_replace=coerce_logical(p.Results.replace);
use_mean_se_for_se_se=coerce_logical(p.Results.mean_se_for_se_se);
save_multi_out=coerce_logical(p.Results.save_multi_out);
save_input_data=coerce_logical(p.Results.save_input_data);
verbose=p.Results.verbose;
opp_arguments=p.Results.opp_arguments;
repeat_samp_prefactor=p.Results.num_samp_rep;
%input taken care of

%if the number of fractions to sample is one
if p.Results.num_samp_frac==1
     sample_frac_vec=mean(p.Results.samp_frac_lims);
%sample over multiple fractions
elseif size(p.Results.samp_frac_lims,2)==2 
sample_frac_vec=linspace(p.Results.samp_frac_lims(1),...
    p.Results.samp_frac_lims(2),p.Results.num_samp_frac)';
%if the lims are a scalar then just do one fraction
elseif isscalar(p.Results.samp_frac_lims)
    sample_frac_vec=p.Results.samp_frac_lims;
end
%overwrite plots if only samplign at one fraction of the data
%if size(sample_frac_vec,1)==1
%   do_plots=false;
%end

n_total=numel(data);
sample_num_vec=floor(sample_frac_vec*n_total);
sample_frac_vec=sample_num_vec/n_total;
%cull anything below min_frac or min_num
mask=sample_num_vec>=min_num & sample_frac_vec>min_frac;
sample_num_vec=sample_num_vec(mask);
sample_frac_vec=sample_frac_vec(mask);
iimax=numel(sample_frac_vec);

%catch when there are no valid sample sizes
if iimax>0
    %prealocate the moments of the distribution
    mean_sub=NaN(numel(sample_frac_vec),1);
    moments_sub=NaN(numel(sample_frac_vec),3);
    repeat_samp=NaN(numel(sample_frac_vec),1);

    %find the output size of the passed estimator function
    %should build in optional argument to specify this
    output_size=nargout(est_fun);
    %this will only take the first output of an an anonymous function
    if output_size==-1
        output_size=1;
    end
    
    out_cell=cell(iimax,repeat_samp_prefactor);
    in_cell=out_cell;
    out_cell_tmp=cell(1,output_size);
    if output_size>1 && save_multi_out
        multi_out=true;
    else
        multi_out=false;
    end

    if verbose>0
        fprintf('Bootstrapping with different sample fractions %04u:%04u',0)
    end
    for ii=1:iimax
        n_sample=sample_num_vec(ii);
        %std means nothing for n<3
        %the finte sample correaction for the no replacements method breaks when n_sample=ntot
        if n_sample>3 && (n_sample<n_total || do_replace)
            %this is alowed to vary
            repeat_samp(ii)=round(repeat_samp_prefactor); %*1/sample_frac
            est_fun_res_sub=NaN(repeat_samp(ii),1);
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
                if multi_out
                    [out_cell_tmp{:}]=est_fun(data_smpl,opp_arguments{:});
                    est_fun_res_sub(jj)=out_cell_tmp{1};
                    out_cell{ii,jj}=out_cell_tmp{2:end};
                else
                    est_fun_res_sub(jj)=est_fun(data_smpl,opp_arguments{:});
                end 
                if save_input_data
                    in_cell{ii,jj}=data_smpl;
                end
            end
            mean_sub(ii)=mean(est_fun_res_sub);
            % biased sample variance of the subset
            moments_sub(ii,1)=moment(est_fun_res_sub,2);
            moments_sub(ii,2)=moment(est_fun_res_sub,3);
            moments_sub(ii,3)=moment(est_fun_res_sub,4);

            %use finte population correction to estimate the population std using sampling without
            %replacements, if do_replace finite_pop_corr=1
            moments_sub(ii,1)=moments_sub(ii,1)/finite_pop_corr;
        end
        if verbose>0, fprintf('\b\b\b\b%04u',ii), end
    end
    
    out=[];

    
    %unbiased sample variance:
    unbias_samp_var=(repeat_samp./(repeat_samp-1)).* moments_sub(:,1);
    std_est_subsamp=sqrt(unbias_samp_var); %TODO calulate unbiased sample standard deviation
    
    %now calulate the corrected sample standard deviation in the anal
    %operation on the whole datset
    est_se_opp=std_est_subsamp...
        .*sqrt(sample_num_vec)./sqrt(n_total);
    out.moments_sub=moments_sub;
    out.opp_frac_est_se(:,1)=sample_frac_vec;
    out.opp_frac_est_se(:,2)=est_se_opp;

    %the mean of the estimation function as a function of data subsample
    %size tells us about the bais of the estimation function
    out.sampled_est.mean.val=mean_sub;
    out.sampled_est.mean.se=std_est_subsamp./sqrt(repeat_samp);
    out.sampled_est.samp_num=sample_num_vec;
    out.sampled_est.samp_repeats=sample_num_vec;

    
    %calculate the estimated sample variance of the estimated sample variance
    if use_mean_se_for_se_se
        est_var_opp=nanmean(unbias_samp_var);
    else
        est_var_opp=abs(moments_sub(:,1));
    end
    est_var_var_opp=(1./sample_num_vec).*(moments_sub(:,3)-...  
                ((sample_num_vec-3)./(sample_num_vec-1)).*(moments_sub(:,1).^2));
    est_std_std_opp=est_se_opp.*0.5.*est_var_var_opp./est_var_opp;
    %i have no idea where these come from
    est_std_std_opp=est_std_std_opp.*(sample_num_vec.^(3.5))...
        .*(repeat_samp.^(-0.5)).*(n_total^(-1));

    out.opp_frac_est_se(:,3)=est_std_std_opp;
    mask=sum(isnan(out.opp_frac_est_se),2)==0 & sum(isinf(out.opp_frac_est_se),2)==0;
    out.opp_frac_est_se=out.opp_frac_est_se(mask,:);
    
    
    out.se_opp_unweighted=nanmean(out.opp_frac_est_se(:,2));
    %calucalte the weighted mean of the estimates SE
    out.se_opp=nansum(out.opp_frac_est_se(:,2)...
                    ./(out.opp_frac_est_se(:,3).^2))./...
                nansum(1./(out.opp_frac_est_se(:,3).^2));
    out.se_se_opp=sqrt(1./sum((out.opp_frac_est_se(:,3).^-2)));
    out.std_se_opp=out.se_se_opp.*sqrt(sum(~isnan(out.opp_frac_est_se(:,2))));
    if numel(out.opp_frac_est_se(:,2))==1
        out.std_se_opp_unweighted=nan;
    else
        out.std_se_opp_unweighted=nanstd(out.opp_frac_est_se(:,2));
    end
    out.se_se_opp_unweighted=out.std_se_opp_unweighted...
        /sqrt(sum(~isnan(out.opp_frac_est_se(:,2))));

    %% try to fit the dependence of mean est fun (subset) so that we can estimate the bias
   if do_mean_fit
    modelfun=@(b,x) b(1)+b(2).*x;
    weights=1./(out.sampled_est.mean.se.^2);
    weights=weights./sum(weights);
    
    beta0=[mean(out.sampled_est.mean.val),0];
    cof_names={'offset','grad'};
    opt = statset('TolFun',1e-10,'TolX',1e-10,...
            'MaxIter',1e4,... %1e4
            'UseParallel',1);
    fitobject=fitnlm(out.sampled_est.samp_num,out.sampled_est.mean.val,modelfun,beta0,...
        'Weights',weights,'options',opt,...
        'CoefficientNames',cof_names);
    out.est_mean_dep_fit=fitobject;
    %%itparam=fitobject.Coefficients;
    %osc_fit.model_coefs(ii,:,:)=[fitparam.Estimate,fitparam.SE];
    is_grad_sig=fitobject.Coefficients.pValue(2)<0.05;
    if ~is_grad_sig
        warning('%s: warning fit to mean of set fun subset shows that the gradient is not withing error of zero',mfilename)
    end
   end
    
    if verbose>0, fprintf('..Done\n'), end
    if do_plots
        if isempty(p.Results.plot_fig_name) || isequal(p.Results.plot_fig_name,'')
            stfig('bootstrap results','add_stack',1);
        else
            stfig(p.Results.plot_fig_name,'add_stack',1);
        end
        subplot(2,1,1)
        errorbar(out.sampled_est.samp_num,out.opp_frac_est_se(:,2),...
        out.opp_frac_est_se(:,3),'ko',...
        'CapSize',0,...
        'MarkerSize',6,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1,1,1]*0.2)
        hold on
        xl=xlim(gca);
        yl=ylim(gca);
        line(xl,[1,1]*out.se_opp,'Color','k','LineWidth',2)
        line(xl,[1,1]*(out.se_opp-out.std_se_opp_unweighted),'Color','m','LineWidth',2)
        line(xl,[1,1]*(out.se_opp+out.std_se_opp_unweighted),'Color','m','LineWidth',2)
        line([1,1]*n_total,yl,'Color','k','LineWidth',2)
        
        legends={'Est SE','mean Est SE','+std Est SE','-std Est SE','data size'};
        if ~isnan(p.Results.true_dist_se)
            legends=[legends,'true dist SE'];
            line(xl,[1,1]*p.Results.true_dist_se,'Color','r','LineWidth',2)
        end
        if ~isnan(p.Results.true_samp_se)
            legends=[legends,'true Samp SE'];
            line(xl,[1,1]*p.Results.true_samp_se,'Color','b','LineWidth',2)
        end
        hold off
        legend(legends)
        xlabel('subsample fraction of whole data set')
        ylabel('est. SE in estimator function on whole data set')
        drawnow
        subplot(2,1,2)
        errorbar(out.sampled_est.samp_num,out.sampled_est.mean.val,out.sampled_est.mean.se,'ko',...
        'CapSize',0,...
        'MarkerSize',6,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1,1,1]*0.2)
        xlabel('subsample fraction of whole data set')
        ylabel('mean est fun of subsample')
        if do_mean_fit
            hold on
            x_plot_fit=col_vec(linspace(min(out.sampled_est.samp_num),max(out.sampled_est.samp_num),1e4));
            [y_plot_fit_val,y_plot_fit_ci]=predict(fitobject,x_plot_fit);
            plot(x_plot_fit,y_plot_fit_val,'r')
            plot(x_plot_fit,y_plot_fit_ci(:,1),'g')
            plot(x_plot_fit,y_plot_fit_ci(:,2),'g')
            hold off
        end
        
        
        yl=ylim(gca);
        line([1,1]*n_total,yl,'Color','k','LineWidth',2)

        %check that the estimated errors look about right
        %hist(abs((boot.se_opp-boot.opp_frac_est_se(:,2))./boot.opp_frac_est_se(:,3)),1e2)
    end
    detailed_data_out=[];
    if save_multi_out, detailed_data_out.out_cell=out_cell; end
    if save_input_data, detailed_data_out.in_cell=in_cell; end
    detailed_data_out.sample_num_vec=sample_num_vec;
    detailed_data_out.sample_frac_vec=sample_frac_vec;
else %if there are no valid sample sizes then set the output to nans
    out.opp_frac_est_se=nan;
    out.se_opp_unweighted=nan;
    %calucalte the weighted mean of the estimates SE
    out.se_opp=nan;
    out.se_se_opp=nan;
    out.std_se_opp=nan;
    detailed_data_out=[];
end

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