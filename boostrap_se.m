function [stats,outs]=boostrap_se(anal_opp,data,varargin)
%funciton calculates the standard error of performing anal_opp(data)
%data is a vector of cells or scalars

%note the moment based error estimate will be downwards biased if
%num_samp_frac is too small

p = inputParser;
is_lims=@(x) isequal(size(x),[1,2]) & isnumeric([1,1]) & x(2)<=1 & x(1)>1e-4;
addOptional(p,'replace',false,@logicalable);
addOptional(p,'num_samp_frac',10,@(x) isnumeric(x) & x>=1);
addOptional(p,'num_samp_rep',10,@(x) isnumeric(x) & x>=1);
addOptional(p,'samp_frac_lims',[1e-3,1],is_lims);
addOptional(p,'plots',false,@logicalable);
addOptional(p,'plot_fig_num',10,@(x) isnumeric(x) & x>=1);
%optional arguments for diagnostic/test plots
addOptional(p,'true_dist_se',nan,@(x) isnumeric(x) & x>0);
addOptional(p,'true_samp_se',nan,@(x) isnumeric(x) & x>0);
addOptional(p,'mean_se_for_se_se',nan,@logicalable);
parse(p,varargin{:});

do_plots=coerce_logical(p.Results.plots);
do_replace=coerce_logical(p.Results.replace);
use_mean_se_for_se_se=coerce_logical(p.Results.mean_se_for_se_se);

%input taken care of

sample_frac_vec=linspace(p.Results.samp_frac_lims(1),...
    p.Results.samp_frac_lims(2),p.Results.num_samp_frac)';
repeat_samp_prefactor=p.Results.num_samp_rep;
%prealocate the moments of the distribution
moments_sub=NaN(numel(sample_frac_vec),3);
repeat_samp=NaN(numel(sample_frac_vec),1);
n_total=numel(data);
sample_num_vec=floor(sample_frac_vec*n_total);
sample_frac_vec=sample_num_vec/n_total;

fprintf('Bootstaping with different sample fractions %03u:%03u',0)
for ii=1:numel(sample_frac_vec)
    n_sample=sample_num_vec(ii);
    %std means nothing for n<3
    %the finte sample correaction for the no replacements method breaks when n_sample=ntot
    if n_sample>3 && (n_sample<n_total || do_replace)
        %the scaling of std(std(x)) is 1/n so lets take more data at smaller sample_frac
        repeat_samp(ii)=round(repeat_samp_prefactor); %*1/sample_frac
        anal_opp_sub=NaN(repeat_samp(ii),1);
        if do_replace
            for jj=1:repeat_samp(ii)
                %calculate the analysis operation on the subset of dat
                anal_opp_sub(jj)=anal_opp(randsample(data,n_sample,true));
            end
            %setimate the population std
            % biased sample variance of the subset
            moments_sub(ii,1)=moment(anal_opp_sub,2);
            moments_sub(ii,2)=moment(anal_opp_sub,3);
            moments_sub(ii,3)=moment(anal_opp_sub,4);
        else
            for jj=1:repeat_samp(ii)
                anal_opp_sub(jj)=anal_opp(randsample(data,n_sample));
            end
            %use finte population correction to estimate the population std using sampling without
            %replacements
            finite_pop_corr=sqrt((n_total-n_sample)/(n_total-1));
            moments_sub(ii)=std(anal_opp_sub)*sqrt(n_sample)/finite_pop_corr;
        end
    end
    fprintf('\b\b\b%03u',ii);
end

%unbiased sample variance:
unbias_samp_var=(repeat_samp./(repeat_samp-1)).* moments_sub(:,1);
%now calulate the corrected sample standard deviation in the anal
%operation on the whole datset
est_se_opp=sqrt(unbias_samp_var)...
    .*sqrt(sample_num_vec)./sqrt(n_total);
stats=[];
stats.moments_sub=moments_sub;
stats.opp_frac_est_se(:,1)=sample_frac_vec;
stats.opp_frac_est_se(:,2)=est_se_opp;
%calculate the estimated sample variance of the estimated sample variance
if use_mean_se_for_se_se
    est_var=nanmean(unbias_samp_var);
else
    est_var=abs(moments_sub(:,1));
end
est_var_std_opp=(1./sample_num_vec).*...
        sqrt(...
            moments_sub(:,3)-...
            ((sample_num_vec-3)./(sample_num_vec-1)).*(moments_sub(:,1).^2)...
            )...
        -2.*est_var;
est_std_std_opp=sqrt(abs(est_var_std_opp))./sqrt(n_total) ; %sqrt(sample_num_vec)



stats.opp_frac_est_se(:,3)=est_std_std_opp;
mask=sum(isnan(stats.opp_frac_est_se),2)==0 & sum(isinf(stats.opp_frac_est_se),2)==0;
stats.opp_frac_est_se=stats.opp_frac_est_se(mask,:);
stats.se_opp_unweighted=nanmean(stats.opp_frac_est_se(:,2));
%calucalte the weighted mean of the estimates SE
stats.se_opp=nansum(stats.opp_frac_est_se(:,2)...
                ./(stats.opp_frac_est_se(:,3).^2))./...
            nansum(1./(stats.opp_frac_est_se(:,3).^2));
stats.se_se_opp=sqrt(1./sum(1./(stats.opp_frac_est_se(:,3).^2)));
stats.std_se_opp=stats.se_se_opp.*sqrt(sum(~isnan(stats.opp_frac_est_se(:,2))));

stats.std_se_opp_unweighted=nanstd(stats.opp_frac_est_se(:,2));
stats.se_se_opp_unweighted=stats.std_se_opp/sqrt(sum(~isnan(stats.opp_frac_est_se(:,2))));

fprintf('..Done\n')
if do_plots
    figure(p.Results.plot_fig_num);
    clf
    set(gcf,'color','w')
    errorbar(stats.opp_frac_est_se(:,1),stats.opp_frac_est_se(:,2),...
    stats.opp_frac_est_se(:,3),'ko',...
    'CapSize',0,...
    'MarkerSize',6,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1,1,1]*0.2)
    hold on
    xl=xlim(gca);
    line(xl,[1,1]*stats.se_opp,'Color','k','LineWidth',2)
    line(xl,[1,1]*(stats.se_opp-stats.std_se_opp_unweighted),'Color','m','LineWidth',2)
    line(xl,[1,1]*(stats.se_opp+stats.std_se_opp_unweighted),'Color','m','LineWidth',2)
    legends={'Est SE','mean Est SE','+std Est SE','-std Est SE'};
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
    xlabel('fraction of whole data set')
    ylabel('est. SE in opp on whole data set')
    
    %check that the estimated errors look about right
    %hist(abs((boot.se_opp-boot.opp_frac_est_se(:,2))./boot.opp_frac_est_se(:,3)),1e2)
end


%
end

function good=logicalable(x)
good=isequal(x,'true')|| isequal(x,'false') || islogical(x) || isequal(x,0) || isequal(x,1);
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