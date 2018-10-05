function [stats,outs]=boostrap_se(anal_opp,data,varargin)
%funciton calculates the standard error of performing anal_opp(data)
%data is a vector of cells or scalars

p = inputParser;
is_lims=@(x) isequal(size(x),[1,2]) & isnumeric([1,1]) & x(2)<=1 & x(1)>1e-4;
addOptional(p,'replace',false,@logicalable);
addOptional(p,'num_samp_frac',10,@(x) isnumeric(x) & x>=1);
addOptional(p,'num_samp_rep',10,@(x) isnumeric(x) & x>=1);
addOptional(p,'samp_frac_lims',[1e-3,1],is_lims);
addOptional(p,'plots',false,@logicalable);
addOptional(p,'plot_fig_num',10,@(x) isnumeric(x) & x>=1);
parse(p,varargin{:});

do_plots=coerce_logical(p.Results.plots);
do_replace=coerce_logical(p.Results.replace);

%input taken care of

sample_frac_vec=linspace(p.Results.samp_frac_lims(1),...
    p.Results.samp_frac_lims(2),p.Results.num_samp_frac);
repeat_samp_prefactor=p.Results.num_samp_rep;
std_frac=NaN(numel(sample_frac_vec),1);
n_total=numel(data);

fprintf('Bootstaping with different sample fractions %03u:%03u',0)
for ii=1:numel(sample_frac_vec)
    sample_frac=sample_frac_vec(ii);
    n_sample=floor(sample_frac*n_total);
    %std means nothing for n<3
    %the finte sample correaction for the no replacements method breaks when n_sample=ntot
    if n_sample>3 && (n_sample<n_total || do_replace)
        %the scaling of std(std(x)) is 1/n so lets take more data at smaller sample_frac
        repeat_samp_scaled=round(repeat_samp_prefactor*1/sample_frac);
        anal_sample=NaN(repeat_samp_scaled,1);
        if do_replace
            for jj=1:repeat_samp_scaled
                anal_sample(jj)=anal_opp(randsample(data,n_sample,true));
            end
            %setimate the population std
            std_frac(ii)=std(anal_sample)*sqrt(n_sample);
        else
            for jj=1:repeat_samp_scaled
                anal_sample(jj)=anal_opp(randsample(data,n_sample));
            end
            %use finte population correction to estimate the population std using sampling without
            %replacements
            finite_pop_corr=sqrt((n_total-n_sample)/(n_total-1));
            std_frac(ii)=std(anal_sample)*sqrt(n_sample)/finite_pop_corr;
        end
    end
    fprintf('\b\b\b%03u',ii);
end
fprintf('..Done\n')
stats=[];
stats.opp_frac_est_se(:,1)=sample_frac_vec;
stats.opp_frac_est_se(:,2)=std_frac/sqrt(numel(data));
stats.se_opp=nanmean(stats.opp_frac_est_se(:,2));
stats.std_se_opp=nanstd(stats.opp_frac_est_se(:,2));
stats.se_se_opp=stats.std_se_opp/sqrt(sum(~isnan(std_frac)));

if do_plots
    figure(p.Results.plot_fig_num);
    clf
    set(gcf,'color','w')
    plot(stats.opp_frac_est_se(:,1),stats.opp_frac_est_se(:,2),'ko-');
    hold on
    xl=xlim(gca);
    line(xl,[1,1]*stats.se_opp,'Color','k','LineWidth',2)
    line(xl,[1,1]*(stats.se_opp-stats.std_se_opp),'Color','m','LineWidth',2)
    line(xl,[1,1]*(stats.se_opp+stats.std_se_opp),'Color','m','LineWidth',2)
    hold off
    legend('Est SE','mean Est SE','±std Est SE')
    xlabel('fraction of whole data set')
    ylabel('est. SE in opp on whole data set')
    
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