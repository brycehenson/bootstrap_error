function [mean_val,ste_mean]=unc_wmean(x,unc)

x=col_vec(x);
unc=col_vec(unc);
if numel(unc)~=1 && numel(x)~=numel(unc)
    error('unc must be the same size as data vec or single element')
end
nan_mask=isnan(x) | isnan(unc);
if sum(nan_mask)~=0
    warning('%s: nan found in data will ignore these points\n',mfilename)
end

x=x(~nan_mask);
unc=unc(~nan_mask);

if numel(x)==0
    warning('%s: no data pts remaining after removing nans will retun nans',mfilename)
    ste_mean=nan;
    mean_val=nan;
else
    if numel(unc)==1
        mean_val=mean(x);
        ste_mean=unc/sqrt(numel(x));
    else
        weight=1./(unc.^2);
        [ste_mean,mean_val]=sewm(x,weight);
    end
end

end