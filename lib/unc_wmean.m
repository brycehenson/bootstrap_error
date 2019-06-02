function [mean_val,ste_mean]=unc_wmean(x,unc)

x=col_vec(x);
unc=col_vec(unc);
if numel(unc)~=1 && numel(x)~=numel(unc)
    error('unc must be the same size as data vec or single element')
end

if numel(unc)==1
    mean_val=mean(x);
    ste_mean=unc/sqrt(numel(x));
else
    weight=1./(unc.^2);
    [ste_mean,mean_val]=sewm(x,weight);
end

end