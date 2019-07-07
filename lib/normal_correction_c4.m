function out=normal_correction_c4(n)
%https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
%https://www.jstor.org/stable/2682923?origin=crossref&seq=1#metadata_info_tab_contents

if nargin==0
    error('no input')
end

if n==1
    error('undefined for n=1')
elseif n>1e7
    c4n=1; %gives an error of 2e-7
else
    c4n=sqrt(2./(n-1)).*(gamma(n/2)./gamma((n-1)./2));
    nan_mask=isnan(c4n);
    if sum(nan_mask)>0 %handle numerical overflow and use the vpa toolbox (slower)
        %warning('using vpa')
        for ii=1:numel(c4n)
            n_vpa=vpa(n(nan_mask));
            c4n(nan_mask)=sqrt(2./(n_vpa-1)).*(gamma(n_vpa./2)./gamma((n_vpa-1)./2));
        end
    end
end
out=c4n;
end
