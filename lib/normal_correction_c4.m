function out=normal_correction_c4(n)
%https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
%https://www.jstor.org/stable/2682923?origin=crossref&seq=1#metadata_info_tab_contents

if n>1e7
    c4n=1; %gives an error of 2e-7
else
    c4n=sqrt(2./(n-1)).*(gamma(n/2)./gamma((n-1)./2));
    if isnan(c4n)%handle numerical overflow and use the vpa toolbox (slower)
        %warning('using vpa')
        n=vpa(n);
        c4n=sqrt(2./(n-1)).*(gamma(n./2)./gamma((n-1)./2));
    end
end
out=c4n;
end
