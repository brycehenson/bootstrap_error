function out=std_c4(x)
%https://en.wikipedia.org/wiki/Unbiased_estimation_of_standard_deviation
%https://www.jstor.org/stable/2682923?origin=crossref&seq=1#metadata_info_tab_contents

x=x(:);
n=numel(x);

out=std(x)/normal_correction_c4(n);
%out=c4n;
end
