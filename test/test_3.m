%what does it do with a cauchy distribution???
%generate radom numbers drawn from a cauchy dist
%https://au.mathworks.com/matlabcentral/answers/243408-how-to-generate-samples-from-a-cauchy-distribution
cvar = @(N) tan((rand(1,N) - 0.5)*pi);
x=cvar(100000000);
histogram(x,linspace(-20,20,1e3))