# bootstrap_error
**Bryce Henson, Dong K. Shin** 
Matlab code for fast masking/selection of ordered vectors based on binary search.
A matlab function that uses bootstrapping to find the standard error in an arbitrary analysis operation.

It only takes a moderate amount of complexity in data analysis before it is difficult to determine the error in the result. Bootstrapping is a powerfull statistical method that performs the analysis repeatedly on smaller subsets of the data in order to *estimate* the error in the final result (using all the data). Further the method is able to work with an analysis operation that only produces meaningfull results when performed with many data points (such as a linear fit)


**[TO BE CHECKED]** The standard error estimate of the bootsrapping procedure is reasonably robust against non gaussian distributions and can be used to detect them.  However the error in the sample
standard error is an inherently biased estimator and any nongaussian distribution will mean that the estimated value will change (decrease?) from the true value with a finite sample size.
**[TO BE CHECKED]** The method also assumes that the data has no correlations/ that the data set are independent samples. What would correlations do here???


The procedure is reasonably simple given some analysis operation A(x) where x is the input
1. select a random sample of the data S  of length n_samp (with replacements) out of all data collected (D, with length n_tot)
2. compute the analysis operation A(S)
3. repeat steps 1 to three many times saving the result of each analysis operation (on the subset)
4. calulate the standard deviation across these results and multiply by sqrt(n_samp)/sqrt(n_tot) to estimate the standard error in A(D).


The above uses random sampling with replacement in order to prevent biasing of the standard error estimate. It is however possible to use the Finite sample correction from [L. Isserlis,On the Value of a Mean as Calculated from a Sample,J. Royal Stat. Soc
Vol. 81, No. 1 (Jan., 1918), pp. 75-81](http://doi.org/10.2307/2340569) to correct for the bias when using random sampling without replacements. Both methods are implemented in this work.


As a test it is advisable to check that there is no trend in the estimated standard error as function of the size of the subsample. Thus the above procedure is repeated at many different fractions of the whole dataset giving the graph below. 
![fig1](/fig1.png)

## Error in the estimated error
- Error in the SE estimate [dist. of sample var.](https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance) and [centeral moments](https://en.wikipedia.org/wiki/Central_moment)
-see derivation folder

## Details
- Finite sample correction from [L. Isserlis,On the Value of a Mean as Calculated from a Sample,J. Royal Stat. Soc
Vol. 81, No. 1 (Jan., 1918), pp. 75-81](http://doi.org/10.2307/2340569)



## Further Reading
- [wiki](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
- [pi bootstrapped](https://pypi.org/project/bootstrapped/)
- [Jason Brownlee,A gentler introduction to the bootstrap method](https://machinelearningmastery.com/a-gentle-introduction-to-the-bootstrap-method/)
- [https://ijpam.eu/contents/2005-21-3/10/10.pdf]
- [https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance]
- [https://ijpam.eu/contents/2009-52-1/5/5.pdf]


## To Do
contributors welcome! There is a lot to do to build this into a powerful tool. Drop me an email. 
- build convergence test for some distributions (normal,uniform,(arb. kurtosis)[https://en.wikipedia.org/wiki/Kurtosis#The_Pearson_type_VII_family]
- check estimated SE of the SE estimate by nesting the bootstrapper
- fix the overestimation of error at small fractions of the dataset
- more Documentation
  - commenting in main function with links
  - organizing the resources in this readme
- allow second output from anal_opp function (to be used as a structure of details about the fit)
  - allow for this second output not to be provided
  - provide all these second outputs in a cell matrix
- write tutorial in more detail
- normalily testing
  - during the bootstrap try and determine if the underlying distribution is normal or not
  - (package)[https://au.mathworks.com/matlabcentral/fileexchange/60147-normality-test-package]



