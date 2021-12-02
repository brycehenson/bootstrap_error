
***[Bryce M. Henson](https://github.com/brycehenson), [Dong K. Shin](https://github.com/spicydonkey), [Kieran F. Thomas](https://github.com/KF-Thomas)***   

**[![Build Status](https://img.shields.io/static/v1.svg?label=CSL&message=software%20against%20climate%20change&color=green?style=flat&logo=github)](https://img.shields.io/static/v1.svg?label=CSL&message=software%20against%20climate%20change&color=green?style=flat&logo=github)
**

A matlab function that uses bootstrapping to find the standard error in an arbitrary analysis operation.
**Status:** This core functionality provided here  **is ready for use in other projects**. Testing is implemented and passing for the core functionality which provides error determination.

It only takes a moderate amount of complexity in data analysis operation(estimation function) before it is difficult to determine the error in the result. Bootstraping/resampling is a powerful statistical method that performs the analysis operation repeatedly on smaller subsets of the data in order to *estimate* the error in the result of the the operation(estimation function) on the full data set. Further the method is able to work with an analysis operation that only produces meaningful results when performed with many data points (such as a (non)linear fit)

The procedure is reasonably simple given some analysis operation (estimation function) **A(x)** (that produces a scalar) and a dataset **D**
1. select a random sample of the data **S** of length **n_samp** (with replacements) out of all data collected (D, with length n_tot)
2. compute the analysis operation **A(S)**
3. repeat steps 1 to 2 many times saving the result of each analysis operation (on the subset)
4. calculate the standard deviation across these results and multiply by **sqrt(n_samp)/sqrt(n_tot)** to estimate the standard error in **A(D)**. (This is known as mean-like scaling.)




As a test it is advisable to check that there is no trend in either: the output of the function, or the estimated standard error, as function of the size of the subset. Thus the above procedure may be repeated at many different fractions of the whole dataset giving the graph below. 

| ![Figure 1](/figs/fig1.png "Fig1") | 
|:--:| 
 **Figure 1**- Bias analysis output. This graph can be used to reveal the bias of the estimation function with sample size and how the error in the result scales with sample size. An estimation function is mean-like if the estimated SE in the operation on the whole data set does not change with subsample fraction. |


## Features
### sampling without replacement
The above uses random sampling with replacement in order to prevent biasing of the standard error estimate. It is however possible to use the Finite sample correction from [L. Isserlis,On the Value of a Mean as Calculated from a Sample,J. Royal Stat. Soc
Vol. 81, No. 1 (Jan., 1918), pp. 75-81](http://doi.org/10.2307/2340569) to correct for the bias when using random sampling without replacements. Both methods are implemented in this work.
### Estimated error in the error
If you are studying the error from some analysis operation then it may be required to know how significant some change in this error is. This is where it is natural to start worrying about the error in the SE. This code provides two estimates, the first assumes a normal distribution (and is unbiased) the second does not and is slightly biased.



## Further Reading
- [wiki](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
- [pi bootstrapped](https://pypi.org/project/bootstrapped/)
- [Jason Brownlee,A gentler introduction to the bootstrap method](https://machinelearningmastery.com/a-gentle-introduction-to-the-bootstrap-method/)
- [https://ijpam.eu/contents/2005-21-3/10/10.pdf]
- [https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance]
- [https://ijpam.eu/contents/2009-52-1/5/5.pdf]
- [dist. of sample var.](https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance) 
- [centeral moments](https://en.wikipedia.org/wiki/Central_moment)
- see derivation folder
- https://ijpam.eu/contents/2005-21-3/10/10.pdf
- https://ijpam.eu/contents/2009-52-1/5/5.pdf

## To Do
contributors welcome! There is a lot to do to build this into a powerful tool. Drop me an email. 
- [x] allow for vctor output of function, producing the estimated error for each element of the vector
- [x] add an option to let the estimator function handle the subsampling
- [x] fix the overestimation of error at small fractions of the dataset
- [x] allow second output from anal_opp function (to be used as a structure of details about the fit)
  - [x] allow for this second output not to be provided
  - [x] provide all these second outputs in a cell matrix
- [ ] error in the error
- [ ] fix how the estimated SE of the SE estimate is calculated
  - [ ] make all work for even a single subsample size
  - [x] check estimated SE of the SE estimate by nesting the bootstrapper
  - [ ] undestand why results are wrong
  - [ ] understand how should treat combining multiple sampling fractions
- [ ] fit a [laurent series](https://en.wikipedia.org/wiki/Laurent_series) to the mean and error dependence
- [ ] more Documentation
  - [ ] careful documentation in function of what each output is
  - [ ] commenting in main function with links
  - [ ] organizing the resources in this readme
- [ ] write tutorial in more detail
- [ ] normalily testing
  - during the bootstrap try and determine if the underlying distribution is normal or not
  - [package](https://au.mathworks.com/matlabcentral/fileexchange/60147-normality-test-package)
  - build convergence test for some distributions (normal,uniform,[arb. kurtosis](https://en.wikipedia.org/wiki/Kurtosis#The_Pearson_type_VII_family)
- [ ] make a nice logo/diagram
- [ ] add to matlab file exchange
  



