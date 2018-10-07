# boostrap_error
A matlab function that uses bootstaping to find the standard error in an arbitrary analysis operation.

It only takes a moderate amount of complexity in data analysis before it is difficult to determine the error in the result. Boostraping is a powerfull statistical method that performs the analysis repeatedly on smaller subsets of the data in order to *estimate* the error in the final result (using all the data). Further the method is able to work with an analysis operation that only produces meaningfull results when performed with many data points (such as a linear fit)


**[TO BE CHECKED]** The standard error estimate of the bootraping procedure is reasonably robust against non gaussian distributions and can be used to detect them.  However the error in the sfinite sample.
tandard error is an inherently biased estimator and any nongaussian distribution will mean that the estimated value will change (decrease?) with a 


**[TO BE CHECKED]** The method assumes that the data has no correlations/ that the data set are independent samples. What would correlations do here???
The procedure is reasonably simple given some analysis operation A(x) where x is the input
1. select a random sample of the data S  of length n_samp (with replacements) out of all data collected (D, with length n_tot)
2. compute the analysis operation A(S)
3. repeat steps 1 to three many times saving the result of each analysis operation (on the subset)
4. calulate the standard deviation across these results and multiply by sqrt(n_samp)/sqrt(n_tot) to estimate the standard error in A(D).


The above uses random sampling with replacement in order to prevent biasing of the standrd error estimate. It is however possible to use the Finite sample correction from [L. Isserlis,On the Value of a Mean as Calculated from a Sample,J. Royal Stat. Soc
Vol. 81, No. 1 (Jan., 1918), pp. 75-81](http://doi.org/10.2307/2340569) to correct for the bias when using random sampling without replacements. Both methods are implemented in this work.


As a test it is advisable to check that there is no trend in the estimated standard error as function of the size of the subsample. Thus the above procedure is repeated at many different fractions of the whole dataset. 


# Details
- Finite sample correction from [L. Isserlis,On the Value of a Mean as Calculated from a Sample,J. Royal Stat. Soc
Vol. 81, No. 1 (Jan., 1918), pp. 75-81](http://doi.org/10.2307/2340569)

![fig1](/fig1.png)

# Further Reading
- [wiki](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
- [pi bootstrapped](https://pypi.org/project/bootstrapped/)
- [Jason Brownlee,A gentler introduction to the bootstrap method](https://machinelearningmastery.com/a-gentle-introduction-to-the-bootstrap-method/)
- [https://ijpam.eu/contents/2005-21-3/10/10.pdf]
- [https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance]
- [https://ijpam.eu/contents/2009-52-1/5/5.pdf]


## To Do
- Error in the SE estimate [dist. of sample var.](https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance) and [centeral moments](https://en.wikipedia.org/wiki/Central_moment)
- Documentation
- allow second output from anal_opp function (to be used as a structure of details about the fit)
  - allow for this second output not to be provided
- provide all these second outputs in a cell matrix
- test for sub datasets
- test biasing for non gaussian distribution
- write out as tutorial
- normalily testing

