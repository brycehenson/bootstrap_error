# boostrap_error
A matlab function that uses bootstaping to find the standard error in an arbitrary analysis operation.

It only takes a moderate amount of complexity in data analysis before it is difficult to determine the error in the result. Boostraping is a powerfull statistical method that performs the analysis repeatedly on smaller subsets of the data in order to *estimate* the error in the final result (using all the data).

**[TO BE CHECKED]** The standard error estimate of the bootraping procedure is reasonably robust against non gaussian distributions and can be used to detect them. 

The procedure is this

# Details
- Finite sample correction from [L. Isserlis,On the Value of a Mean as Calculated from a Sample,J. Royal Stat. Soc
Vol. 81, No. 1 (Jan., 1918), pp. 75-81](http://doi.org/10.2307/2340569)

![fig1](/fig1.png)

# Further Reading
- [wiki](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
- [pi bootstrapped](https://pypi.org/project/bootstrapped/)
- [Jason Brownlee,A gentler introduction to the bootstrap method](https://machinelearningmastery.com/a-gentle-introduction-to-the-bootstrap-method/)


## To Do
- Documentation
- allow second output from anal_opp function (to be used as a structure of details about the fit)
  - allow for this second output not to be provided
- provide all these second outputs in a cell matrix
- test for sub datasets
- test biasing for non gaussian distribution
- write out as tutorial
- normalily testing
