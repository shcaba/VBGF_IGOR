#Starting parameter values
#init value #phase
518     1     #Linf
0.2     2     #k
0       1     #t0

1                           #Number of reads to use

#STANDARD NON-LINEAR FIT
1                           #IF STANDARD NON-LINEAR FIT: Regression (1); Robust Regression (2)
#lower bound    #upper bound    #init value #phase
0     1     0.7     -2
1                                  #IF STANDARD NON-LINEAR FIT: Fit to primary read (1); Fit to Avg. age (2); Fit length vs. age (3); Fit to lengths with all age reads (4)
0.1     -2

-999   #Standard deviation of ages