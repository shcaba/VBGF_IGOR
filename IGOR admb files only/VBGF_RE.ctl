#Starting parameter values
#init value #phase
912     1    #Linf
0.04     2          #k
0     1        #t0

3                           #Number of reads to use
0.150063212712715                             #CV for the random effects error term
2                           #Age likelihood: Normal (1), Exponential (2), Gamma (3)
2                            #Random effects phase

#init value #phase
0.1     2   #Variance on epsilon error term

#If NORMAL
#init value #phase
5     -2     #Mean value of sampling distribution
1     -2     ##Variance of sampling distribution

#If EXPONENTIAL
#init value #phase
0.2     2         #Rate change of the exponential (Z)

#If GAMMA
#init value #phase
5     -2         #Shape parameter
3     -2         #Scale parameter