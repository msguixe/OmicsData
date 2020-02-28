"""statisticstools module
   contains functions for statistical analysis"""
    
def mean(numbers):
    """Calculate the mean of a string"""
    return float(sum(numbers)) / max(len(numbers), 1)

def median(numbers):
    """Calculate the median of a string"""
    length = len(numbers)
    if length < 1:
        return print("The median cannot be calculated for just one value")
    if length % 2 == 1:
        return sorted(numbers)[length//2]
    else:
        return sum(sorted(numbers)[length//2-1:length//2+1])/2.0

def sd(numbers, ddof=0):
    """Calculate the standard deviation of a string"""
    import statistics
    return statistics.stdev(numbers)

def normality_evaluation(x):
    """evaluates the normality of a distribution"""
    from scipy import stats
#   x=stats.norm.rvs(loc=5, scale=3, size=100)
    stats.shapiro(x)
    stats.shapiro(x)[1] > 0.05
    return (str(stats.shapiro(x)[1]))
def ttest_one(x,y):
    """stimates the t-test for one sample"""
    from scipy import stats
    pval=str(stats.ttest_1samp(x,y)[1])
    return  pval
def wilcoxontest_one(x,y):
    """stimates the wilcoxon for one sample"""
    from scipy.stats import wilcoxon
    import numpy as np
    x=np.array(x)
    pval = str(wilcoxon(x-y)[1])
    return  pval


def ttest_equal(x,y):
    """stimates the t-test for two samples with equal variances"""
    from scipy import stats
    pval = str(stats.ttest_ind(x,y, equal_var = True)[1])
    return  pval

def ttest_unequal(x,y):
    """stimates the t-test for two samples with equal variances"""
    from scipy import stats
    pval = str(stats.ttest_ind(x,y, equal_var = False)[1])
    return  pval

def wilcoxon(x,y):
    """stimates the Wilcoxon test for independent samples"""
    from scipy import stats
    pval = str(stats.mannwhitneyu(x,y)[1])
    return  pval

def ttest_related(x,y):
    """stimates the t-test for NON independent samples"""
    from scipy import stats
    pval = str(stats.ttest_rel(x,y)[1])
    return  pval

def wilcoxon_related(x,y):
    """stimates the Wilcoxon test for NON independent samples"""
    from scipy import stats
    pval = str(stats.wilcoxon(x,y)[1])
    return  pval

def bartlett(x,y):
    """stimates the bartlett test for equality of variance"""
    from scipy import stats
    pval = str(stats.bartlett(x,y)[1])
    return  pval


def levene(x,y):
    """stimates the levene test for equality of variance"""
    from scipy import stats
    pval = str(stats.levene(x,y)[1])
    return  pval



def prop_test(x,y):
    """Chi-square test of independence of variables in a contingency table"""
    import numpy as np
    from scipy import stats
    data=[x,y]
    pval=str(stats.chi2_contingency(data)[1])
    return  pval
    

