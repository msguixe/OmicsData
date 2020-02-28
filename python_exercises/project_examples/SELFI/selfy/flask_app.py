#################################################################
#
# SELFY 1.0 SIMPLE STATISTICAL APP FOR EVERYONE
#
#################################################################

#/usr/bin/env python3
#Python version: Python 3.6.4
#Operative sistem :Ubuntu 16.04.4
#Author Ruth Barral-Arca
#####################################################################################################################
'''-SELFY is a Web app that allow to perform the most common statistical tests and figures.
   -SELFY is designed thinking about people that are neither confortable with R nor  statistics , but wish to perform
    quick and easy statistical analysis. For example check if the value of a citokine is related with the severity level
    and plot it in a boxplot
   -SELFY Input is data directly copied from an excel file column.WARNING!! The decimal symbol must be dot NOT COMMA
   -SELFY output :p-values, plots ,r-square etc
'''
######################################################################################################################
from flask import Flask, render_template, request
from flask import request
from werkzeug import secure_filename
import statisticstools

# My app 
app = Flask(__name__)

@app.route('/')
def main():
    return render_template('index.html')

#################################################################
#
# An analysis of one continuous variable
#
#################################################################

@app.route('/continuous.html')
def transform():
    return render_template('continuous.html')

@app.route('/continuous_results', methods=['POST'])
def continuous_results():
#statistical analysis of a continuous variable
    input_variable =request.form['input_variable']
    input_variable = input_variable.splitlines()
    input_variable = [float(i) for i in input_variable]
    transform = request.form['continuous']

    if transform =='mean': #calculates the mean
        output= statisticstools.mean( input_variable)
    elif transform == 'median': #calculates the median
        output = statisticstools.median(input_variable)
    elif transform == 'sd':#calculates the standard deviation
        output = statisticstools.sd(input_variable)
    elif transform == 'iqr': #calculates the interquartile range
        import numpy as np
        iqr = np.subtract(*np.percentile(input_variable, [75, 25]))
        output = iqr

    return render_template('continuous_results.html', **locals())

#################################################################
#
# correlation
#
#################################################################
@app.route('/twocontinuous.html')
def transform8():
    return render_template('twocontinuous.html')

@app.route('/twocontinuous_results', methods=['POST'])
def twocontinuous_results():
#calculates the correlation coeficient
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]
    input_variable2 =request.form['input_variable2']
    input_variable2 = input_variable2.splitlines()
    input_variable2= [float(i) for i in input_variable2]
    transform8 = request.form['twocontinuous']
    if transform8 =="pearson":
        from scipy.stats.stats import pearsonr
        pearsonr(input_variable1,input_variable2)
        output= str(pearsonr(input_variable1,input_variable2)[0])
        pval= str(pearsonr(input_variable1,input_variable2)[1])
    elif transform8 =="spearman":
        from scipy.stats.stats import spearmanr
        spearmanr(input_variable1,input_variable2)
        output= str(spearmanr(input_variable1,input_variable2)[0])
        pval= str(spearmanr(input_variable1,input_variable2)[1])
    print(output) #works 
    return render_template('twocontinuous_results.html', **locals())
    
#################################################################
#
# plots
#
#################################################################

##scatterplot
@app.route('/dotplots.html')
def transform2():
    return render_template('dotplots.html')

@app.route('/dotplots_result', methods=['POST'])
def dotplot_results():
#creates a scatterplot
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    input_variable2 =request.form['input_variable2']
    input_variable2 = input_variable2.splitlines()
    input_variable2= [float(i) for i in input_variable2]

    import matplotlib.pyplot as plt
    import io
    import base64 
    img = io.BytesIO()
    plt.scatter(input_variable1,input_variable2)
    plt.savefig(img, format='png')
    img.seek(0)
    plt.close() 
    plot_url = base64.b64encode(img.getvalue()).decode()
    return'<img src="data:image/png;base64,{}">'.format(plot_url)

##Histogram
@app.route('/histogram.html')
def transform3():
    return render_template('histogram.html')

@app.route('/histogram_result', methods=['POST'])
def histogram_results():
#creates an histogram
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    import seaborn as sns
    import numpy as np
    import pandas as pd
    from scipy import stats, integrate
    import matplotlib.pyplot as plt
    from pylab import savefig
    import io
    import base64 
    img = io.BytesIO()   
    svm=sns.distplot(input_variable1, kde=False, rug=True)
    figure = svm.get_figure()    
    figure.savefig(img, format='png')
    figure.clf()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    return'<img src="data:image/png;base64,{}">'.format(plot_url)

##Boxplot
@app.route('/Boxplot.html')
def transform4():
    return render_template('Boxplot.html')

@app.route('/boxplot_result', methods=['POST'])
def boxplot_results():
#creates a boxplot of two variables
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    input_variable2 =request.form['input_variable2']
    input_variable2 = input_variable2.splitlines()
    input_variable2= [float(i) for i in input_variable2]
    data = [input_variable1, input_variable2]

    import matplotlib.pyplot as plt
    import io
    import base64 
    img = io.BytesIO()
    
    plt.boxplot(data)
    plt.savefig(img, format='png')
    plt.close() 
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    return'<img src="data:image/png;base64,{}">'.format(plot_url)

##Boxplot
@app.route('/boxplot_one.html')
def boxpot_one():
    return render_template('boxplot_one.html')
@app.route('/boxplot_one_results', methods=['POST'])
def boxplot_results_one():
#creates a boxplot of one variable
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    import matplotlib.pyplot as plt
    import io
    import base64 
    img = io.BytesIO()
    
    plt.boxplot(input_variable1)
    plt.savefig(img, format='png')
    plt.close() 
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    return'<img src="data:image/png;base64,{}">'.format(plot_url)

##Boxplot for six variables
@app.route('/Boxplot6.html')
def boxplotsix():
    return render_template('Boxplot6.html')

@app.route('/boxplotsix_results', methods=['POST'])
def box6_results():
#creates a boxplot of 6 variables (it was an special request of a selfy user)
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    input_variable2 =request.form['input_variable2']
    input_variable2 = input_variable2.splitlines()
    input_variable2= [float(i) for i in input_variable2]

    input_variable3 =request.form['input_variable3']
    input_variable3 = input_variable3.splitlines()
    input_variable3 = [float(i) for i in input_variable3]

    input_variable4 =request.form['input_variable4']
    input_variable4 = input_variable4.splitlines()
    input_variable4= [float(i) for i in input_variable4]

    input_variable5 =request.form['input_variable5']
    input_variable5 = input_variable5.splitlines()
    input_variable5= [float(i) for i in input_variable5]
   
    input_variable6 =request.form['input_variable6']
    input_variable6 = input_variable6.splitlines()
    input_variable6= [float(i) for i in input_variable6]

    data = [input_variable1,input_variable2,input_variable3,input_variable4,input_variable5, input_variable6]
    import matplotlib.pyplot as plt
    import io
    import base64 
    img = io.BytesIO()
    
    plt.boxplot(data)
    plt.savefig(img, format='png')
    plt.close() 
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    return'<img src="data:image/png;base64,{}">'.format(plot_url)

##Density
@app.route('/Density.html')
def density_one():
    return render_template('Density.html')
@app.route('/Density_results', methods=['POST'])
def Density_results():
#creates a density plot
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    import seaborn as sns
    import numpy as np
    import pandas as pd
    from scipy import stats, integrate
    import matplotlib.pyplot as plt
    from pylab import savefig
    import io
    import base64 
    img = io.BytesIO()
    svm=sns.distplot(input_variable1, hist=False, rug=True)
    figure = svm.get_figure()    
    figure.savefig(img, format='png')
    figure.clf()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    return'<img src="data:image/png;base64,{}">'.format(plot_url)
##Bars

@app.route('/bars.html')
def transform5():
    return render_template('bars.html')

@app.route('/bars_result', methods=['POST'])
def bar_results():
#creates a barplot
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    input_variable2 =request.form['input_variable2']
    input_variable2 = input_variable2.splitlines()
    import seaborn as sns
    import io
    import base64
    img = io.BytesIO()

    svm=sns.barplot(x=input_variable2,y=input_variable1)
    figure = svm.get_figure()
    figure.savefig(img, format='png')
    figure.clf()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    return'<img src="data:image/png;base64,{}">'.format(plot_url)
##Pie
@app.route('/pie.html')
def pie():
    return render_template('pie.html')

@app.route('/pie_result', methods=['POST'])
def pie_results():
#creates a pie plot
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    input_variable2 =request.form['input_variable2']
    input_variable2 = input_variable2.splitlines()
    
    import matplotlib.pyplot as plt
    import io
    import base64 
    img = io.BytesIO()
    
    plt.pie(input_variable1, labels=input_variable2)
    plt.savefig(img, format='png')
    plt.close() 
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode()
    return'<img src="data:image/png;base64,{}">'.format(plot_url)

#################################################################
#
# statistical tests one variable
#
#################################################################

@app.route('/statistical_test1.html')
def test1():
    return render_template('statistical_test1.html')

@app.route('/statistical_test1_result', methods=['POST'])
def test1_results():
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]
    test1 = request.form['statistical_test1']
    if test1=="normality":# shapiro test
        output = statisticstools.normality_evaluation(input_variable1) 
        print(output) 
        result=render_template('statistical_test1_result.html', **locals())
    elif test1=="mean": # t test
        constant =float(request.form['input_variable2'])
        output = statisticstools.ttest_one(input_variable1,constant) 
        print(output)
        result=render_template('statistical_test1_result2.html', **locals())
    elif test1=="nonmean":# wilcoxon test
        constant =float(request.form['input_variable2'])
        output = statisticstools.wilcoxontest_one(input_variable1,constant) 
        print(output)
        result=render_template('statistical_test1_result3.html', **locals())
    return result
#################################################################
#
# statistical tests two variables
#
#################################################################

@app.route('/statistical_test2.html')
def test2():
    return render_template('statistical_test2.html')

@app.route('/statistical_test2_result', methods=['POST'])
def test2_results():
    input_variable1 =request.form['input_variable1']
    input_variable1 = input_variable1.splitlines()
    input_variable1 = [float(i) for i in input_variable1]

    input_variable2 =request.form['input_variable2']
    input_variable2 = input_variable2.splitlines()
    input_variable2= [float(i) for i in input_variable2]
    test2 = request.form['statistical_test2']

    if test2=="2meanequal": #t test with equal means
        output = statisticstools.ttest_equal(input_variable1,input_variable2) 
        print(output) 
        result=render_template('2meanequal.html', **locals())
    elif test2=="2meannonequal": # t test unequal means
        output = statisticstools.ttest_unequal(input_variable1,input_variable2)  
        print(output)
        result=render_template('2meannonequal.html', **locals())
    elif test2=="2nonmean":# wilcoxon test for the mean
        output = statisticstools.wilcoxon(input_variable1,input_variable2)  
        print(output)
        result=render_template('2nonmean.html', **locals())
    elif test2=="meanpaired": #t test paired samples
        output = statisticstools.ttest_related(input_variable1,input_variable2)  
        print(output)
        result=render_template('2meanequal.html', **locals())
    elif test2=="nonmeanpaired": #wilcoxon test paired samples
        output = statisticstools.wilcoxon_related(input_variable1,input_variable2)  
        print(output)
        result=render_template('2nonmean.html', **locals())
    elif test2=="variances": #bartlett test for equality of variances
        output = statisticstools.bartlett(input_variable1,input_variable2)  
        print(output)
        result=render_template('variances.html', **locals())
    elif test2=="nonvariances":#levene test for equality of variances
        output = statisticstools.levene(input_variable1,input_variable2)  
        print(output)
        result=render_template('nonvariances.html', **locals())
    elif test2=="proportions": #test for equality of proportions
        output = statisticstools.prop_test(input_variable1,input_variable2)  
        print(output)
        result=render_template('proportions.html', **locals())
    return result

#################################################################
#
# Video tutorial
#
#################################################################
@app.route('/tutorial.html')
def video():
#go to the video hosted in youtube
    return render_template('tutorial.html')


#################################################################
#
# contact
#
#################################################################
@app.route('/email.html')
def email():
#go to my contact page and allow users to send me suggestions about the program
    return render_template('email.html')

# Start the app (let's keep debug=True during debugging)
app.run(debug=True)
