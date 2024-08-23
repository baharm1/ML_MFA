import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score
import scipy

def plt_learning_curve(history, filename = None):
    plt.title('Learning Curve')
    plt.xlabel('Epoch')
    plt.ylabel('Loss (MSE)')
    plt.plot(history.history['loss'], label = 'Train')
    plt.plot(history.history['val_loss'], label = 'Validation')
    plt.legend()
    if filename != None:
        plt.savefig(''.join([filename, '.pdf']), format = 'pdf', bbox_inches = 'tight')
    plt.show()
    
    
def plt_target_pred(target, pred, istrain = False, filename = None):
    fig, ax = plt.subplots(figsize = (5, 5))
    plt.scatter(target, pred, s = 1)
    b, a = np.polyfit(target, pred, deg = 1)

    # Create sequence of 100 numbers from 0 to 100 
    xseq = np.linspace(0, 1, num = 100)
    # Plot regression line
    ax.plot(xseq, a + b * xseq, color = "k", lw = 2.5);
    coefficient_of_dermination = r2_score(target, pred)
    res = scipy.stats.linregress(target, pred)
    if istrain:
        title = 'Train'
    else:
        title = 'Test'
    plt.title(f'R2 = {coefficient_of_dermination :.2f}, y = {res.slope:.3f}x + {res.intercept:.3f}')
    plt.xlabel(' '.join([title, 'Target']))
    plt.ylabel(' '.join([title, 'Prediction']))
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    if filename != None:
        plt.savefig(''.join([filename, '.pdf']), format = 'pdf', bbox_inches = 'tight')
    plt.show()  

    
def plot_regression(x, y, filename = None):
    '''
    Plot the predicted vs ground truth figure. Include R2 and fit line.
    '''
    # Calculate the point density
    plt.rcParams.update({'font.size': 15})
    xy = np.vstack((x,y))
    z = scipy.stats.gaussian_kde(xy)(xy)

    fig, ax = plt.subplots(figsize = (5.5, 5))
    den = ax.scatter(x, y, c = z, s = 1, cmap = 'plasma_r')

    # regression
    xseq = np.linspace(x.min(), x.max(), num = 100)
    res = scipy.stats.linregress(x, y)
    ax.plot(xseq, res.intercept + res.slope * xseq, color = "k", lw = 1.5)
    print(res.pvalue)
    plt.title(f" (Pearson's r = {res.rvalue:.2f})", fontsize = 15)
        
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    ax.set_xlabel('Test Target', fontsize=15)
    ax.set_ylabel('Test Prediction', fontsize=15)

    fig.patch.set_facecolor('white')
    fig.colorbar(den, label = 'Density')
    if filename != None:
        plt.savefig(''.join([filename, '.pdf']), format = 'pdf', bbox_inches = 'tight')
        
    plt.show()
    