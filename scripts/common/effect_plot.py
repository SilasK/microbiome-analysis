import numpy as np
import matplotlib
import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns

def plotting_Setup(sns_contex='paper',font_scale=1.5, save_dpi=150):
    import matplotlib
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_context(sns_contex,font_scale=font_scale)
    matplotlib.rcParams['savefig.dpi']=save_dpi
    matplotlib.rcParams['pdf.fonttype']=42


def saveplot(name,figurefolder='../Figures/',formats=['.pdf','.png'],SAVEPLOT=True):
    matplotlib.rcParams['pdf.fonttype']=42 # to save as vector format
    if SAVEPLOT:
        if not os.path.exists(figurefolder):
            os.makedirs(figurefolder)
        if len(formats)==0:
            raise Exception('no format specified')

        for format in formats:
            plt.savefig(os.path.join(figurefolder,name+format),bbox_inches='tight')


def format_p_value(p,Stars=False):
    if not Stars:
        return 'P = {:.1g}'.format(p).replace('0.','.')
    else:
        if p<0.001:
            return '***'
        elif p<0.01:
            return '**'
        if p < 0.05:
            return '*'
        else:
            return 'n.s.'

def effect_plot(Values,P_values,colors=None,width=0.7,Labels=None, Title= None):

    sns.set_palette(colors)

    from matplotlib.font_manager import FontProperties
    font=FontProperties(size='x-small')


    #assert Values.isnull().any().any()==False

    if colors is None:
        colors= sns.color_palette('deep',n_colors= Values.shape[1])


    ##
    d1,d2= Values.shape

    assert P_values.shape[1] == d2
    Values= Values.loc[reversed(Values.index)]
    P_values= P_values.loc[Values.index]


    ax= Values.plot.barh(width=width,figsize=(10,d1*d2/width/10+1))

#

    Horizontal_alignement={1:'left',-1:'right'}

    for i in range(d1):
        for j in range(d2):

            y= Values.iloc[i,j]
            text= format_p_value(P_values.iloc[i,j],Stars=True)
            if not text =='n.s.':
                ax.annotate(text, xy=(y,i+width/d2*(0.5+j)- width/2),textcoords='offset points',xytext=(np.sign(y),-1),
                    horizontalalignment= Horizontal_alignement[np.sign(y)], verticalalignment='center',font_properties=font)

    if Labels is not None:
        if type(Labels) is pd.Series:
            ax.set_yticklabels(Labels.loc[Values.index])
        else:
            ax.set_yticklabels(Labels)
    if Title is not None:
        ax.set_title(Title)

    handles,labels = ax.get_legend_handles_labels()

    ax.legend(reversed(handles),reversed(labels))

    return ax



def aldex_plot(D):

    f,axe = plt.subplots(1,2,figsize=(15,5))

    sns.set_palette(['orange','darkgrey'])

    sns.scatterplot(data=D,x='rab.all',y='effect',hue='label',ax=axe[0],hue_order=['Significant','Not significant'])

    x= np.array(axe[0].get_xlim())
    axe[0].plot(x,[1,1],'k--')
    axe[0].plot(x ,[-1,-1],'k--')
    axe[0].set_xlabel("Median $\log_2$ relative abundance")
    axe[0].set_ylabel("Effect size [$\log_2$]")
    axe[0].get_legend().remove()
    sns.scatterplot(data=D,x='diff.win',y='diff.btw',hue='label',ax=axe[1],hue_order=['Significant','Not significant'])

    x= np.array(axe[1].get_xlim())
    axe[1].plot(x,x,'k--')
    axe[1].plot(x ,-x,'k--')
    axe[1].set_xlabel("Median $\log_2$ win condition diff")
    axe[1].set_ylabel("Median $\log_2$ btw condition diff")
    axe[1].legend(*axe[1].get_legend_handles_labels(),bbox_to_anchor=(1, 1))

#plt.savefig(name)
