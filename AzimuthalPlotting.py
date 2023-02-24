import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import pickle

# --- Parameters
pickleFile = 'data/AzimuthalData.pkl'
figpath ='./figs/'
cases =['min30deg', '00deg', '30deg']
tools =['BEM','FVW','AMR']
colors=pl.cm.tab20b(np.linspace(0,1,10))
fontsize = 14
plt.rc('font', family='serif')
plt.rc('font', size=14)



# --- Load data
dfsAzi = pickle.load(open(pickleFile,'rb'))


# --- utils
def no_unit(s):
    s=s.replace('_[',' [')
    iu=s.rfind(' [')
    if iu>1:
        return s[:iu]
    else:
        return s

# from prettyplotlib.utils import remove_chartjunk
def remove_chartjunk(ax, spines, grid=None, ticklabels=None, show_ticks=False,
                     xkcd=False):
    '''
    Removes "chartjunk", such as extra lines of axes and tick marks.
    If grid="y" or "x", will add a white grid at the "y" or "x" axes,
    respectively
    If ticklabels="y" or "x", or ['x', 'y'] will remove ticklabels from that
    axis
    '''
    all_spines = ['top', 'bottom', 'right', 'left', 'polar']
    for spine in spines:
        # The try/except is for polar coordinates, which only have a 'polar'
        # spine and none of the others
        try:
            ax.spines[spine].set_visible(False)
        except KeyError:
            pass

    # For the remaining spines, make their line thinner and a slightly
    # off-black dark grey
    if not xkcd:
        for spine in set(all_spines).difference(set(spines)):
            # The try/except is for polar coordinates, which only have a
            # 'polar' spine and none of the others
            try:
                ax.spines[spine].set_linewidth(0.5)
            except KeyError:
                pass
                # ax.spines[spine].set_color(almost_black)
                # ax.spines[spine].set_tick_params(color=almost_black)
                # Check that the axes are not log-scale. If they are, leave
                # the ticks because otherwise people assume a linear scale.
    x_pos = set(['top', 'bottom'])
    y_pos = set(['left', 'right'])
    xy_pos = [x_pos, y_pos]
    xy_ax_names = ['xaxis', 'yaxis']

    for ax_name, pos in zip(xy_ax_names, xy_pos):
        axis = ax.__dict__[ax_name]
        # axis.set_tick_params(color=almost_black)
        #print 'axis.get_scale()', axis.get_scale()
        if show_ticks or axis.get_scale() == 'log':
            # if this spine is not in the list of spines to remove
            for p in pos.difference(spines):
                #print 'p', p
                axis.set_tick_params(direction='out')
                axis.set_ticks_position(p)
                #                axis.set_tick_params(which='both', p)
        else:
            axis.set_ticks_position('none')

    if grid is not None:
        for g in grid:
            assert g in ('x', 'y')
            ax.grid(axis=grid, color='white', linestyle='-', linewidth=0.5)

    if ticklabels is not None:
        if type(ticklabels) is str:
            assert ticklabels in set(('x', 'y'))
            if ticklabels == 'x':
                ax.set_xticklabels([])
            if ticklabels == 'y':
                ax.set_yticklabels([])
        else:
            assert set(ticklabels) | set(('x', 'y')) > 0
            if 'x' in ticklabels:
                ax.set_xticklabels([])
            elif 'y' in ticklabels:
                ax.set_yticklabels([])



# --- Plot
labels=['-30$^\circ$','0$^\circ$','30$^\circ$']

lab_alm_ff = [('Rotor\nTorque [kN-m]', 'RotTorq_[kN-m]'),
             ('OoP Blade-Root\nBending Moment [kN-m]]', 'RootMOoP1_[kN-m]'),
             ('Fore/Aft Tower-Base\nBending Moment [kN-m]', 'TwrBsMyt_[kN-m]'),
             ('Tower-Base\nYaw Moment [kN-m]', 'TwrBsMzt_[kN-m]')
             ]

for k,var in enumerate(lab_alm_ff):
    fig = plt.figure(figsize=(8,4))
    axbig = fig.add_subplot(111)
#     #axbig.set_ylabel(label)
    axbig.set_xlabel("Azimuthal Angle [deg]",labelpad=5)
    axbig.spines['top'].set_color('none')
    axbig.spines['bottom'].set_color('none')
    axbig.spines['left'].set_color('none')
    axbig.spines['right'].set_color('none')
    axbig.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    axbig.tick_params(axis='both', which='both',length=0)
    ax  = {}
    ax[0] = fig.add_subplot(1,2,1)
    ax[1] = fig.add_subplot(1,2,2)
    ax[0].tick_params(axis='both', which='both',length=0)
    ax[1].tick_params(axis='both', which='both',length=0)
    for j,yaw in enumerate(['min30deg','00deg','30deg']):
        label, ffname = var
        psi = dfsAzi['BEM'][yaw].index.values
        BEM = dfsAzi['BEM'][yaw][ffname]
        FVW = dfsAzi['FVW'][yaw][ffname]
        AMR = dfsAzi['AMR'][yaw][ffname]
        BEM_PDiff = (BEM-AMR)/((BEM+AMR)/2.)*100.
        FVW_PDiff = (FVW-AMR)/((FVW+AMR)/2.)*100.

        ax[0].plot(psi,BEM, color = colors[2*j+1],ls = '--')
        ax[0].plot(psi,FVW, color = colors[2*j+1],ls = '-',label = labels[j])
        ax[0].plot(psi,AMR, color = colors[2*j+1],ls = '-.')
        
        ax[1].plot(psi, BEM_PDiff, color=colors[2*j+1],ls='--')
        ax[1].plot(psi, FVW_PDiff, color=colors[2*j+1],ls='-',label=labels[j])

        ax[1].set_ylabel('Percent Difference [%]',labelpad=10,fontsize=fontsize)

        if yaw == '30deg':
            ax[1].plot(psi ,BEM_PDiff*np.nan,color=colors[2*j+1], ls='--',label='BEM')
            ax[1].plot(psi ,FVW_PDiff*np.nan,color=colors[2*j+1], ls='-',label='OLAF')
            ax[1].plot(psi ,FVW_PDiff*np.nan,color=colors[2*j+1], ls='-.',label='SOWFA')

    ax[1].legend(ncol=2, fontsize=fontsize-2, handlelength=2, frameon=True,loc='upper center')
    ax[0].set_ylabel(label)
    ax[1].set_ylabel('Percent Difference [%]',labelpad=10,fontsize=fontsize)
    ax[0].grid()
    ax[1].grid()
    ax[0].spines['top'].set_color('none')
    ax[0].spines['right'].set_color('none')
    ax[1].spines['top'].set_color('none')
    ax[1].spines['right'].set_color('none')
    remove_chartjunk(ax[0], ['top', 'right'])
    remove_chartjunk(ax[1], ['top', 'right'])
    fname = os.path.join(figpath,'Bar02_{0}_Azim_Yaw2.pdf'.format(no_unit(ffname)))
    print(fname)
    fig.tight_layout()
    fig.savefig(fname, bbox_to_inches='tight',dpi=500)



# --- Individual plots
labels=['-30$^\circ$','0$^\circ$','30$^\circ$']

label, ffname = ('OoP Blade-Root\nBending Moment [kN-m]]', 'RootMOoP1_[kN-m]')


for j,yaw in enumerate(['min30deg','00deg','30deg']):
    fig = plt.figure(figsize=(8,4))
    axbig = fig.add_subplot(111)
    axbig.set_xlabel("Azimuthal Angle [deg]",labelpad=5)
    axbig.spines['top'].set_color('none')
    axbig.spines['bottom'].set_color('none')
    axbig.spines['left'].set_color('none')
    axbig.spines['right'].set_color('none')
    axbig.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    axbig.tick_params(axis='both', which='both',length=0)
    axbig.set_title('yaw = ' + labels[j])
    ax  = {}
    ax[0] = fig.add_subplot(1,2,1)
    ax[1] = fig.add_subplot(1,2,2)
    ax[0].tick_params(axis='both', which='both',length=0)
    ax[1].tick_params(axis='both', which='both',length=0)
    psi = dfsAzi['BEM'][yaw].index.values
    BEM = dfsAzi['BEM'][yaw][ffname]
    FVW = dfsAzi['FVW'][yaw][ffname]
    AMR = dfsAzi['AMR'][yaw][ffname]
    BEM_PDiff = (BEM-AMR)/((BEM+AMR)/2.)*100.
    FVW_PDiff = (FVW-AMR)/((FVW+AMR)/2.)*100.
    ax[0].plot(psi,BEM, color = colors[2*j+1],ls = '--')
    ax[0].plot(psi,FVW, color = colors[2*j+1],ls = '-' )
    ax[0].plot(psi,AMR, color = colors[2*j+1],ls = '-.')
    ax[1].plot(psi, BEM_PDiff, color=colors[2*j+1],ls='--')
    ax[1].plot(psi, FVW_PDiff, color=colors[2*j+1],ls='-')
    ax[1].set_ylabel('Percent Difference [%]',labelpad=10,fontsize=fontsize)

    ax[1].plot(psi ,BEM_PDiff*np.nan,color=colors[2*j+1], ls='--',label='BEM')
    ax[1].plot(psi ,FVW_PDiff*np.nan,color=colors[2*j+1], ls='-',label='OLAF')
    ax[1].plot(psi ,FVW_PDiff*np.nan,color=colors[2*j+1], ls='-.',label='SOWFA')

    ax[1].legend(ncol=1, fontsize=fontsize-2, handlelength=2, frameon=True,loc='upper center')
    ax[0].set_ylabel(label)
    ax[1].set_ylabel('Percent Difference [%]',labelpad=10,fontsize=fontsize)
    ax[0].grid()
    ax[1].grid()
    ax[0].spines['top'].set_color('none')
    ax[0].spines['right'].set_color('none')
    ax[1].spines['top'].set_color('none')
    ax[1].spines['right'].set_color('none')
    remove_chartjunk(ax[0], ['top', 'right'])
    remove_chartjunk(ax[1], ['top', 'right'])
    fname = os.path.join(figpath,'{}_Azim_Yaw{}.pdf'.format(no_unit(ffname),yaw))       
    print(fname)
    fig.tight_layout()
    fig.savefig(fname, bbox_to_inches='tight',dpi=500)



plt.show()
