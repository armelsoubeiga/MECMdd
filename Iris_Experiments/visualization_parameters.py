from math import gamma
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
from matplotlib import cm
from pypalettes import load_cmap

# Load the data
df = pd.read_csv('Iris/cri_parameters_10.csv')
betas = np.unique(df['beta'])
gammas = np.unique(df['gamma'])
etas = np.unique(df['eta'])

betaslab = [f'{b:.1f}' for b in betas]
gammaslab = [f'{g:.1f}' for g in gammas]
etaslab = [f'{e:.1f}' for e in etas]

num_gradient_steps = 100
cmap = load_cmap("BluetoOrangeRed_14")



def plot_3d_bars(df, algo_name, i):
    fig = plt.figure(figsize=(18, 5)) #(18, 6)
    
    
    # First subplot: beta vs gamma
    ax1 = fig.add_subplot(131, projection='3d')
    df_algo = df[df['algorithme'] == algo_name]

    xpos, ypos = np.meshgrid(betas, gammas) 
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)

    Z = np.array([df_algo[(df_algo['beta'] == b) & (df_algo['gamma'] == g)]['cri'].mean() 
              for b in betas for g in gammas])
    dz = Z
    max_dz = dz.max()

    dx = 0.2 ; dy = 0.22  
    xgap1 = 2.6;  ygap1 = 2.9

    for idx in range(len(xpos)):
        for j in range(num_gradient_steps):
            z_current = zpos[idx] + j * dz[idx] / num_gradient_steps
            dz_current = dz[idx] / num_gradient_steps
            color_value = (z_current + dz_current) / max_dz
            color = cmap(color_value)
            ax1.bar3d(xpos[idx]*xgap1, ypos[idx]*ygap1, z_current, dx, dy, dz_current, color=color)
    
    ax1.set_xlabel(r'$\beta$', fontsize=10)
    ax1.set_ylabel('$\gamma$', rotation=45, fontsize=10)
    ax1.set_zlabel('CRI', labelpad=0.5, rotation=90, fontsize=10)
    
    ax1.set_xticks(betas*xgap1)
    ax1.set_xticklabels(betaslab, fontsize=8, rotation=45, verticalalignment='baseline', horizontalalignment='right')
    ax1.set_yticks(gammas*ygap1)
    ax1.set_yticklabels(gammaslab, fontsize=8, rotation=-45, verticalalignment='baseline', horizontalalignment='left')
    ax1.tick_params(axis='z', labelsize=8)

    ax1.view_init(elev=25, azim=-50)  
    norm = Normalize(vmin=0, vmax=dz.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar1 = plt.colorbar(sm, ax=ax1, pad=0.11, fraction=0.02, format='%.2f', shrink=0.8)
    cbar1.ax.tick_params(labelsize=8)





    # Second subplot: beta vs eta
    ax2 = fig.add_subplot(132, projection='3d')
    xpos, ypos = np.meshgrid(betas, etas)
    xpos = xpos.flatten() 
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)

    Z = np.array([df_algo[(df_algo['beta'] == b) & (df_algo['eta'] == e)]['cri'].mean() 
              for b in betas for e in etas])
    dz = Z
    max_dz = dz.max()

    dx = 0.12 ; dy = 1  
    xgap1 = 1.5;  ygap1 = 1.5

    for idx in range(len(xpos)):
        for j in range(num_gradient_steps):
            z_current = zpos[idx] + j * dz[idx] / num_gradient_steps
            dz_current = dz[idx] / num_gradient_steps
            color_value = (z_current + dz_current) / max_dz
            color = cmap(color_value)
            ax2.bar3d(xpos[idx]*xgap1, ypos[idx]*ygap1, z_current, dx, dy, dz_current, color=color)

    ax2.set_xlabel(r'$\beta$')
    ax2.set_ylabel('$\eta$')
    ax2.set_zlabel('CRI', labelpad=0.5, rotation=90, fontsize=10)
    
    ax2.set_xticks(betas*xgap1)
    ax2.set_xticklabels(betaslab, fontsize=8, rotation=45, verticalalignment='baseline', horizontalalignment='right')
    ax2.set_yticks(etas*ygap1)
    ax2.set_yticklabels(etaslab, fontsize=8, rotation=-45, verticalalignment='baseline', horizontalalignment='left')
    ax2.tick_params(axis='z', labelsize=8)

    ax2.view_init(elev=25, azim=-50)  
    norm = Normalize(vmin=0, vmax=dz.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar2 = plt.colorbar(sm, ax=ax2, pad=0.11, fraction=0.02, format='%.2f', shrink=0.8)
    cbar2.ax.tick_params(labelsize=8)
    plt.title(f'({chr(97+i)}) {algo_name}', fontsize=20, y=-0.15, ha='center')



    # Third subplot: gamma vs eta
    ax3 = fig.add_subplot(133, projection='3d')
    xpos, ypos = np.meshgrid(gammas, etas)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)

    Z = np.array([df_algo[(df_algo['gamma'] == g) & (df_algo['eta'] == e)]['cri'].mean() 
              for g in gammas for e in etas])
    dz = Z
    max_dz = dz.max()

    dx = 0.1 ; dy = 1  
    xgap1 = 1.3;  ygap1 = 1.3

    for idx in range(len(xpos)):
        for j in range(num_gradient_steps):
            z_current = zpos[idx] + j * dz[idx] / num_gradient_steps
            dz_current = dz[idx] / num_gradient_steps
            color_value = (z_current + dz_current) / max_dz
            color = cmap(color_value)
            ax3.bar3d(xpos[idx]*xgap1, ypos[idx]*ygap1, z_current, dx, dy, dz_current, color=color)

    ax3.set_xlabel('$\gamma$')
    ax3.set_ylabel('$\eta$')
    ax3.set_zlabel('CRI', labelpad=0.5, rotation=90, fontsize=10)
    
    ax3.set_xticks(gammas*xgap1)
    ax3.set_xticklabels(gammaslab, fontsize=8, rotation=45, verticalalignment='baseline', horizontalalignment='right')
    ax3.set_yticks(etas*ygap1)
    ax3.set_yticklabels(etaslab, fontsize=8, rotation=-45, verticalalignment='baseline', horizontalalignment='left')
    ax3.tick_params(axis='z', labelsize=8)

    ax3.view_init(elev=25, azim=-50)  
    norm = Normalize(vmin=0, vmax=dz.max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar3 = plt.colorbar(sm, ax=ax3, pad=0.11, fraction=0.02, format='%.2f', shrink=0.8)
    cbar3.ax.tick_params(labelsize=8)




    plt.subplots_adjust(wspace=0.11)
    #plt.suptitle(f'({chr(97+i)}) {algo_name}', fontsize=20, y=-0.1, ha='center')
    plt.tight_layout()
    plt.savefig(f'{algo_name}.png', format='png', dpi=300)
    plt.show()




# Run the function for each algorithm
df['algorithme'] = df['algorithme'].replace({
    'ls': 'MECMdd-RWL-S',
    'lp': 'MECMdd-RWL-P',
    'gs': 'MECMdd-RWG-S',
    'gp': 'MECMdd-RWG-P'
})
algorithms = ['MECMdd-RWL-S', 'MECMdd-RWL-P', 'MECMdd-RWG-S', 'MECMdd-RWG-P']
for i, algo in enumerate(algorithms):
    plot_3d_bars(df, algo, i)



