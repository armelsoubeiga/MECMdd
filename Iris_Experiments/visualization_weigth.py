import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from pypalettes import load_cmap

df_l = pd.read_csv('Iris/df_l_weight.csv', index_col=False)
df_g = pd.read_csv('Iris/df_g_weight.csv', index_col=False)

num_gradient_steps = 20
# colors = ['#1228a4', '#8e92cd', '#e6dac5', '#d4a10b']
# cmap = LinearSegmentedColormap.from_list("custom_palette", colors, N=256)
cmap = load_cmap("BluetoOrangeRed_14")

# ls viz
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((1.5, 0.5, 1.5))
ls = df_l['ls']
xpos = np.array([1, 2, 3, 4, 5, 6, 7] * 2)
ypos = np.array([1] * 7 + [1.2] * 7) 
zpos = np.zeros_like(ls)
dx = 0.7
dy = 0.1
dz = ls
co = ['Petal','Sepal']
ro = ['$\omega_1$','$\omega_2$','$\omega_1,\omega_2$','$\omega_3$',
             '$\omega_1,\omega_3$','$\omega_2,\omega_3$','$\Omega$'] 
max_ls = ls.max()
for idx in range(len(xpos)):
    for j in range(num_gradient_steps):
        z_current = zpos[idx] + j * dz[idx] / num_gradient_steps
        dz_current = dz[idx] / num_gradient_steps

        color_value = (z_current + dz_current) / max_ls
        color = cmap(color_value)
        ax.bar3d(xpos[idx]-0.25, ypos[idx]-0.08, z_current, dx, dy, dz_current, color=color)

ax.set_xlabel('Clusters', fontsize=12)
ax.set_ylabel('View', fontsize=12, rotation=45)
ax.set_zlabel('View Weight', labelpad=0.5, fontsize=12, rotation=90)
ax.set_xticks([1, 2, 3, 4, 5, 6, 7])
ax.set_yticks([1, 1.2])
ax.set_xticklabels(ro, fontsize=9, verticalalignment='baseline', horizontalalignment='right')
ax.set_yticklabels(co, fontsize=9, verticalalignment='baseline', horizontalalignment='left')
ax.tick_params(axis='x', labelsize=9)
ax.tick_params(axis='y', labelsize=9)
ax.tick_params(axis='z', labelsize=9)

# ax.zaxis.set_rotate_label(False)
#ax.view_init(elev=25, azim=-30)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=ls.min(), vmax=ls.max()))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.02, fraction=0.022, format='%.1f', shrink=0.8)
cbar.ax.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig("ls.png", format='png', dpi=300)
plt.show()



# lp viz
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((1.5, 0.5, 1.5))
lp = df_l['lp']
xpos = np.array([1, 2, 3, 4, 5, 6, 7] * 2)
ypos = np.array([1] * 7 + [1.2] * 7) 
zpos = np.zeros_like(lp)
dx = 0.7
dy = 0.1
dz = lp
co = ['Petal','Sepal']
ro = ['$\omega_1$','$\omega_2$','$\omega_1,\omega_2$','$\omega_3$',
             '$\omega_1,\omega_3$','$\omega_2,\omega_3$','$\Omega$'] 
max_lp = lp.max()
for idx in range(len(xpos)):
    for j in range(num_gradient_steps):
        z_current = zpos[idx] + j * dz[idx] / num_gradient_steps
        dz_current = dz[idx] / num_gradient_steps

        color_value = (z_current + dz_current) / max_lp
        color = cmap(color_value)
        ax.bar3d(xpos[idx]-0.25, ypos[idx]-0.08, z_current, dx, dy, dz_current, color=color)

ax.set_xlabel('Clusters', fontsize=12)
ax.set_ylabel('View', fontsize=12, rotation=45)
ax.set_zlabel('View Weight', labelpad=0.5, fontsize=12, rotation=90)
ax.set_xticks([1, 2, 3, 4, 5, 6, 7])
ax.set_yticks([1, 1.2])
ax.set_xticklabels(ro, fontsize=9, verticalalignment='baseline', horizontalalignment='right')
ax.set_yticklabels(co, fontsize=9, verticalalignment='baseline', horizontalalignment='left')
ax.tick_params(axis='x', labelsize=9)
ax.tick_params(axis='y', labelsize=9)
ax.tick_params(axis='z', labelsize=9)

# ax.zaxis.set_rotate_label(False)
#ax.view_init(elev=25, azim=-30)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=lp.min(), vmax=lp.max()))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.02, fraction=0.022, format='%.1f', shrink=0.8)
cbar.ax.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig("lp.png", format='png', dpi=300)
plt.show()





#gs viz
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((0.5, 0.8, 1))
gs = df_g['gs']
xpos = [1.5,1.5]
ypos = [1, 1.2] 
zpos = np.zeros_like(gs)
dx = 0.03; dy = 0.1
dz = gs
co = ['Petal','Sepal']
max_gs = gs.max()

for idx in range(len(xpos)):
    for j in range(num_gradient_steps):
        z_current = zpos[idx] + j * dz[idx] / num_gradient_steps
        dz_current = dz[idx] / num_gradient_steps

        color_value = (z_current + dz_current) / max_gs
        color = cmap(color_value)
        ax.bar3d(xpos[idx]+0.03, ypos[idx]-0.08, z_current, dx, dy, dz_current, color=color)

ax.set_ylabel('View', fontsize=12, rotation=45, labelpad=5)
ax.set_zlabel('View Weight', labelpad=10, fontsize=12, rotation=90)
ax.set_xticks([1.5,1.5])
ax.set_yticks([1, 1.2])
ax.set_xticklabels('')
ax.set_yticklabels(co, fontsize=9, verticalalignment='baseline', horizontalalignment='left')

ax.tick_params(axis='y', labelsize=9)
ax.tick_params(axis='z', labelsize=9)

# ax.zaxis.set_rotate_label(False)
ax.view_init(elev=25, azim=-50)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=gs.min(), vmax=gs.max()))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.02, fraction=0.022, format='%.1f')
cbar.ax.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig("gs.png", format='png', dpi=300)
plt.show()








#gp viz
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((0.5, 0.8, 1))
gp = df_g['gp']
xpos = [1.5,1.5]
ypos = [1, 1.2] 
zpos = np.zeros_like(gp)
dx = 0.03; dy = 0.1
dz = gp
co = ['Petal','Sepal']
max_gp = gp.max()

for idx in range(len(xpos)):
    for j in range(num_gradient_steps):
        z_current = zpos[idx] + j * dz[idx] / num_gradient_steps
        dz_current = dz[idx] / num_gradient_steps

        color_value = (z_current + dz_current) / max_gp
        color = cmap(color_value)
        ax.bar3d(xpos[idx]+0.03, ypos[idx]-0.08, z_current, dx, dy, dz_current, color=color)

ax.set_ylabel('View', fontsize=12, rotation=45, labelpad=5)
ax.set_zlabel('View Weight', labelpad=10, fontsize=12, rotation=90)
ax.set_xticks([1.5,1.5])
ax.set_yticks([1, 1.2])
ax.set_xticklabels('')
ax.set_yticklabels(co, fontsize=9, verticalalignment='baseline', horizontalalignment='left')

ax.tick_params(axis='y', labelsize=9)
ax.tick_params(axis='z', labelsize=9)

# ax.zaxis.set_rotate_label(False)
ax.view_init(elev=25, azim=-50)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=gp.min(), vmax=gp.max()))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.02, fraction=0.022, format='%.1f')
cbar.ax.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig("gp.png", format='png', dpi=300)
plt.show()