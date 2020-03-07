#!/usr/bin/env python3

import numpy as np
import pandas as pd
import seaborn as sn
import argparse as ap
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

font = {'size'   : 20}
plt.rc('font', **font)

parser = ap.ArgumentParser(description='further assess training and test split, specifically quiescent proportions')
parser.add_argument('--trainSet', action = 'store', nargs=1, type=str, required = True, help='npz file with train set state proportions')
parser.add_argument('--testSet', action = 'store', nargs=1, type=str, required=True, help='npz file with test set state proportions')
args = parser.parse_args()
trainSet = args.trainSet[0]
testSet = args.testSet[0]

train_npz = np.load(trainSet, allow_pickle=True)
train_props = train_npz['props'].astype(np.float64)
train_cell_index = train_npz['cellIndex']
test_npz = np.load(testSet, allow_pickle=True)
test_props = test_npz['props'].astype(np.float64)
#test_cell_index = test_npz['cellIndex'] they match; I checked
cellN = train_cell_index.shape[0]

index = pd.MultiIndex.from_product([['train', 'test'], train_cell_index], names=['set', 'celltype'])
df = pd.DataFrame(columns = index)

for i, ct in enumerate(train_cell_index):
    df.loc[:,('train',ct)] = train_props[:,i,0]
    df.loc[:,('test',ct)] = np.hstack((test_props[:,i,0], np.array([np.nan]*(train_props.shape[0]-test_props.shape[0])).reshape((-1,))))

colors = ['dodgerblue', 'darkorange', 'saddlebrown', 'coral', 'darkslategray', 'purple', 'olive', 'rosybrown','black', 'khaki','grey', 'tan', 'darkgoldenrod', 'turquoise', 'mediumvioletred', 'cyan', 'navy', 'indigo', 'magenta', 'deeppink']


'''Violin plots'''
locations = []
datasets_to_plot = []
labels = []
label_locs = []
next_start = 0
for i, ct in enumerate(train_cell_index):
    locations.append(next_start)
    locations.append(next_start + 1.5)

    datasets_to_plot.append(df.loc[:,('train', ct)])
    datasets_to_plot.append(df.loc[:,('test', ct)][~np.isnan(df.loc[:,('test', ct)])])

    label_locs.append((next_start + next_start + 1.5)/2)
    labels.append(ct)

    next_start = locations[-1] + 3

fig, ax = plt.subplots(figsize=(28,7))
ax.set_title('Proportion of Region in Quiescent State')
parts = ax.violinplot(datasets_to_plot, positions=locations, showmeans = True, showextrema= True)
#ax.boxplot(datasets_to_plot, positions=locations)
ax.set_ylabel('proportion quiescent')
ax.set_xlabel('cell type')
ax.set_xticks(label_locs)
ax.set_xticklabels(labels)
for i, pc in enumerate(parts['bodies']):
    if i%2 == 0:
        pc.set_facecolor('C0')
        pc.set_edgecolor('C0')
    elif i%2 == 1:
        pc.set_facecolor('orange')
        pc.set_edgecolor('orange')
    pc.set_alpha(1)
for partname in ['cmeans', 'cmins', 'cbars', 'cmaxes']:
    parts[partname].set_edgecolor('black')
    parts[partname].set_linewidth(1)

train_patch = mpatches.Patch(color='C0', label = 'Train Region')
test_patch = mpatches.Patch(color='orange', label='Test Region')
plt.legend(handles =[train_patch, test_patch], loc='upper left')
fig.savefig('prop_in_quiescent.png')
plt.close(fig)

'''Ridge plots'''
fig, (ax1, ax2) = plt.subplots(2,1)

ridge_df = pd.DataFrame(index = index, columns=np.arange(20))
ridge_df2 = pd.DataFrame(index = index, columns=np.arange(20))

for i, ct in enumerate(train_cell_index):
    y_tr, bins_tr, patches = ax1.hist(train_props[:,i,0], bins=20)
    y_tr /= np.sum(y_tr)
    ridge_df.loc[('train', ct), :] = y_tr
    bincenters_tr = 0.5*(bins_tr[1:]+bins_tr[:-1])
    ridge_df2.loc[('train', ct), :] = bincenters_tr
    y_t, bins_t, patches = ax2.hist(test_props[:,i,0], bins=20)
    y_t /= np.sum(y_t)
    ridge_df.loc[('test', ct), :] = y_t
    bincenters_t = 0.5*(bins_t[1:]+bins_t[:-1])
    ridge_df2.loc[('test', ct), :] = bincenters_t

ridge_df = pd.concat([ridge_df[col] for col in ridge_df])
ridge_df = ridge_df.reset_index()
ridge_df = ridge_df.rename(columns={0: "y"})

ridge_df2 = pd.concat([ridge_df2[col] for col in ridge_df2])
ridge_df2 = ridge_df2.reset_index()
ridge_df2 = ridge_df2.rename(columns={0: "proportion of region in quiescent state" })

ridge_df["proportion of region in quiescent state"] = ridge_df2["proportion of region in quiescent state"]

pal = sn.cubehelix_palette(10, rot=-.25, light=.7)
g = sn.FacetGrid(ridge_df, col="set", row="celltype", hue="celltype", aspect=15, height=.5, palette = pal, sharex=True)

g.map(plt.plot, "proportion of region in quiescent state", "y")
g.map(plt.axhline, y=0, lw=2, clip_on=False)

def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight='bold', color=color, ha='left', va='center', transform=ax.transAxes)

g.map(label, "proportion of region in quiescent state")
g.fig.subplots_adjust(hspace=0.4, wspace=0.2)
g.set_titles("")
g.set(yticks=[])
g.despine(bottom=True, left=True)
axes = g.axes.flatten()
axes[0].set_title("Train Regions")
axes[1].set_title("Test Regions")
# plt.xlim(0,1)
sn.set(font_scale=0.6)
plt.tight_layout()
g.fig.savefig('ridge_line_prop_in_quiescent.png')
plt.close(g.fig)


'''overlapping histograms and lines. Yuck looking'''
fig1, (ax1, ax2) = plt.subplots(1,2)
fig2, (ax3, ax4) = plt.subplots(1,2)
fig1.suptitle('proportion of region in quiescent state')
ax1.set_ylabel('number of regions')
ax1.set_xlabel('proportion quiescent')
ax1.set_xlim(0,1)
ax1.set_title('Train Set')
ax2.set_xlabel('proportion of region in quiescent state')
ax2.set_title('Test Set')
ax2.set_xlim(0,1)
fig2.suptitle('proportion of region in quiescent state')
ax3.set_ylabel('number of regions')
ax3.set_xlabel('proportion quiescent')
ax3.set_title('Train Set')
ax3.set_xlim(0,1)
ax4.set_xlabel('proportion of region in quiescent state')
ax4.set_title('Test Set')
ax4.set_xlim(0,1)
for i, (dataset_train, dataset_test) in enumerate(zip(train_datasets, test_datasets)):
    y_tr, bins_tr, patches = ax1.hist(dataset_train, bins=20, color=colors[i%len(colors)])
    bincenters_tr = 0.5*(bins_tr[1:]+bins_tr[:-1])
    y_t, bins_t, patches = ax2.hist(dataset_test, bins=20, color=colors[i%len(colors)])
    bincenters_t = 0.5*(bins_t[1:]+bins_t[:-1])
    ax3.plot(bincenters_tr, y_tr, color=colors[i%len(colors)])
    ax4.plot(bincenters_t, y_t, color=colors[i%len(colors)])
#plt.legend(celltypes)
fig1.savefig("quiescent_split_hist.png")
fig2.savefig('quiescent_split_line.png')
plt.close(fig1)
plt.close(fig2)
