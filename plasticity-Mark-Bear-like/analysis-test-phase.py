# %% [markdown]
# # Plasticity of Grating Representations

# %%
# general python modules for scientific analysis
import sys, pathlib, os
import numpy as np

sys.path += ['../physion/src'] # add src code directory for physion
from physion.utils import plot_tools as pt

from physion.analysis.read_NWB import Data,\
    scan_folder_for_NWBfiles

from physion.analysis.episodes.build import EpisodeData

# %%
folder = os.path.join(os.path.expanduser('~'), 
                'DATA', 'Taddy', 'PN_shGrid1-2026', 'NWBs')

notebook_folder =\
    os.path.join(os.path.expanduser('~'), 
        'Documents', 'Notebook', 'Projects', 'Taddy-GluN3', 'figs')

if not os.path.isdir(os.path.join(folder, 'temp')):
    os.mkdir(os.path.join(folder, 'temp'))

# %%
learningRESPS = np.load('learning-resps.npy', allow_pickle=True).item()



# %%
dataset = scan_folder_for_NWBfiles(\
        os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026'),
            for_protocol='Testing-Novel-Familiar-Grating',
            )

quantity = 'dFoF'

# %%

dFoF_parameters = dict(\
    roi_to_neuropil_fluo_inclusion_factor=0., # no factor here
    neuropil_correction_factor = 0.5,
    # with_computed_neuropil_fact=True, # no factor here
    method_for_F0 = 'sliding_percentile',
    percentile=10., # percent
    sliding_window = 5*60, # seconds
)

quantities = ['dFoF']

def analyze_subject(filename,
                    subject,
                    resp_window=[-2,5],
                    shift=.5):

    fig, ax = pt.figure(wspace=0.3, left=2., right=2.,
                        ax_scale=(1.8,1.5))

    data = Data(filename)
    data.build_dFoF(**dFoF_parameters, verbose=True)

    resp = {}
    for i, key, name, color in zip(range(2),\
        ['familiar', 'novel'], 
        ['protocol-1', 'protocol-2'], ['k', 'C2']):
        ep = EpisodeData(data, 
                        quantities=quantities,
                        protocol_name=name)
        baseline = ep.dFoF[:,:,ep.t<0].mean()
        ax.plot(ep.t, -baseline+ep.dFoF.mean(axis=(0,1)), color=color)
        pt.annotate(ax, 'n=%i ROIs ' % data.nROIs, (0, 0.3), ha='right')
        pt.annotate(ax, os.path.basename(filename).replace('.nwb', ''),
                    (0,0), ha='right', fontsize=5)
        pt.annotate(ax, key+i*'\n', (1, 0.4), color=color)

        resp[key] = ep.dFoF.mean(axis=1) # average over ROIs

    pt.draw_bar_scales(ax, Ybar=0.1, Ybar_label='0.1$\\Delta$F/F ', 
                       Xbar=1e-3)
    ylim = ax.get_ylim()
    ax.fill_between([0, ep.time_duration[0]], ylim[0], ylim[1], color='k', alpha=.1, lw=0)
    ax.set_ylim(ylim)
    pt.set_plot(ax, ['bottom'])
    ax.set_xlabel(' time (s) ')
    ax.set_title(subject)

    return ep.t, resp, fig

RESPS = {'Grid1':{'familiar':[], 'novel':[]},
         'Scramble':{'familiar':[], 'novel':[]}}
for s, subject in enumerate(\
                np.unique(dataset['subjects'])):

    subject_files = dataset['files'][subject==dataset['subjects']]

    virus=\
         dataset['viruses'][subject==dataset['subjects']][0].replace('CamKII-Cre+sh', '')

    t, resp, fig = analyze_subject(subject_files[0],
                    '%s -- %s ' % (subject, virus))
    fig_name = 'markBear-plasticity-Test-trial-average-%s.svg' % subject
    pt.save(fig, notebook_folder, fig_name)
    for key in ['familiar', 'novel']:
        RESPS[virus][key].append(resp[key])

for virus in RESPS:
    for key in ['familiar', 'novel']:
        RESPS[virus][key] = np.array(RESPS[virus][key])
RESPS['t'] = t

# %%
# -- print for markdown notebook
for s, subject in enumerate(\
                np.unique(dataset['subjects'])):
    fig_name = 'markBear-plasticity-Test-trial-average-%s.svg' % subject
    print('![](figs/%s)     ' % fig_name) # for notebook !

# %%
from scipy import stats
fig, ax = pt.figure(left=1.1, ax_scale=(2,2), right=2.)
# fig2, ax2 = pt.figure(left=1.1, ax_scale=(1.5,1.5))

pre_window = [-2, 0]
post_window= [0, 4]

for v, virus in enumerate(['Grid1', 'Scramble']):
    points = np.zeros((RESPS[virus]['familiar'].shape[0],2))
    for i, key, color in zip(range(2), ['familiar', 'novel'], ['grey', 'C2']):
        pre_cond = (RESPS['t']>pre_window[0]) & (RESPS['t']<pre_window[1])
        baseline = np.mean(RESPS[virus][key][:,:,pre_cond], axis=(1,2))
        post_cond = (RESPS['t']>post_window[0]) & (RESPS['t']<post_window[1])
        resp = ((RESPS[virus][key][:,:,post_cond].mean(axis=(1,2))).T-baseline.T).T
        # average responses
        ax.bar([i+3*v], [resp.mean()], yerr=[stats.sem(resp)], color = color)
        points[:,i] = resp
    # individual responses
    ax.plot(np.arange(2)*0.8+0.1+3*v, points.T, 'ko-', ms=1, lw=0.4)
    # now relative plot
    # ax2.bar([v], [np.mean((points[:,1]-points[:,0])/points[:,1])*100.], 
    #        yerr=[stats.sem(resp)], color = 'tab:blue')

for i, key, color in zip(range(2), ['familiar', 'novel'], ['grey', 'C2']):
    pt.annotate(ax, key+i*'\n', (1,0), color=color)
pt.set_plot(ax, 
            xticks=[0.5, 3.5],
            xticks_labels=\
            ['%s\n(N=%i)' % (key, RESPS[key]['familiar'].shape[0]) for key in ['Grid1', 'Scramble']],
            ylabel='$\\delta$ $\\Delta$F/F')
pt.save(fig, notebook_folder, 
        'markBear-plasticity-Test-summary.svg',
        transparent=True)

# %%

for v, virus in enumerate(['Grid1', 'Scramble']):

    fig, AX = pt.figure(axes=(10,RESPS[virus]['familiar'].shape[0]),
                        ax_scale=(0.8,0.9),
                        wspace=0.2, hspace=0.4, left=2., right=2.)
    fig.suptitle(virus, fontsize=10)
    for n in range(RESPS[virus]['familiar'].shape[0]):

        for i, key, color in zip(range(2),\
                    ['familiar', 'novel'], ['k', 'C2']):

            for k in range(10):
                AX[n][k].plot(RESPS['t'], 
                              RESPS[virus][key][n, 10*k:10*k+10, :].mean(axis=0),
                              color=color)
                AX[n][k].axis('off')

            pt.annotate(AX[n][0], 'mouse %i' % (n+1), (0,0.5),
                        ha='right', va='top')
        pt.draw_bar_scales(AX[n][0], Ybar=0.4, Ybar_label='0.4$\\Delta$F/F', Xbar=1, Xbar_label='1s')
        pt.set_common_ylims(AX[n])
        for k in range(10):
            pt.annotate(AX[0][k], '%i-%i' % (10*k,10*(k+1)-1), (1,1),
                        fontsize=5, ha='right', va='top')
    fig_name = 'markBear-plasticity-Test-trial-dynamics-%s.svg' % virus
    pt.save(fig, notebook_folder, fig_name,
            transparent=True)
    print('![](figs/%s)     ' % fig_name) # for notebook !

# %%
