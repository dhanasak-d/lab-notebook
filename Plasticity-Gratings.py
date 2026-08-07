# %% [markdown]
# # Plasticity of Grating Representations

# %%
# general python modules for scientific analysis
import sys, pathlib, os
import numpy as np

sys.path += ['physion/src'] # add src code directory for physion
from physion.utils import plot_tools as pt

from physion.analysis.read_NWB import Data,\
    scan_folder_for_NWBfiles

from physion.analysis.episodes.build import EpisodeData
from physion.dataviz.episodes.trial_average\
              import plot as plot_trial_average

# %%
folder = os.path.join(os.path.expanduser('~'), 
                'DATA', 'Taddy', 'PN_shGrid1-2026', 'NWBs')

if not os.path.isdir(os.path.join(folder, 'temp')):
    os.mkdir(os.path.join(folder, 'temp'))

dataset = scan_folder_for_NWBfiles(\
        os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026'),
            for_protocols=['Learning-Familiar-Grating-45deg'],
            )

quantity = 'dFoF'

# %%

dFoF_parameters = dict(\
    roi_to_neuropil_fluo_inclusion_factor=0., # no factor here
    # neuropil_correction_factor = 0.5,
    with_computed_neuropil_fact=True, # no factor here
    method_for_F0 = 'sliding_percentile',
    percentile=20., # percent
    sliding_window = 5*60, # seconds
)

quantities = ['dFoF']

def analyze_subject(filenames,
                    subject,
                    resp_window=[-2,8],
                    shift=.5):

    fig, AX = pt.figure(axes=(2,1), 
                        wspace=0.3, left=2.,
                        ax_scale=(1.3,3.))

    resp = []
    for day, filename in enumerate(filenames):

        color = pt.copper(1-day/len(filenames))
        data = Data(filename)
        data.build_dFoF(**dFoF_parameters, verbose=True)
        ep = EpisodeData(data, 
                        quantities=quantities,
                        protocol_id=0)


        baseline = ep.dFoF[:,:,ep.t<0].mean()
        cond = (ep.t>-3) & (ep.t<10)
        AX[0].plot(ep.t[cond], 
            -shift*day-baseline+\
                ep.dFoF.mean(axis=(0,1))[cond],
                color=color)
        cond = (ep.t>90) & (ep.t<103)
        AX[1].plot(ep.t[cond], 
            -shift*day-baseline+\
                ep.dFoF.mean(axis=(0,1))[cond],
            color=color)
        AX[0].annotate('day %i (n=%i ROIs) ' % (day+1, data.nROIs),
                    (-3, -shift*day), ha='right',
                    xycoords='data', color=color)
        AX[0].annotate('\n'+os.path.basename(filename).replace('.nwb', '')+10*' ',
                    (-3, -shift*day), ha='right',
                    xycoords='data', color=color,
                    fontsize=4, va='top')

        resp_cond = (ep.t>resp_window[0]) &\
                        (ep.t<resp_window[1])
        resp.append(\
            ep.dFoF[:,:,(ep.t>resp_window[0]) &\
                (ep.t<resp_window[1])].mean(axis=1))

    pt.set_common_ylims(AX)
    pt.draw_bar_scales(AX[1], Ybar=0.4, Ybar_label='0.4$\Delta$F/F', Xbar=1e-3)
    for ax in AX:
        pt.set_plot(ax, ['bottom'])
    AX[0].set_xlabel(50*' '+' time (s)     --       stim in [0,100]s')
    AX[0].set_title(subject)

    return ep.t[resp_cond], resp, fig

RESPS = {'Grid1':[], 'Scramble':[]}
for s, subject in enumerate(\
                np.unique(dataset['subjects'])):

    subject_files = dataset['files'][subject==dataset['subjects']]
    rec_dates = dataset['dates'][subject==dataset['subjects']]
    virus=\
         dataset['viruses'][subject==dataset['subjects']][0].replace('CamKII-Cre+sh', '')

    t, resp, fig = analyze_subject(subject_files[np.argsort(rec_dates)], 
                    '%s -- %s ' % (subject, virus))
    pt.save(fig, 'Desktop/plasticity', '%s.png' % (s+1),
            transparent=False)
    RESPS[virus].append(resp)

RESPS['t'] = t
for v, virus in enumerate(RESPS):
    RESPS[virus] = np.array(RESPS[virus])

# %%
from scipy import stats
fig, ax = pt.figure(left=1.1)

pre_window = [-2, 0]
post_window= [0, 4]

for v, virus in enumerate(['Grid1', 'Scramble']):
    if len(RESPS[virus])>0:
        # RESPS[virus].shape = (subject, days, repeats, time)
        pre_cond = (RESPS['t']>pre_window[0]) & (RESPS['t']<pre_window[1])
        baseline = RESPS[virus][:,:,:,pre_cond].mean(axis=-1)
        post_cond = (RESPS['t']>post_window[0]) & (RESPS['t']<post_window[1])
        resp = (RESPS[virus][:,:,:,post_cond].T-baseline.T).T
        # average responses
        for i in range(resp.shape[1]):
            ax.bar([i+v*(1+resp.shape[1])], 
                [resp[:,i,:].mean()],
                yerr=[stats.sem(resp[:,i,:].mean(axis=(1,2)))],
                color = pt.copper(1-i/resp.shape[1]))
        # individual responses
        for j in range(resp.shape[0]):
            ax.plot(v*(1+resp.shape[1])+\
                    np.arange(resp.shape[1]),
                    resp[j,:,:,:].mean(axis=(1,2)), 
                    'k-', lw=0.2)
        
pt.set_plot(ax, 
            xticks=[1, 9],
            xticks_labels=\
            ['%s\n(N=%i)' % (key, len(RESPS[key])) for key in ['Grid1', 'Scramble']],
            ylabel='$\delta$ $\Delta$F/F')
pt.save(fig, 'Desktop/plasticity', 'summary.png',
        transparent=False)

# %%
from scipy import stats
from scipy.ndimage import gaussian_filter1d

pre_window = [-2, 0]
post_window= [0, 8]

shift, tshift = 0.6, 12
smoothing = 5
for v, virus in enumerate(['Grid1', 'Scramble']):
    if len(RESPS[virus])>0:
        # RESPS[virus].shape = (subject, days, repeats, time)
        pre_cond = (RESPS['t']>pre_window[0]) & (RESPS['t']<pre_window[1])
        baseline = RESPS[virus][:,:,:,pre_cond].mean(axis=-1)
        post_cond = (RESPS['t']>post_window[0]) & (RESPS['t']<post_window[1])
        resp = (RESPS[virus][:,:,:,:].T-baseline.T).T

        for rec in range(resp.shape[0]):

            fig, ax = pt.figure(ax_scale=(2,3), left=1.1)
            ax.set_title('%s - mouse %i' % (virus, rec+1))
            ax.axis('off')
            for day in range(resp.shape[1]):
                for repeat in range(resp.shape[2]):

                    ax.plot(RESPS['t']+repeat*tshift, 
                        -shift*day+\
                            gaussian_filter1d(resp[rec, day, repeat, :], smoothing),
                        color = pt.copper(1-day/resp.shape[1]))

                    if repeat==0:
                        ax.annotate('day %i ' % (day+1),
                                    (-3, -shift*day), ha='right',
                                    xycoords='data', 
                                    color = pt.copper(1-day/resp.shape[1]))
                    if day==0:
                        ax.annotate('#%i' % (repeat+1),
                                    (tshift*repeat, -shift*(resp.shape[1]-.5)), 
                                    va='top', ha='center', xycoords='data')

            pt.draw_bar_scales(ax, 
                               loc='top-right',
                    Ybar=0.4, Ybar_label='0.4$\Delta$F/F', 
                    Xbar=4, Xbar_label='4s', color='k')
        
# pt.save(fig, 'Desktop/plasticity', '.png',
#         transparent=False)

# %%

pre_window = [-2, 0]
post_window= [0, 8]

shift, tshift = 0.6, 12
smoothing = 50
for v, virus in enumerate(['Grid1', 'Scramble']):
    if len(RESPS[virus])>0:
        # RESPS[virus].shape = (subject, days, repeats, time)
        pre_cond = (RESPS['t']>pre_window[0]) & (RESPS['t']<pre_window[1])
        baseline = RESPS[virus][:,:,:,pre_cond].mean(axis=-1)
        post_cond = (RESPS['t']>post_window[0]) & (RESPS['t']<post_window[1])
        resp = gaussian_filter1d(\
            (RESPS[virus][:,:,:,:].T-baseline.T).T,
            smoothing)

        fig, ax = pt.figure(ax_scale=(2,3), left=1.1)
        ax.set_title('%s (N=%i mice)' % (virus, resp.shape[0]))
        ax.axis('off')
        for day in range(resp.shape[1]):
            for repeat in range(resp.shape[2]):

                pt.plot(RESPS['t']+repeat*tshift, 
                    -shift*day+\
                        resp[:, day, repeat, :].mean(axis=0),
                    sy=stats.sem(resp[:, day, repeat, :], axis=0),
                    color = pt.copper(1-day/resp.shape[1]),
                    ax=ax)

                if repeat==0:
                    ax.annotate('day %i ' % (day+1),
                                (-3, -shift*day), ha='right',
                                xycoords='data', 
                                color = pt.copper(1-day/resp.shape[1]))
                if day==0:
                    ax.annotate('#%i' % (repeat+1),
                                (tshift*repeat, -shift*(resp.shape[1]-.5)), 
                                va='top', ha='center', xycoords='data')

        pt.draw_bar_scales(ax, 
                            loc='top-right',
                Ybar=0.5, Ybar_label='0.5$\Delta$F/F', 
                Xbar=10, Xbar_label='10s', color='k')
        
# pt.save(fig, 'Desktop/plasticity', '.png',
#         transparent=False)

# %%


