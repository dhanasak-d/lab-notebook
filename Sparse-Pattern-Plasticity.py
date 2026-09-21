# %% [markdown]
# # Plasticity of Sparse Temporal Sequences

# %%
# general python modules for scientific analysis
import sys, pathlib, os
import numpy as np

sys.path += ['physion/src'] # add src code directory for physion
from physion.utils import plot_tools as pt

from physion.analysis.read_NWB import Data,\
    scan_folder_for_NWBfiles

from physion.analysis.episodes.build import EpisodeData

pt.set_style('dark')

# %%
folder = os.path.join(os.path.expanduser('~'), 
                'DATA', 'Taddy', 'PN_shGrid1-2026', 'NWBs')

notebook_folder =\
    os.path.join(os.path.expanduser('~'), 
        'OneDrive - ICM', 'Lab-Notebook', 'Projects', 'Taddy-GluN3', 'figs')

if not os.path.isdir(os.path.join(folder, 'temp')):
    os.mkdir(os.path.join(folder, 'temp'))

# %%
dataset = scan_folder_for_NWBfiles(\
        os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026'),
            for_protocol='plasticity-60',
            )

# %%

sorted = {}
for v, virus in enumerate(\
        np.unique(dataset['viruses'])):

    print(' - %s:' % virus)
    virus_cond = (virus==dataset['viruses'])

    sorted[virus] = {}
    for s, subject in enumerate(\
            np.unique(dataset['subjects'][virus_cond])):

        sorted[virus][subject] = []
        subject_cond = dataset['subjects'][virus_cond]==subject

        subject_files = [os.path.basename(f) for f in \
            dataset['files'][virus_cond][subject_cond]]

        print('     - %s:' % subject)
        for s in subject_files:
            print('         - [%s](../../Data/%s.md)' %\
                (s.replace('.nwb',''), s.replace('.nwb','')))
        for s in dataset['files'][virus_cond][subject_cond]:
            sorted[virus][subject].append(s)
    print('\n')

# %% [markdown]
# # Analysis Per Session

# %%

dFoF_parameters = dict(\
    roi_to_neuropil_fluo_inclusion_factor=0., # no factor here
    neuropil_correction_factor = 0.5,
    # with_computed_neuropil_fact=True, 
    method_for_F0 = 'sliding_percentile',
    percentile=10., # percent
    sliding_window = 5*60, # seconds
)

quantities = ['dFoF']
# %%

def single_rec(filename,
                nsplit = 3,
                window=[0.3, 0.6]):

    data = Data(filename)
    data.build_dFoF(**dFoF_parameters, verbose=False)
    Ep = EpisodeData(data, prestim_duration=3, 
                protocol_id=0, quantities=['dFoF'])

    fig, AX = pt.figure(axes=(3,1), ax_scale=(1.3, 1.1), top=1.5, right=1.5)
    ntrials = Ep.dFoF.shape[0]
    dn = int(ntrials/nsplit)
    colors = [pt.spring(i/(nsplit-1)) for i in range(nsplit)]
    resps = {'mB':[], 'norm':[], 'raw':[]}
    for i in range(nsplit):
        trials = np.arange(i*dn, (i+1)*dn)
        resp = Ep.dFoF[trials,:,:].mean(axis=(0,1))
        AX[0].plot(Ep.t, resp, color=colors[i])
        respmB = resp-resp[Ep.t<0].mean()
        AX[1].plot(Ep.t, respmB, color=colors[i])
        respN = respmB/respmB[(Ep.t>window[0]) & (Ep.t<window[1])].mean()
        AX[2].plot(Ep.t, respN, color=colors[i])

        pt.annotate(AX[-1], (i+1)*'\n'+'%i-%i' % (trials[0], trials[-1]),
                    (1,1), va='top', color=colors[i])

        cond = (Ep.t>window[1]) & (Ep.t<Ep.time_duration[0])
        resps['raw'].append(resp[cond].mean())
        resps['mB'].append(respmB[cond].mean())
        resps['norm'].append(respN[cond].mean())

    pt.annotate(AX[-1], 'trials:', (1,1), va='top')
    pt.set_plot(AX[0], xlabel='time (s)', title='raw $\\Delta$F/F')
    AX[1].plot(Ep.t, 0*Ep.t, ':')
    pt.set_plot(AX[1], xlabel='time (s)', title='$\\Delta$F/F - baseline')
    AX[2].plot(Ep.t, 1+0*Ep.t, ':')
    AX[2].plot(Ep.t, 0*Ep.t, ':')
    pt.set_plot(AX[2], xlabel='time (s)', title='norm to early [%.1f,%.1f] $\\Delta$F/F' % tuple(window))
    pt.annotate(fig, 'n=%i ROIs ' % data.nROIs, (0,0))
    pt.timestamp(fig)

    return fig, AX, resps

single_rec(dataset['files'][0])
 
# %%
nsplit = 3
resps = {}
for virus in np.unique(dataset['viruses']):
    virus_cond = (virus==dataset['viruses'])
    v = virus.split('+sh')[1] # short name for virus
    resps[v] = {'raw':[], 'mB':[], 'norm':[]}
    for s, subject in enumerate(\
        np.unique(dataset['subjects'][virus_cond])):
        subject_cond = dataset['subjects'][virus_cond]==subject
        subject_files = dataset['files'][virus_cond][subject_cond]
        # resps[virus][subject] = []
        for f in subject_files:
            fn = os.path.basename(f).replace('.nwb','')
            fig_name = 'sparse-plasticity-%s-%s.svg' % (subject, fn)
            print('![](figs/%s)' % fig_name) # for notebook !
            print('\n')
            if 1:
                fig, ax, r = single_rec(f, nsplit=nsplit)
                pt.annotate(fig, '%s - %s - %s' % (v, subject, fn), 
                (0,1), va='top')
                pt.save(fig, notebook_folder, fig_name)
            # resps[virus][subject].append(r)
            for key in r:
                resps[v][key].append(r[key])
# %%
from scipy import stats
window = [0.3, 0.6]
colors = [pt.spring(i/(nsplit-1)) for i in range(nsplit)]

fig, AX = pt.figure((3,1), ax_scale=(1.3,1.))
for v, virus in enumerate(resps):
    for k, key in enumerate(resps[virus]):
        resps[virus][key] = np.array(resps[virus][key])
        for s in range(nsplit):
            pt.scatter([4*v+s],
                   [resps[virus][key][:,s].mean()],
                   sy=[stats.sem(resps[virus][key][:,s])],
                   color=colors[s], ax=AX[k])
        if v==1:
            pt.set_plot(AX[k], title='%s $\\Delta$F/F' % key,
                        xticks=np.arange(2)*4+1,
                        xticks_labels=['%s\n(N=%i)' %\
                             (k, len(resps[virus]['raw'])) for k in resps])
for i in range(nsplit):
    dn = 20
    trials = np.arange(i*dn, (i+1)*dn)
    pt.annotate(AX[-1], (i+1)*'\n'+'%i-%i' % (trials[0], trials[-1]),
                (1,1), va='top', color=colors[i])

pt.timestamp(fig)
fig_name = 'sparse-plasticity-summary.svg'
print('![](figs/%s)' % fig_name) # for notebook !
pt.save(fig, notebook_folder, fig_name)

# %%
# # Analysis Across Days

# %%
def multiple_recs(filenames,
                nsplit = 3,
                window=[0.3, 0.6]):

    fig, AX = pt.figure(axes=(3,1), ax_scale=(1.3, 1.1), top=1.5, right=1.5)
    colors = [pt.spring(i/(len(filenames)-1)) for i in range(len(filenames))]

    resps = {'mB':[], 'norm':[], 'raw':[]}
    for i, filename in enumerate(filenames):
        data = Data(filename)
        data.build_dFoF(**dFoF_parameters, verbose=False)
        Ep = EpisodeData(data, prestim_duration=3, 
                    protocol_id=0, quantities=['dFoF'])

        resp = Ep.dFoF.mean(axis=(0,1))
        AX[0].plot(Ep.t, resp, color=colors[i])
        respmB = resp-resp[Ep.t<0].mean()
        AX[1].plot(Ep.t, respmB, color=colors[i])
        respN = respmB/respmB[(Ep.t>window[0]) & (Ep.t<window[1])].mean()
        AX[2].plot(Ep.t, respN, color=colors[i])

        pt.annotate(AX[-1], (i+1)*'\n'+'day %i' % (i+1),
                    (1,1), va='top', color=colors[i])

        cond = (Ep.t>window[1]) & (Ep.t<Ep.time_duration[0])
        resps['raw'].append(resp[cond].mean())
        resps['mB'].append(respmB[cond].mean())
        resps['norm'].append(respN[cond].mean())

    pt.annotate(AX[-1], 'days:', (1,1), va='top')
    pt.set_plot(AX[0], xlabel='time (s)', title='raw $\\Delta$F/F')
    AX[1].plot(Ep.t, 0*Ep.t, ':')
    pt.set_plot(AX[1], xlabel='time (s)', title='$\\Delta$F/F - baseline')
    AX[2].plot(Ep.t, 1+0*Ep.t, ':')
    AX[2].plot(Ep.t, 0*Ep.t, ':')
    pt.set_plot(AX[2], xlabel='time (s)', title='norm to early [%.1f,%.1f] $\\Delta$F/F' % tuple(window))
    pt.timestamp(fig)

    return fig, AX, resps

multiple_recs(sorted['CamKII-Cre+shGrid1']['016M'])
# %%
resps = {}

for virus in sorted:
    resps[virus] = {}
    for mouse in sorted[virus]:
        fig, AX, resp = multiple_recs(sorted[virus][mouse])
        pt.annotate(fig, 'mouse: %s - %s ' % (
                    mouse, virus.split('+sh')[1]), (0,1), va='top')
        fig_name = 'sparse-plasticity-across-days-%s.svg' % mouse
        print('![](figs/%s)     \n' % fig_name) # for notebook !
        pt.save(fig, notebook_folder, fig_name)
        for key in resp:
            if key in resps[virus]:
                resps[virus][key].append(resp[key])
            else:
                resps[virus][key] = [resp[key]]

# %%
from scipy import stats
window = [0.3, 0.6]
colors = [pt.spring(i/(nsplit-1)) for i in range(nsplit)]

fig, AX = pt.figure((3,1), ax_scale=(1.3,1.), right=2.2)
for v, virus in enumerate(resps):
    for k, key in enumerate(resps[virus]):
        resps[virus][key] = np.array(resps[virus][key])
        for s in range(nsplit):
            pt.scatter([4*v+s],
                   [resps[virus][key][:,s].mean()],
                   sy=[stats.sem(resps[virus][key][:,s])],
                   color=colors[s], ax=AX[k])
        if v==1:
            pt.set_plot(AX[k], title='%s $\\Delta$F/F' % key,
                        xticks=np.arange(2)*4+1,
                        xticks_labels=['%s\n(N=%i)' %\
                             (k.split('+sh')[1], len(resps[virus]['raw'])) for k in resps])

for i in range(3):
    pt.annotate(AX[-1], (i+1)*'\n'+'day %i' % (i+1),
                (1,1), va='top', color=colors[i])
pt.timestamp(fig)
fig_name = 'sparse-plasticity-across-summary.svg'
print('![](figs/%s)' % fig_name) # for notebook !
pt.save(fig, notebook_folder, fig_name)
# %%
