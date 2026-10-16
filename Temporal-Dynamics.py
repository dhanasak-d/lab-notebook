# %% [markdown]
# # Build Deconvolved Responses across Conditions

# %%
import os

import os, sys , shutil 
import multiprocessing
import numpy as np

sys.path += ['physion/src']

from physion.analysis.read_NWB\
                import scan_folder_for_NWBfiles, Data
from physion.analysis.episodes.build import EpisodeData
import physion.utils.plot_tools as pt

parallelized, debug = False, False 

dFoF_parameters = dict(\
        roi_to_neuropil_fluo_inclusion_factor=0., # no factor here
        neuropil_correction_factor = 0.5,
        method_for_F0 = 'sliding_percentile',
        percentile=20., # percent
        sliding_window = 5*60, # seconds
        )

folder = os.path.join(os.path.expanduser('~'), 
                'DATA', 'Taddy', 'PN_shGrid1-2026', 'NWBs')

notebook_folder =\
    os.path.join(os.path.expanduser('~'), 
        'OneDrive - ICM', 'Lab-Notebook', 'Projects', 'Taddy-GluN3', 'figs')

if not os.path.isdir(os.path.join(folder, 'temp')):
    os.mkdir(folder.replace('NWBs', 'temp'))

pt.set_style('dark')

# %%
 
PROTOCOLS = [\
    'tuning-low-contrast',
    'tuning-mid-contrast',
    'tuning-high-contrast',
    'contrast-sensitivity']

dataset = scan_folder_for_NWBfiles(\
        os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026'),
            for_protocols=PROTOCOLS
            )
# %%
def process_file(filename):

    # to be a valid datafile:
    nMIN_ROIs = 4

    TAU_DECONVOLUTION = 0.8

    # statistical test for visually-evoked-responses
    stat_test_props=dict(interval_pre=[-1,0],
                         interval_post=[0.,1],                                   
                         test='ttest',                                            
                         sign='positive')

    response_significance_threshold=1e-3

    print('    analyzing file: %s  [...] ' % filename)
    data = Data(filename, verbose=False)
    data.build_dFoF(**dFoF_parameters, verbose=False)
    data.build_Deconvolved(Tau=TAU_DECONVOLUTION)

    quantities = ['Deconvolved']
    if 'Running-Speed' in data.nwbfile.acquisition:
        quantities += ['running']

    if data.nROIs>=nMIN_ROIs:

        for protocol in PROTOCOLS:

            # process episodes
            Episodes = EpisodeData(data, 
                                    quantities=quantities,
                                    protocol_name=protocol,
                                    verbose=debug,
                                    dt_sampling=10)

            significant = []
            for r in range(data.nROIs):
                    cell_resp = Episodes.pre_post_statistics(stat_test_props=stat_test_props,
                                                             response_significance_threshold=\
                                                                response_significance_threshold,
                                                             response_args=dict(quantity='Deconvolved',
                                                                                index=r),
                                                            verbose=debug)
                    significant.append(cell_resp['significant'])

            stim_cond = (Episodes.t>0) & (Episodes.t<Episodes.time_duration[0])
            Response = {
                    't':Episodes.t,
                    'Deconvolved':Episodes.Deconvolved[:,:,:],
                    'running':Episodes.running[:,stim_cond].mean(axis=1),
                    'significant':np.array(significant, dtype=bool),
                    'nROIs_original': data.original_nROIs,
                    'nROIs_final': data.nROIs,
                    'subject':data.nwbfile.subject.subject_id,
            }

            fn = os.path.join(folder.replace('NWBs', 'temp'),
                                 'tempResponse-%s-%s.npy' %\
                                      (protocol, 
                                       os.path.basename(filename).replace('.nwb','')))
            np.save(fn, Response)
            print('      [v] --> %s - included, n=%i ROIs ' % (protocol, data.nROIs))
            print('                ', fn)

process_file(dataset['files'][0])
# %%
for i, f in enumerate(dataset['files']):
    print(dataset['viruses'][i])
    process_file(f)

# %%
# check running
RUNNING_THRESHOLD = 0.1
running = {}

fig, ax = pt.figure(ax_scale=(1.,1.2))
for v, virus in enumerate(\
        np.unique(dataset['viruses'])[::-1]):
    running[virus] = []
    virus_cond = (virus==dataset['viruses'])

    for protocol in PROTOCOLS:
        summary[virus][protocol] = {'frac_running':[]}
        for filename in dataset['files'][virus_cond]:
            fn = os.path.join(folder.replace('NWBs', 'temp'),
                                    'tempResponse-%s-%s.npy' %\
                                        (protocol, 
                                        os.path.basename(filename).replace('.nwb','')))
            if os.path.isfile(fn):
                responses = np.load(fn, allow_pickle=True).item()
                run_ep = responses['running'].flatten()>RUNNING_THRESHOLD
                running[virus].append(
                    100.*np.sum(run_ep)/len(run_ep)
                )
    ax.bar([v],
           [np.mean(running[virus])],
           yerr=[stats.sem(running[virus])],
           color='C%i' % v, label=virus)
pt.set_plot(ax, xticks=[], ylabel='% running ep. ')
ax.legend(loc=(1,1))
# %%
protocol = 'tuning-high-contrast'
filename = dataset['files'][0]

fig, AX = pt.figure(axes=(3,1), )
for ax, protocol in zip(AX, PROTOCOLS):
    fn = os.path.join(folder.replace('NWBs', 'temp'),
                            'tempResponse-%s-%s.npy' %\
                                (protocol, 
                                os.path.basename(filename).replace('.nwb','')))
    responses = np.load(fn, allow_pickle=True).item()
    significant = np.sum(responses['significant'], axis=1)
    ax.plot(responses['t'], responses['Deconvolved'][:,significant,:].mean(axis=(0,1)))
pt.set_common_ylims(AX)

# %%
from scipy.ndimage import gaussian_filter1d

def build_summary(over='sessions',
                  smoothing=7):

    summary = {'over':over}
    for v, virus in enumerate(\
            np.unique(dataset['viruses'])):

        print(' - %s:' % virus)
        summary[virus] = {}
        virus_cond = (virus==dataset['viruses'])

        #
        protocol='contrast-sensitivity' 
        summary[virus][protocol] = {'frac_resp':[]}
        for filename in dataset['files'][virus_cond]:
            fn = os.path.join(folder.replace('NWBs', 'temp'),
                                    'tempResponse-%s-%s.npy' %\
                                        (protocol, 
                                        os.path.basename(filename).replace('.nwb','')))
            if os.path.isfile(fn):
                responses = np.load(fn, allow_pickle=True).item()
                significant = responses['significant']
                summary[virus][protocol]['frac_resp'].append(
                    100.*np.sum(significant, axis=0)/significant.shape[0]
                )
        summary[virus][protocol]['frac_resp'] =\
                np.array(summary[virus][protocol]['frac_resp'])


        for protocol in ['tuning-low-contrast', 
                         'tuning-mid-contrast',
                         'tuning-high-contrast']:
            summary[virus][protocol] = {
                'frac_resp':[], 
                'resp':[]}


            for filename in dataset['files'][virus_cond]:
                fn = os.path.join(folder.replace('NWBs', 'temp'),
                                        'tempResponse-%s-%s.npy' %\
                                            (protocol, 
                                            os.path.basename(filename).replace('.nwb','')))
                if os.path.isfile(fn):
                    responses = np.load(fn, allow_pickle=True).item()
                    significant = np.sum(responses['significant'], axis=1)
                    summary[virus][protocol]['frac_resp'].append(
                        100.*np.sum(significant>0)/len(significant)
                    )
                    if summary['over']=='sessions':
                        summary[virus][protocol]['resp'].append(
                            gaussian_filter1d(
                                responses['Deconvolved'].mean(axis=(0,1)),
                                smoothing))
                    elif summary['over']=='ROIs':
                        for roi in range(responses['Deconvolved'].shape[1]):
                            summary[virus][protocol]['resp'].append(
                                gaussian_filter1d(
                                    responses['Deconvolved'][:,roi,:].mean(axis=0),
                                    smoothing))
                    else:
                        print('choose between ROIs and sessions')

            summary[virus][protocol]['resp'] =\
                    np.array(summary[virus][protocol]['resp'])
            summary[virus][protocol]['frac_resp'] =\
                    np.array(summary[virus][protocol]['frac_resp'])
    summary['t'] = responses['t']
    return summary
summary = build_summary('sessions')

# %%

from scipy import stats
fig, ax = pt.figure(ax_scale=(1.4,1.1), right=3)
inset = pt.inset(ax, [1.4, 0., 0.3, 1.])
contrast = np.linspace(0.05, 1., 8)
viruses = np.unique(dataset['viruses'])[::-1]
protocol = 'contrast-sensitivity'
slopes = {}
for v, virus in enumerate(viruses):
    slopes[virus] = []
    for c in range(len(contrast)):
        ax.bar([.04*v+contrast[c]],
            [np.mean(summary[virus][protocol]['frac_resp'][:,c])],
            yerr=[stats.sem(summary[virus][protocol]['frac_resp'][:,c])],
            color='C%i' % v, width=.06)
    resp = summary[virus][protocol]['frac_resp']
    for r in range(resp.shape[0]):
        slopes[virus].append(np.polyfit(contrast, resp[r,:], 1)[0])

    # pt.violin(slopes[virus], x=[v],
    #            color='C%i' % v, ax=inset)
    pt.scatter([v], [np.mean(slopes[virus])], 
               sy=[stats.sem(slopes[virus])],
               color='C%i' % v, ax=inset)

inset.annotate(
    pt.from_pval_to_star(
        stats.ttest_ind(slopes['CamKII-Cre+shScramble'],
               slopes['CamKII-Cre+shGrid1']).pvalue),
               (0.5,inset.get_ylim()[1]), ha='center')

pt.set_plot(inset,
            xticks=[], xlim=[-0.5,1.5],
            ylabel='slope\n(%resp./contr.)')
pt.set_plot(ax,
            xticks=[contrast[0], 1],
            xticks_labels=['5%', '100%'],
            xlabel='contrast', ylabel='% responsive')


# %%
from scipy import stats

def make_fig(summary,
        peak_window = [0, 0.5],
        inset_lims = [0.1,0.9],
        baseline_window = [-0.7,-0.2]):

    fig, AX = pt.figure(axes=(3,1), 
                        wspace=2.,
                        ax_scale=(0.8,1.1), right=3)

    peak_cond =\
          (summary['t']>peak_window[0]) & (summary['t']<peak_window[1])

    baseline_cond = (summary['t']>baseline_window[0]) &\
            (summary['t']<baseline_window[1])

    viruses = np.unique(dataset['viruses'])[::-1]
    insets = [] 
    for v, virus in enumerate(viruses):
        for p, protocol in enumerate(PROTOCOLS[:3]):
            baseline =\
                np.mean(summary[virus][protocol]['resp'][:,baseline_cond], axis=1)
            # remove baseline
            resp = (summary[virus][protocol]['resp'].T-baseline.T).T
            # pt.plot(summary['t'], np.mean((resp_mB.T/norm.T).T, axis=0),
            pt.plot(summary['t'], 100.*np.mean(resp, axis=0),
                sy=stats.sem(resp, axis=0),
                color='C%i' % v, ax=AX[p])
            # normalize
            if v==0:
                insets.append(pt.inset(AX[p], [0.6,0.55,0.7,0.6]))

            print(p, insets[p])
            norm = np.max(resp.mean(axis=0))
            pt.plot(summary['t'], np.mean(resp/norm, axis=0),
                sy=stats.sem(resp/norm, axis=0),
                color='C%i' % v, ax=insets[p])

    pt.set_common_xlims(AX, lims=[-1.5,4])
    pt.set_common_ylims(AX)

    for ax in insets:
        ax.plot(inset_lims, [0,0], ':', lw=0.4, color='C0')
        ax.annotate('0', (inset_lims[0], 0), ha='right', va='center')
        ax.plot(inset_lims, [1,1], ':', lw=0.4, color='C0')
        ax.annotate('1', (inset_lims[0], 1), ha='right', va='center')
        ax.plot(inset_lims[0]+np.zeros(2), [0,1], '-', color='C0')
        ax.plot(inset_lims[1]-np.arange(2)*0.3, [-0.1,-0.1], '-', color='C0')
        ax.annotate('0.3s', (inset_lims[1], -0.2), ha='right', va='top')
        ax.axis('off')
    pt.set_common_ylims(insets)
    pt.set_common_xlims(insets, inset_lims)

    for v, virus in enumerate(viruses):
        if summary['over']=='sessions':
            averaging = '(N=%i sessions)' % summary[virus]['tuning-low-contrast']['resp'].shape[0]
        elif summary['over']=='ROIs':
            averaging = '(n=%i ROIs)' % summary[virus]['tuning-low-contrast']['resp'].shape[0]
        pt.annotate(ax, 2*v*'\n'+virus.split('+')[1]+ '\n%s' % averaging, 
                (1.3,1.), va='top', color='C%i' % v)

    for ax, title in zip(AX, ['20%', '60%', '100%']):
        pt.set_plot(ax, 
                    yticks_labels=([] if ax!=AX[0] else None),
                    ylabel=('$\delta$ deconv. (10$^2$ a.u.)     ' if ax==AX[0] else ''),
                    xlabel='time (s)')
        pt.annotate(ax, 'c=%s' % title, (0,1), ha='center')

    return fig, AX

summary = build_summary('sessions')
make_fig(summary)
summary = build_summary('ROIs')
make_fig(summary)

# %%
from scipy import stats

summary = build_summary('sessions')

fig, ax = pt.figure(ax_scale=(1.3,1.1), right=10)
inset = pt.inset(ax, [1.6, 0., 0.3, 1.])
inset2 = pt.inset(ax, [2.4, 0., 0.3, 1.])

viruses = np.unique(dataset['viruses'])[::-1]
contrasts = [20, 60 , 100]
slopes, means = {}, {}
for v, virus in enumerate(viruses):
    slopes[virus] = [
        np.polyfit(contrasts, 
        [summary[virus][protocol]['frac_resp'][n] for protocol in PROTOCOLS[:3]], 1)[0]
        for n in range(summary[virus][protocol]['frac_resp'].shape[0])]
    means[virus] = [
        np.mean([summary[virus][protocol]['frac_resp'][n] for protocol in PROTOCOLS[:3]])\
            for n in range(summary[virus][protocol]['frac_resp'].shape[0])]
    for p, protocol in enumerate(PROTOCOLS[:3]):
        ax.bar([v+3*p],
            [np.mean(summary[virus][protocol]['frac_resp'])],
            yerr=[stats.sem(summary[virus][protocol]['frac_resp'])],
            color='C%i' % v)
    
    pt.scatter([v], [np.mean(slopes[virus])], 
               sy=[stats.sem(slopes[virus])],
               color='C%i' % v, ax=inset)
    pt.scatter([v], [np.mean(means[virus])], 
               sy=[stats.sem(means[virus])],
               color='C%i' % v, ax=inset2)

inset.annotate(
    'p=%.3f' % stats.ttest_ind(slopes['CamKII-Cre+shScramble'],
                        slopes['CamKII-Cre+shGrid1']).pvalue,
               (0.5,inset.get_ylim()[1]), ha='center')
inset2.annotate(
    'p=%.3f' % stats.ttest_ind(means['CamKII-Cre+shScramble'],
                        means['CamKII-Cre+shGrid1']).pvalue,
    # pt.from_pval_to_star(
    #     stats.ttest_ind(slopes['CamKII-Cre+shScramble'],
    #            slopes['CamKII-Cre+shGrid1']).pvalue),
               (0.5,inset2.get_ylim()[1]), ha='center')

pt.set_plot(inset, xticks=[], xlim=[-0.5,1.5],
            ylabel='slope\n(%resp./contr.)')
pt.set_plot(inset2, xticks=[], xlim=[-0.5,1.5],
            ylabel='mean %resp.')

ylim = ax.get_ylim()
for p, protocol in enumerate(PROTOCOLS[:3]):
    pt.annotate(ax,
        pt.from_pval_to_star(
            stats.ttest_ind(
                *[summary[virus][protocol]['frac_resp'] for virus in viruses]
                ).pvalue), (0.5+3*p, ylim[1]), ha='center', xycoords='data')
for v, virus in enumerate(viruses):
    if summary['over']=='sessions':
        averaging = '(N=%i sessions)' % summary[virus]['tuning-low-contrast']['resp'].shape[0]
    elif summary['over']=='ROIs':
        averaging = '(n=%i ROIs)' % summary[virus]['tuning-low-contrast']['resp'].shape[0]

    pt.annotate(inset2, 2*v*'\n'+virus.split('+')[1]+ '\n%s' % averaging, 
            (1,1), va='top', color='C%i' % v)

pt.set_plot(ax,
            xticks=0.5+np.arange(3)*3,
            xticks_labels=['20%', '60%', '100%'],
            xlabel='contrast', ylabel='% responsive')


# %%
if __name__=='__main__':

    from physion.assembling.dataset import read_spreadsheet

    cpus = multiprocessing.cpu_count()-1 # leaving 1 cpu for the rest

    # temporary folder for parallelization
    os.makedirs(os.path.join(summary_folder, 'temp'), exist_ok=True)

    Nstart = 0
    Nend = len(datasets)

    for n in range(Nstart, Nend):

        c = list(datasets.keys())[n]

        table = datasets[c]['datafolder'].replace('NWBs', 'DataTable.xlsx')

        dataset_table, subjects_table, analysis =\
                read_spreadsheet(table, get_metadata_from='table')
        print()
        print()
        print('=================================================================')
        print('-----------------------------------------------------------------')
        print('------- %i) computing : %s ' % (n+1, c))
        print('-----------------------------------------------------------------')
        print()

        DATASET = scan_folder_for_NWBfiles(datasets[c]['datafolder'])
        
        # FILTER
        # 1) protocol type: orientation tuning
        cond = np.array([np.sum(['8orientation' in p for p in protocols])\
                        for protocols in DATASET['protocols']], dtype=bool)
        # 2) age condition
        if datasets[c]['age_interval'] is not None:
            cond = cond &\
                (DATASET['ages']>=datasets[c]['age_interval'][0]) &\
                (DATASET['ages']<=datasets[c]['age_interval'][1])


        if len(DATASET['files'][cond])>nMIN_DATAFILES:

            if parallelized:
                ################################################
                ###    parallelization here !   #################
                ################################################

                nruns = int(len(DATASET['files'][cond])/cpus)+1

                for r in range(nruns):
                    i0 = r*cpus
                    imax = np.min([i0+cpus, len(DATASET['files'][cond])]) 
                    print(' - running set of files %i:%i' % (i0, imax))

                    # start the processes
                    procs = []
                    for i in range(i0,imax):
                        proc = multiprocessing.Process(\
                                            target=process_file, 
                                            args=(DATASET['files'][cond][i], i, c))
                        procs.append(proc)
                        proc.start()

                    # complete the processes
                    for proc in procs:
                        proc.join()

            else:
                #####################################
                ###### UN-PARALLELIZED VERSION ######
                for i, f in enumerate(DATASET['files'][cond]):
                    process_file(f, i, c)
                #####################################

            # now that we have stored all datafile outputs
            Responses = []
            for i, f in enumerate(DATASET['files'][cond]):

                if os.path.isfile(os.path.join(summary_folder, 'temp', 
                                              'tempResponse-%s-%i.npy' % (c, i))):
                    Response = np.load(os.path.join(summary_folder, 'temp', 
                                                'tempResponse-%s-%i.npy' % (c, i)),
                                        allow_pickle=True).item()
                    Responses.append(Response)

            # # saving data
            np.save(os.path.join(summary_folder, 'Deconvolved_%s.npy' % c), 
                    Responses)

        else:
            print()
            print('   [!!]   DATASET NOT LARGE ENOUGH   [!!] ')
            print('               only N=%i sessions available' %\
                                        len(DATASET['files'][cond]))
            print('   [!!]   DATASET not analyzed       [!!] ')
            print()

        print('-----------------------------------------------------------------')
        print('=================================================================')
    # shutil.rmtree(os.path.join(summary_folder, 'temp'))

# %%
if False:

    import sys, os
    sys.path += ['physion/src']
    import physion.utils.plot_tools as pt
    from scipy import stats
    from physion.dataviz.episodes.temporal_dynamics\
          import plot_response_dynamics

    summary_folder = os.path.join(os.path.expanduser('~'), 
                                'CURATED', 'Cibele', 'summary')
 
    fig, ax = plot_response_dynamics(\
                        ['PV-cells_WT_Adult_V1_contrast-1.0', 
                         'PV-cells_WT_Adult_V1_contrast-0.5'],
                         colors=['tab:red', 'lightpink'],
                        # average_by='ROIs',
                        average_by='subjects',
                        path=summary_folder)
    
    # %%
