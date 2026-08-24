# %% [markdown]
# # Contrast Sensitivity

# %%
# general python modules for scientific analysis
import sys, pathlib, os
import numpy as np

sys.path += ['physion/src'] # add src code directory for physion
from physion.utils import plot_tools as pt

from physion.analysis.read_NWB import Data,\
    scan_folder_for_NWBfiles
from physion.analysis.episodes.build import EpisodeData
from physion.dataviz.episodes.trial_average import plot as plot_trial_average

PROTOCOL = 'contrast-sensitivity'

folder = os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026', 'NWBs', 'Orientation-Contrast')
if not os.path.isdir(os.path.join(folder, 'temp')):
    os.mkdir(os.path.join(folder, 'temp'))

# %%
dataset = scan_folder_for_NWBfiles(
        folder, for_protocols=[PROTOCOL])

# %%
dFoF_parameters = dict(\
    roi_to_neuropil_fluo_inclusion_factor= 0.0, # no factor here
    neuropil_correction_factor = 0.5,
    method_for_F0 = 'sliding_percentile',
    percentile=5., # percent
    sliding_window = 5*60, # seconds
)

quantity = 'dFoF'
# %%

def cell_sensitivity_example_fig(filename,
                                 quantity=quantity,
                                stat_test_props = dict(interval_pre=[-1,0], 
                                                       interval_post=[0,1],
                                                       test='ttest',
                                                       sign='positive'),
                                response_significance_threshold = 0.01,
                                Nsamples = 10, # how many cells we show
                                 color='k',
                                seed=10):
    
    np.random.seed(seed)
    
    data = Data(filename)
    data.build_dFoF(**dFoF_parameters)
    if quantity=='Deconvolved':
        data.build_Deconvolved()

    ep = EpisodeData(data,
                           quantities=[quantity],
                           protocol_name=PROTOCOL,
                           verbose=False)
    print(ep.varied_parameters)
    fig, AX = pt.plt.subplots(Nsamples, 
                              len(ep.varied_parameters['contrast']), 
                          figsize=(7,7))
    pt.plt.subplots_adjust(right=0.75, left=0.1, top=0.97, bottom=0.05, wspace=0.1, hspace=0.8)
    
    for Ax in AX:
        for ax in Ax:
            ax.axis('off')

    for i, r in enumerate(np.random.choice(np.arange(data.nROIs), 
                                           min([Nsamples, data.nROIs]), replace=False)):

        # SHOW trial-average
        plot_trial_average(ep,
                           quantity=quantity,
                           column_key='contrast',
                           color=color,
                           Ybar=1 if quantity=='dFoF' else 0,
                           Ybar_label='1dF/F' if quantity=='dFoF' else 'a.u.',
                           Xbar=1., Xbar_label='1s',
                           index=r,
                           with_stat_test=True,
                           stat_test_props=stat_test_props,
                           with_screen_inset=True,
                           AX=[AX[i]], no_set=False)
        AX[i][0].annotate('roi #%i  ' % (r+1), (0,0), 
                          ha='right', xycoords='axes fraction')

        # SHOW summary angle dependence
        inset = pt.inset(AX[i][-1], (2.2, 0.2, 1.2, 0.8))

        contrasts, y, sy, responsive_contrasts = [], [], [], []
        responsive = False

        for c, contrast in enumerate(ep.varied_parameters['contrast']):

            stats = ep.stat_test_for_evoked_responses(episode_cond=\
                                            ep.find_episode_cond(key='contrast',
                                                                       value=contrast),
                                                            response_args=dict(quantity=quantity,
                                                                               index=r),
                                                            **stat_test_props)

            contrasts.append(contrast)
            y.append(np.mean(stats.y-stats.x))    # means "post-pre"
            sy.append(np.std(stats.y-stats.x))    # std "post-pre"

            if stats.significant(threshold=response_significance_threshold):
                responsive = True
                responsive_contrasts.append(contrast)

        pt.scatter(contrasts, np.array(y), 
                   sy=np.array(sy), ax=inset, ms=1, lw=1, color=color)
        inset.plot(contrasts, 0*np.array(contrasts), 'k:', lw=0.5)
        if quantity=='Deconvolved':
            inset.set_ylabel('$\\delta$ Deconv.    ', fontsize=7)
        else:
            inset.set_ylabel('$\\delta$ $\\Delta$F/F     ', fontsize=7)
        inset.set_xticks([0,1])
        #inset.set_xticklabels(['%i'%a if (i%2==0) else '' for i, a in enumerate(contrasts)], fontsize=7)
    inset.set_xlabel('contrast', fontsize=7)

    return fig

# iSession = 0 # session index
fig = cell_sensitivity_example_fig(dataset['files'][-1])

# %%
data = Data(dataset['files'][-1])
# %%

from physion.analysis.protocols.contrast_sensitivity import\
                        compute_sensitivity_per_cells

# to be a valid dataset:
nMIN_DATAFILES = 2

def process_file(filename, i, c, quantity):

    # to be a valid datafile:
    nMIN_ROIs = 4

    # statistical test for visually-evoked-responses
    stat_test_props=dict(interval_pre=[-1.,-0.0],
                        interval_post=[0.0, 1.0],                                   
                        test='ttest',                                            
                        sign='positive')

    response_significance_threshold=5e-2

    print('%i) ' % (i+1), 'analyzing file: %s  [...] ' % filename)
    data = Data(filename, verbose=False)
    protocol_name=[p for p in data.protocols if PROTOCOL in p][0]
    data.build_dFoF(**dFoF_parameters, verbose=False)

    if data.nROIs>=nMIN_ROIs:

        if quantity=='Deconvolved':
            data.build_Deconvolved()

        # try:
        if True:
            Episodes = EpisodeData(data, 
                                    quantities=[quantity],
                                    protocol_name=protocol_name, 
                                    verbose=False)

            Sensitivity = compute_sensitivity_per_cells(data, Episodes, 
                                                        quantity=quantity,
                                                        stat_test_props=stat_test_props, 
                                                        response_significance_threshold = response_significance_threshold, 
                                                        angle = 0)

            Sensitivity['datafile'] = filename
            Sensitivity['nROIs_original'] = data.original_nROIs
            Sensitivity['nROIs_final'] = data.nROIs
            Sensitivity['subject'] = data.nwbfile.subject.subject_id

            np.save(os.path.join(folder, 'temp', 
                                 'Sensitivity-%s-%s-%i.npy' % (quantity, c, i)),
                    Sensitivity)
            print('      [v] --> included, n=%i ROIs ' % data.nROIs)
        # except BaseException as be:
        #     print('                        [-------------------------------]')
        #     print(be)
        #     print()
        #     print('      [X] --> discarded, problem in datafile, CHECK [!!]')
        #     print('                        [-------------------------------]')

    else:
        print('      [X] --> discarded, n=%i ROIs ' % data.nROIs)

# %%

groups = {
    'shRNA':{'virus':'CamKII-Cre+shGrid1'},
    'scramble':{'virus':'CamKII-Cre+shScramble'}
}

parallelized = False

for c in groups:

    # FILTER
    # 1) protocol type: contrast sensitivity
    cond = np.array([np.sum([PROTOCOL in p for p in protocols])\
                    for protocols in dataset['protocols']], dtype=bool)
    # 2) virus
    cond = cond & (dataset['viruses']==groups[c]['virus'])

    if len(dataset['files'][cond])>nMIN_DATAFILES:

        if parallelized:

            ################################################
            ###    parallelization here !   #################
            ################################################
            import multiprocessing
            cpus = multiprocessing.cpu_count()
            nruns = int(len(dataset['files'][cond])/cpus)+1

            for r in range(nruns):
                i0 = r*cpus
                imax = np.min([i0+cpus, len(dataset['files'][cond])]) 
                print(' - running set of files %i:%i' % (i0, imax))

                # start the processes
                procs = []
                for i in range(i0,imax):
                    proc = multiprocessing.Process(\
                                        target=process_file, 
                                        args=(dataset['files'][cond][i], 
                                                i, c, quantity))
                    procs.append(proc)
                    proc.start()

                # complete the processes
                for proc in procs:
                    proc.join()
        else:
            #####################################
            ###### UN-PARALLELIZED VERSION ######
            for i, f in enumerate(dataset['files'][cond]):
                process_file(f, i, c, quantity)
            #####################################

        # now that we have stored all datafile outputs
        Sensitivities = []
        for i, f in enumerate(dataset['files'][cond]):

            if os.path.isfile(os.path.join(folder, 'temp', 
                                            'Sensitivity-%s-%s-%i.npy' % (quantity, c, i))):
                Sensitivity = np.load(os.path.join(folder, 'temp', 
                                            'Sensitivity-%s-%s-%i.npy' % (quantity, c, i)),
                                    allow_pickle=True).item()
                Sensitivities.append(Sensitivity)

        # # saving data
        np.save(os.path.join(folder, 'Sensitivities_%s_%s.npy' % (quantity, c)), 
                Sensitivities)




# %%
from physion.analysis.protocols.contrast_sensitivity\
        import plot_contrast_sensitivity, plot_contrast_responsiveness


fig, ax = plot_contrast_sensitivity(\
                        ['%s_scramble' % quantity, 
                         '%s_shRNA' % quantity],
                          average_by='ROIs',
                        path=os.path.join(folder))

fig, ax = plot_contrast_sensitivity(\
                        ['%s_scramble' % quantity, 
                         '%s_shRNA' % quantity],
                          average_by='sessions',
                        path=os.path.join(folder))

# %%