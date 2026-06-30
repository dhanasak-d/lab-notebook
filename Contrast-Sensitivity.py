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


dataset = scan_folder_for_NWBfiles(\
        os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026'),
            for_protocols=['contrast-sensitivity']
            )

# %%
dFoF_parameters = dict(\
    roi_to_neuropil_fluo_inclusion_factor=0., # no factor here
    neuropil_correction_factor = 0.7,
    method_for_F0 = 'sliding_percentile',
    percentile=5., # percent
    sliding_window = 5*60, # seconds
)

quantity = 'dFoF'

# %%
dataset
filename = dataset['files'][0]
data = Data(filename)
data.build_dFoF(**dFoF_parameters)
for p in data.protocols[1:]:
    print(p)

    ep = EpisodeData(data, protocol_name=p)
    print(ep.varied_parameters)
# ep = EpisodeData(data, verbose=True, 
#                  protocol_name=['contrast-sensitivity'])
# %%


def cell_sensitivity_example_fig(EPISODES,
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

    EPISODES = EpisodeData(data,
                           quantities=[quantity],
                           protocol_name='contrast-sensitivity',
                           verbose=False)

    fig, AX = pt.plt.subplots(Nsamples, len(EPISODES.varied_parameters['contrast']), 
                          figsize=(7,7))
    plt.subplots_adjust(right=0.75, left=0.1, top=0.97, bottom=0.05, wspace=0.1, hspace=0.8)
    
    for Ax in AX:
        for ax in Ax:
            ax.axis('off')

    for i, r in enumerate(np.random.choice(np.arange(data.nROIs), 
                                           min([Nsamples, data.nROIs]), replace=False)):

        # SHOW trial-average
        plot_trial_average(EPISODES,
                           quantity=quantity,
                           column_key='contrast',
                           color=color,
                           Ybar_label=1 if quantity=='dFoF' else None,
                           Ybar_label='1dF/F' if quantity=='dFoF' else 'a.u.',
                           Xbar=1., Xbar_label='1s',
                           roiIndex=r,
                           with_stat_test=True,
                           stat_test_props=stat_test_props,
                           with_screen_inset=False,
                           AX=[AX[i]], no_set=False)
        AX[i][0].annotate('roi #%i  ' % (r+1), (0,0), 
                          ha='right', xycoords='axes fraction')

        # SHOW summary angle dependence
        inset = pt.inset(AX[i][-1], (2.2, 0.2, 1.2, 0.8))

        contrasts, y, sy, responsive_contrasts = [], [], [], []
        responsive = False

        for c, contrast in enumerate(EPISODES.varied_parameters['contrast']):

            stats = EPISODES.stat_test_for_evoked_responses(episode_cond=\
                                            EPISODES.find_episode_cond(key='contrast',
                                                                       value=contrast),
                                                            response_args=dict(quantity=quantity,
                                                                               roiIndex=r),
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
fig = cell_sensitivity_example_fig(ep)
# plt.show()
# %%


# %%
# to be a valid dataset:
nMIN_DATAFILES = 2

def process_file(filename, i, c, quantity):

    # to be a valid datafile:
    nMIN_ROIs = 4

    # CELL-dependent calcium pre-processing params 
    dFoF_parameters = get_dFoF_params(c)

    # statistical test for visually-evoked-responses
    if quantity=='Deconvolved':
        stat_test_props=dict(interval_pre=[-1.,-0.0],
                            interval_post=[0.0, 1.0],                                   
                            test='ttest',                                            
                            sign='positive')
    else:
        stat_test_props=dict(interval_pre=[-1.,0],
                            interval_post=[1.,2.],                                   
                            test='ttest',                                            
                            sign='positive')


    response_significance_threshold=5e-2

    print('%i) ' % (i+1), 'analyzing file: %s  [...] ' % filename)
    data = Data(filename, verbose=False)
    protocol_name=[p for p in data.protocols if '8contrast' in p][0]
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

            np.save(os.path.join(summary_folder, 'temp', 
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

DATASET = scan_folder_for_NWBfiles(folder)

# FILTER
# 1) protocol type: contrast sensitivity
cond = np.array([np.sum(['8contrast' in p for p in protocols])\
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
                                    args=(DATASET['files'][cond][i], 
                                            i, c, quantity))
                procs.append(proc)
                proc.start()

            # complete the processes
            for proc in procs:
                proc.join()
    else:
        #####################################
        ###### UN-PARALLELIZED VERSION ######
        for i, f in enumerate(DATASET['files'][cond]):
            process_file(f, i, c, quantity)
        #####################################

    # now that we have stored all datafile outputs
    Sensitivities = []
    for i, f in enumerate(DATASET['files'][cond]):

        if os.path.isfile(os.path.join(summary_folder, 'temp', 
                                        'Sensitivity-%s-%s-%i.npy' % (quantity, c, i))):
            Sensitivity = np.load(os.path.join(summary_folder, 'temp', 
                                        'Sensitivity-%s-%s-%i.npy' % (quantity, c, i)),
                                allow_pickle=True).item()
            Sensitivities.append(Sensitivity)

    # # saving data
    np.save(os.path.join(summary_folder, 'Sensitivities_%s_%s.npy' % (quantity, c)), 
            Sensitivities)
