# %% [markdown]
# # Orientation Tuning


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

# %%
PROTOCOLS = [\
    'tuning-low-contrast',
    'tuning-mid-contrast',
    'tuning-high-contrast']

folder = os.path.join(os.path.expanduser('~'), 
                'DATA', 'Taddy', 'PN_shGrid1-2026', 'NWBs')

if not os.path.isdir(os.path.join(folder, 'temp')):
    os.mkdir(os.path.join(folder, 'temp'))

dataset = scan_folder_for_NWBfiles(\
        os.path.join(os.path.expanduser('~'), 
            'DATA', 'Taddy', 'PN_shGrid1-2026'),
            for_protocols=PROTOCOLS
            )

quantity = 'dFoF'

# %%
dFoF_parameters = dict(\
    roi_to_neuropil_fluo_inclusion_factor=0., # no factor here
    neuropil_correction_factor = 0.7,
    method_for_F0 = 'sliding_percentile',
    percentile=5., # percent
    sliding_window = 5*60, # seconds
)

# %%
def process_file(filename, i, c, PROTOCOL, quantity):

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
    data.build_dFoF(**dFoF_parameters, verbose=False)

    # FIX: Temporarily disable running_speed to bypass the resampling unpack bug
    quantities = ['dFoF']
    # if 'Running-Speed' in data.nwbfile.acquisition:
    #     quantities += ['running_speed']

    if data.nROIs>=nMIN_ROIs:

        if quantity=='Deconvolved':
            data.build_Deconvolved()

        # try:
        if True:
            Episodes = EpisodeData(data, 
                                    quantities=quantities,
                                    protocol_name=PROTOCOL,
                                    verbose=False)

            Tuning = compute_tuning_response_per_cells(data, Episodes, 
                                                        quantity=quantity,
                                                        stat_test_props = stat_test_props, 
                                                        response_significance_threshold =\
                                                            response_significance_threshold, 
                                                        contrast = Episodes.contrast[0],
                                                        verbose=True)
            Tuning['datafile'] = filename
            Tuning['nROIs_original'] = data.original_nROIs
            Tuning['nROIs_final'] = data.nROIs
            Tuning['nROIs_responsive'] = np.sum(Tuning['significant_ROIs'])
            Tuning['subject'] = data.nwbfile.subject.subject_id

            np.save(os.path.join(folder, 'temp', 
                                 'Tuning-%s-%i.npy' % (c, i)),
                    Tuning)
            print('      [v] --> included, n=%i ROIs ' % data.nROIs)
        # except BaseException as be:
        #     print('                        [-------------------------------]')
        #     print(be)
        #     print()
        #     print(filename)
        #     print('nROIs=%i' % data.nROIs, ', protocols=%s' % data.protocols) 
        #     print(Episodes.varied_parameters)
        #     print('      [X] --> discarded, problem in datafile, CHECK [!!]')
        #     print('                        [-------------------------------]')

    else:
        print('      [X] --> discarded, n=%i ROIs ' % data.nROIs)

# %%
from physion.analysis.protocols.orientation_tuning import\
    compute_tuning_response_per_cells

groups = {
    'shRNA':{'virus':'CamKII-Cre+shGrid1'},
    'scramble':{'virus':'CamKII-Cre+shScramble'}
}

parallelized = False
nMIN_DATAFILES = 2

for g in groups:

    for PROTOCOL in PROTOCOLS:

        c = '%s_%s' % (g, PROTOCOL)

        # FILTER
        # 1) protocol type: contrast sensitivity
        cond = np.array([np.sum([PROTOCOL in p for p in protocols])\
                        for protocols in dataset['protocols']], dtype=bool)
        # 2) virus
        cond = cond & (dataset['viruses']==groups[g]['virus'])

        if len(dataset['files'][cond])>nMIN_DATAFILES:

            if parallelized:
                ################################################
                ###    parallelization here !   #################
                ################################################
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
                                            args=(dataset['files'][cond][i], i, c, 
                                                  PROTOCOL, quantity))
                        procs.append(proc)
                        proc.start()

                    # complete the processes
                    for proc in procs:
                        proc.join()
            else:
                #####################################
                ###### UN-PARALLELIZED VERSION ######
                for i, f in enumerate(dataset['files'][cond]):
                    process_file(f, i, c, PROTOCOL, quantity)
                #####################################

            # now that we have stored all datafile outputs
            Tunings = []
            for i, f in enumerate(dataset['files'][cond]):

                if os.path.isfile(os.path.join(folder, 'temp', 
                                                'Tuning-%s-%i.npy' % (c, i))):
                    Tuning = np.load(os.path.join(folder, 'temp', 
                                                'Tuning-%s-%i.npy' % (c, i)),
                                        allow_pickle=True).item()
                    Tunings.append(Tuning)

            # # saving data
            np.save(os.path.join(folder, 'Tunings_%s_%s.npy' % (c, quantity)), Tunings)

        else:
            print()
            print('   [!!]   DATASET NOT LARGE ENOUGH   [!!] ')
            print('               only N=%i sessions available' %\
                                        len(dataset['files'][cond]))
            print('   [!!]   DATASET not analyzed       [!!] ')
            print()

        print('-----------------------------------------------------------------')
        print('=================================================================')
# shutil.rmtree(os.path.join(folder, 'temp'))

# %%
from physion.analysis.protocols.orientation_tuning\
    import plot_orientation_tuning_curve, plot_selectivity

for PROTOCOL in PROTOCOLS:
    
    fig, ax = plot_orientation_tuning_curve(\
                            ['shRNA_%s_%s' % (PROTOCOL, quantity),
                             'scramble_%s_%s' % (PROTOCOL, quantity)],
                            average_by='sessions',
                            path=folder)
    fig.suptitle(PROTOCOL)

    fig, ax = plot_orientation_tuning_curve(\
                            ['shRNA_%s_%s' % (PROTOCOL, quantity),
                             'scramble_%s_%s' % (PROTOCOL, quantity)],
                            average_by='ROIs',
                            path=folder)
    fig.suptitle(PROTOCOL)

    fig, ax = plot_selectivity(\
                            ['shRNA_%s_%s' % (PROTOCOL, quantity),
                             'scramble_%s_%s' % (PROTOCOL, quantity)],
                            #   average_by='ROIs',
                            #  using='fit',
                            path=folder)
    fig.suptitle(PROTOCOL)

# %%

