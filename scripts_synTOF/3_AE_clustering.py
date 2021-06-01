##%% libraries --------------------------------------------------------------------------------
# to get reproducible results
def disabling_blas():
    n_threads = 1
    import os
    os.environ['OMP_NUM_THREADS'] = str(n_threads)
    os.environ['OPENBLAS_NUM_THREADS'] = str(n_threads)
    os.environ['MKL_NUM_THREADS'] = str(n_threads)
    os.environ['VECLIB_MAXIMUM_THREADS'] = str(n_threads)
    os.environ['NUMEXPR_NUM_THREADS'] = str(n_threads)
    # tf.config.threading.set_inter_op_parallelism_threads(1)

disabling_blas() # to run purely on cpu


def pretrain(x_train, identifier, dims, i):
    # this function pretrains the AEs before attaching clustering layer to it, identifier and i are used for generating model name tag.
    # dims are used for node sizes in each layer
    # set reproducibility
    import tensorflow as tf
    from tensorflow.keras.callbacks import EarlyStopping
    from tensorflow.keras.initializers import glorot_normal
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)
    import utils_test
    # for reproducibility
    disabling_blas()
    seed_value = 42*i
    cb = EarlyStopping(monitor='r_square', min_delta=0.0025, patience=1, \
        verbose=0, mode='max', baseline=None, restore_best_weights=False)  
    utils_test.reproducibility(seed_value)
    # initialize model with specified dimensions (dims)
    init = glorot_normal(seed=seed_value)
    dims_a = [x_train.shape[-1]] + dims[0]
    dims_b = [x_train.shape[-1]] + dims[1]
    ae1 = utils_test.autoencoder_(dims_a, init=init, uniqueID='0')
    ae3 = utils_test.autoencoder_(dims_b, uniqueID='2', init=init)
    opt = tf.keras.optimizers.Adagrad(learning_rate=0.1)
    # compile models
    ae1.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae3.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    # start training
    ae1.fit(x=x_train, y=x_train, batch_size=2**10, epochs=20000, callbacks=[cb], shuffle=True)
    ae3.fit(x=x_train, y=x_train, batch_size=2**10, epochs=20000, callbacks=[cb], shuffle=True)
    save_dir = '../results_ae'
    # save models
    ae1.save_weights(save_dir + '/ae1_' + identifier + '_' + str(i) + '.h5')
    ae3.save_weights(save_dir + '/ae3_' + identifier + '_' + str(i) + '.h5')


def fit_megaAE(x_train, identifier, dims, n_clusters_list, i):
    # This function combine the 2 AEs together, attach the clustering layer, and train the model for clustering
    # set reproducibility
    import tensorflow as tf
    from glob import glob
    from tensorflow.keras.layers import concatenate
    from tensorflow.keras.models import Model
    from tensorflow.keras import backend as K
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)
    # import other scripts
    from utils_test import clustering2K
    import utils_test
    # for reproducibility
    disabling_blas()
    # build the exact same AE models
    dims_a = [x_train.shape[-1]] + dims[0]
    dims_b = [x_train.shape[-1]] + dims[1]
    save_dir = '../results_ae'
    ae1_wt = glob(save_dir + '/ae1_' + identifier +'_*')[i]
    ae3_wt = glob(save_dir + '/ae3_' + identifier +'_*')[i]
    n_clusters = n_clusters_list[i]
    ae1 = utils_test.autoencoder_(dims_a, uniqueID='0')
    ae3 = utils_test.autoencoder_(dims_b, uniqueID='2')
    opt = tf.keras.optimizers.Adagrad(learning_rate=0.1)
    ae1.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae3.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    # initialize weights using the pretrained results
    ae1.load_weights(ae1_wt)
    ae3.load_weights(ae3_wt)
    # build mega AE
    merged_hidden = concatenate([ae1.get_layer(name='encoder_' + '03').output,
                                 ae3.get_layer(name='encoder_' + '23').output])
    encoder = Model(inputs=[ae1.input, ae3.input], outputs=merged_hidden)
    # add clustering layer
    clustering_layer = utils_test.ClusteringLayer(n_clusters, name='clustering')(merged_hidden)
    megaAE = Model(inputs=[ae1.input, ae3.input],
                   outputs=[clustering_layer, ae1.output, ae3.output])
    megaAE.compile(loss={'clustering': 'kld', 
                        'decoder_' + '00': 'mse',
                        'decoder_' + '20': 'mse'},
                   loss_weights=[0.5, 1/4, 1/4],
                   optimizer='Adam')
    # run clustering
    cl = clustering2K(model=megaAE, encoder=encoder, x=x_train, n_clusters=n_clusters, tol=0.03, batch_size=2**10, update_interval=140*5)
    megaAE.save(save_dir + '/megaAE_' + identifier + '_' + str(i) + '.h5')
    return cl #, sample.to_list()


def automated_cluster(x_train, identifier, dims):
    # This function outputs the optimal number of clusters for x_train. Identifier and i are used for finding the trained model name
    # set reproducibility
    import tensorflow as tf
    import numpy as np
    from glob import glob
    from tensorflow.keras.layers import concatenate
    from tensorflow.keras.models import Model
    from tensorflow.keras import backend as K
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)
    # import other scripts
    from importlib import reload
    import utils_test
    reload(utils_test)
    from utils_test import get_cluster_num
    # for reproducibility
    disabling_blas()
    n_clusters = []
    x_train = x_train.to_numpy()
    dims_a = [x_train.shape[-1]] + dims[0]
    dims_b = [x_train.shape[-1]] + dims[1]
    save_dir = '../results_ae'
    for i in range(len(glob(save_dir + '/ae1_' + identifier +'_*'))):
        print('Working on best cluster number for rep {}'.format(i))
        # load saved model
        ae1_wt = glob(save_dir + '/ae1_' + identifier +'_*')[i]
        ae3_wt = glob(save_dir + '/ae3_' + identifier +'_*')[i]
        ae1 = utils_test.autoencoder_(dims_a, uniqueID='0')
        ae3 = utils_test.autoencoder_(dims_b, uniqueID='2')
        opt = tf.keras.optimizers.Adagrad(learning_rate=0.1)
        ae1.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
        ae3.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
        ae1.load_weights(ae1_wt)
        ae3.load_weights(ae3_wt)
        # define layers
        merged_hidden = concatenate([ae1.get_layer(name='encoder_' + '03').output,
                                    ae3.get_layer(name='encoder_' + '23').output])
        encoder = Model(inputs=[ae1.input, ae3.input], outputs=merged_hidden)
        # get concat layer name
        layer_names = [layer.name for layer in encoder.layers]
        concat_ind = np.max(np.where(['concatenate' in layer for layer in layer_names]))
        concat_layer = layer_names[concat_ind]
        get_output = K.function([encoder.input], [encoder.get_layer(concat_layer).output])
        # get hidden rep. from each AE
        h = get_output([x_train, x_train, x_train])[0]
        # get cluster number from the concatenated hidden rep.
        n_clusters.append(int(get_cluster_num(h, maxK=40, subsampling_n=50000, 
                           plot_dir='rss_plots/distortions_' + identifier + '_rep' + str(i) + '.png')))
    return n_clusters


def predict(identifier, x_train, i):
    # This function outputs predicted clusters for x_train. Identifier and i are used for finding the trained model name
    # set reproducibility
    import tensorflow as tf
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)
    from utils_test import ClusteringLayer
    disabling_blas()
    # load saved model
    save_dir = '../results_ae/'
    megaAE = tf.keras.models.load_model(save_dir + '/megaAE_' + identifier + '_' + str(i) + '.h5', 
             custom_objects={'ClusteringLayer': ClusteringLayer})
    # predict outputs
    q, _, _ = megaAE.predict([x_train, x_train])
    cl_pred = q.argmax(1)
    return cl_pred #, sample



def get_hidden(x_train, identifier, i):
    # This function outputs the hidden representation for x_train. Identifier and i are used for finding the trained model name
    # set reproducibility
    import tensorflow as tf
    import numpy as np
    from tensorflow.keras.models import Model
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)
    from utils_test import ClusteringLayer
    import utils_test
    disabling_blas()
    # load saved model
    save_dir = '../results_ae/'
    megaAE = tf.keras.models.load_model(save_dir + '/megaAE_' + identifier + '_' + str(i) + '.h5', 
             custom_objects={'ClusteringLayer': ClusteringLayer})
    # get name of each layer
    layer_names = [layer.name for layer in megaAE.layers]
    concat_ind = np.max(np.where(['concatenate' in layer for layer in layer_names]))
    concat_layer = layer_names[concat_ind]
    # use name to define output section of the model
    encoder = Model(inputs=megaAE.input, outputs=megaAE.get_layer(name=concat_layer).output)
    hidden = encoder.predict([x_train, x_train])
    return hidden #, sample











###################################
# Main
###################################

# import libs
import os
import numpy as np
import pandas as pd
import flowkit as fk
from glob import glob
import multiprocessing
from joblib import Parallel, delayed


# define running parameters
num_cores = multiprocessing.cpu_count() # use all cores
reps = 10 # run 10 times
sess = 1 # just a naming of session
fcs_path = '../raw_data/max_events/fcs/' # path to your syntof fcs files
files = np.sort(glob(fcs_path + '*_LowNo*.fcs')) # grab desired files
files = np.array([x for x in files if (('HF14-017' not in x) & ('HF14-083' not in x) & ('HF14-025' not in x))]) # remove non-true LowNo files
identifier = 'allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd' # naming model identifier
dims = [[512, 256, 128, 10], [512, 256, 128, 5]] # node size in each layer of AE1 and AE3


# load pre-synaptic fcs files
fcs_list = []
sample = []
for file in files:
    ff = fk.Sample(file)
    events = ff.get_orig_events()
    sample = sample + ([file.split('_')[-1]] * events.shape[0])
    fcs_list.append(events)

df = pd.DataFrame(np.vstack(fcs_list), columns=ff.pnn_labels)
sample = pd.Series(sample)
# excluding non-phenotypic markers
excludedPro = ['b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin']
df = df.loc[:, ~df.columns.isin(excludedPro)]
df = df.drop(['NET'], axis=1) # drop NET because of low quality


# # load post-synaptic fcs files
# fcs_path_post = '../raw_data/max_events/fcs_post_synap/'
# files_post = np.sort(glob(fcs_path_post + '*_LowNo*.fcs'))
# # load post-files
# fcs_list = []
# sample = []
# for file in files_post:
#     ff = fk.Sample(file)
#     events = ff.get_orig_events()
#     sample = sample + ([file.split('_')[-1]] * events.shape[0])
#     fcs_list.append(events)

# df_post = pd.DataFrame(np.vstack(fcs_list), columns=ff.pnn_labels)
# sample = pd.Series(sample)
# excludedPro = ['b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
#                 'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin']
#             #    'PARKIN', 'TMEM230_C20orf30', 'DJ-1_PARK7', 'GBA1'] #possible
# df_post = df_post.loc[:, ~df_post.columns.isin(excludedPro)]
# df_post = df_post.drop(['NET'], axis=1)


## %% --------------------------------------------------------------------------------
x_train = np.array(df)#.sample(frac=0.02, random_state=0))
# x_train_post = np.array(df_post)


# run pretrain (10x in parallel)
Parallel(n_jobs=num_cores)(delayed(pretrain)(pd.DataFrame(x_train), identifier, dims, i) for i in range(reps))
# run getting optimal cluster numbers
n_clusters_list = automated_cluster(pd.DataFrame(x_train), identifier, dims)
# run clustering in parallel
res_ = Parallel(n_jobs=reps)(delayed(fit_megaAE)(pd.DataFrame(x_train), identifier, dims, n_clusters_list, i) for i in range(reps))








# predict clusters and export to csv ---------------------------------------------------------------------
import pandas as pd
import flowkit as fk


def get_predict(files, identifier, reps, post=False):
    # This function loads wanted data and call prediction function of the model with tag "identifier". 
    # It outputs cluster prediction for each event (along with the event itself)
    excludedPro = ['b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                   'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin']
    # load files
    fcs_list = []
    sample_pred = []
    for file in files:
        ff = fk.Sample(file)
        events = ff.get_orig_events()
        if post:
            sample_pred = sample_pred + (['_'.join(['post', file.split('/')[4].split('_')[0], file.split('_')[5], 
                                                    file.split('_')[6], file.split('_')[-1]])] * events.shape[0])
        else:
            sample_pred = sample_pred + (['_'.join(['pre', file.split('/')[4].split('_')[0], file.split('_')[3], 
                                                    file.split('_')[4], file.split('_')[-1]])] * events.shape[0])
        fcs_list.append(events)
    df = pd.DataFrame(np.vstack(fcs_list), columns=ff.pnn_labels)
    sample_pred = pd.Series(sample_pred)#.sample(frac=0.02, random_state=0)
    df = df.loc[:, ~df.columns.isin(excludedPro)]
    df = df.drop(['NET'], axis=1)
    # get predictions
    x_train = np.array(df)#.sample(frac=0.02, random_state=0).reset_index(drop=True))
    res = Parallel(n_jobs=reps)(delayed(predict)(identifier, pd.DataFrame(x_train), i) for i in range(reps))
    cl_pred = [res[i] for i in range(len(res))]
    cl_pred = pd.DataFrame(np.column_stack(cl_pred))
    to_R = pd.concat([cl_pred, pd.DataFrame(sample_pred).reset_index(drop=True).rename(columns={0:'sample'})], axis=1)
    return to_R


# get prediction of presynaptic in different groups
fcs_path = '../raw_data/max_events/fcs/'

files = np.sort(glob(fcs_path + '*_LowNo*.fcs'))
identifier_pred = 'predLowNo' + '_maxK40_' + identifier
to_R = get_predict(files, identifier_pred, reps)
to_R.to_csv('R_py_exchange/presynTOF_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')


files = np.sort(glob(fcs_path + '*_LBD*.fcs'))
identifier_pred = 'predLBD' + '_maxK40_' + identifier
to_R = get_predict(files, identifier_pred, reps)
to_R.to_csv('R_py_exchange/presynTOF_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')


files = np.sort(glob(fcs_path + '*_PHAD*.fcs'))
identifier_pred = 'predPHAD' + '_maxK40_' + identifier
to_R = get_predict(files, identifier_pred, reps)
to_R.to_csv('R_py_exchange/presynTOF_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')


# get prediction of postsynaptic in different groups
fcs_path = '../raw_data/max_events/fcs_post_synap/'

files = np.sort(glob(fcs_path + '*_LowNo*.fcs'))
identifier_pred = 'predLowNo' + '_maxK40_' + identifier
to_R = get_predict(files, identifier_pred, reps, post=True)
to_R.to_csv('R_py_exchange/postsynTOF_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')


files = np.sort(glob(fcs_path + '*_LBD*.fcs'))
identifier_pred = 'predLBD' + '_maxK40_' + identifier
to_R = get_predict(files, identifier_pred, reps, post=True)
to_R.to_csv('R_py_exchange/postsynTOF_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')


files = np.sort(glob(fcs_path + '*_PHAD*.fcs'))
identifier_pred = 'predPHAD' + '_maxK40_' + identifier
to_R = get_predict(files, identifier_pred, reps, post=True)
to_R.to_csv('R_py_exchange/postsynTOF_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')


# # get prediction of GFAP- EAAT1- presynaptic
# fcs_path = '../raw_data/max_events/fcs_GFAPnegEAAT1neg/'

# files = np.sort(glob(fcs_path + '*_LowNo*.fcs'))
# identifier_pred = 'predLowNo' + '_maxK40_' + identifier
# to_R = get_predict(files, identifier_pred, reps)
# to_R.to_csv('R_py_exchange/presynTOFGFAPnegEAAT1neg_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')


# files = np.sort(glob(fcs_path + '*_LBD*.fcs'))
# identifier_pred = 'predLBD' + '_maxK40_' + identifier
# to_R = get_predict(files, identifier_pred, reps)
# to_R.to_csv('R_py_exchange/presynTOFGFAPnegEAAT1neg_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')


# files = np.sort(glob(fcs_path + '*_PHAD*.fcs'))
# identifier_pred = 'predPHAD' + '_maxK40_' + identifier
# to_R = get_predict(files, identifier_pred, reps)
# to_R.to_csv('R_py_exchange/presynTOFGFAPnegEAAT1neg_AdamMegaAE' + identifier_pred + '_sess_' + str(sess) + '.csv')



# get hidden and export to R -----------------------------------------------------------------
fcs_path = '../raw_data/max_events/fcs/'
excludedPro = ['b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin']

files = np.sort(glob(fcs_path + '*_LowNo*.fcs'))

# load files
fcs_list = []
sample_pred = []
for file in files:
    ff = fk.Sample(file)
    events = ff.get_orig_events()
    sample_pred = sample_pred + (['_'.join([file.split('/')[4].split('_')[0], file.split('_')[3], file.split('_')[-1]])] * events.shape[0])
    fcs_list.append(events)

df = pd.DataFrame(np.vstack(fcs_list), columns=ff.pnn_labels)
sample_pred = pd.Series(sample_pred)#.sample(frac=0.02, random_state=0).reset_index(drop=True)
excludedPro = ['b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin']
            #    'PARKIN', 'TMEM230_C20orf30', 'DJ-1_PARK7', 'GBA1'] #possible
df = df.loc[:, ~df.columns.isin(excludedPro)]
df = df.drop(['NET'], axis=1)

## %% --------------------------------------------------------------------------------
x_train = np.array(df)#.sample(frac=0.02, random_state=0))


res = Parallel(n_jobs=reps)(delayed(get_hidden)(pd.DataFrame(x_train), identifier, i) for i in range(reps))
hidden = [res[i] for i in range(len(res))]
hidden_ = pd.DataFrame(np.column_stack(hidden))

to_R = pd.concat([hidden_, pd.DataFrame(sample_pred).rename(columns={0:'sample'})], axis=1)
to_R.to_csv('R_py_exchange/hidden_' + identifier + '_sess_' + str(sess) + '.csv')