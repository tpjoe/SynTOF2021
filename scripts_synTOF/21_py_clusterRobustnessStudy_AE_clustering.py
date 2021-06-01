"""
This script is similar to script #3, but with the purpose of doing leave one out clustering to 
investigate the robustness of the clusters.
"""


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

disabling_blas()


def pretrain(x_train, identifier, dims, i):
    # set reproducibility
    import tensorflow as tf
    from glob import glob
    import csv
    import numpy as np
    import random
    import pandas as pd
    from tensorflow.keras.callbacks import EarlyStopping
    from tensorflow.keras.initializers import glorot_normal, glorot_uniform, he_normal, lecun_normal
    from tensorflow.keras.layers import concatenate
    from tensorflow.keras.models import Model
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)
    from importlib import reload
    disabling_blas()
    seed_value = 42*i
    cb = EarlyStopping(monitor='r_square', min_delta=0.0025, patience=1, \
        verbose=0, mode='max', baseline=None, restore_best_weights=False)   
    from utils_test import clustering2, clustering3, clustering3_u
    import utils_test
    utils_test.reproducibility(seed_value)
    init = glorot_normal(seed=seed_value)
    dims_a = [x_train.shape[-1]] + dims[0]
    dims_b = [x_train.shape[-1]] + dims[1]
    ae1 = utils_test.autoencoder_(dims_a, init=init, uniqueID='0')
    ae3 = utils_test.autoencoder_(dims_b, uniqueID='2', init=init) 
    opt = tf.keras.optimizers.Adagrad(learning_rate=0.1)
    ae1.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae3.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae1.fit(x=x_train, y=x_train, batch_size=2**10, epochs=20000, callbacks=[cb], shuffle=True)
    ae3.fit(x=x_train, y=x_train, batch_size=2**10, epochs=20000, callbacks=[cb], shuffle=True)
    save_dir = '../results_ae'
    ae1.save_weights(save_dir + '/ae1_' + identifier + '_' + str(i) + '.h5')
    ae3.save_weights(save_dir + '/ae3_' + identifier + '_' + str(i) + '.h5')



def predict(identifier, x_train, i):
    # set reproducibility
    import tensorflow as tf
    import csv
    import numpy as np
    import random
    from tensorflow.keras.layers import concatenate
    from tensorflow.keras.models import Model
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)
    import utils_test
    from utils_test import clustering2, clustering3, ClusteringLayer, clustering2K
    disabling_blas()
    save_dir = '../results_ae/'
    megaAE = tf.keras.models.load_model(save_dir + '/megaAE152_' + identifier + '_' + str(i) + '.h5', 
             custom_objects={'ClusteringLayer': ClusteringLayer})
    q, _, _ = megaAE.predict([x_train, x_train])
    cl_pred = q.argmax(1)
    return cl_pred #, sample



def fit_predict(x_train, identifier, dims, n_clusters_list, i):
    # set reproducibility
    import tensorflow as tf
    import csv
    import numpy as np
    import random
    import pandas as pd
    import flowkit as fk
    from glob import glob
    from rpy2 import robjects
    import rpy2.robjects.packages as rpackages
    from tensorflow.keras.callbacks import EarlyStopping
    from tensorflow.keras.initializers import glorot_normal, glorot_uniform, he_normal, lecun_normal
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
    from utils_test import clustering2, ClusteringLayer, get_cluster_num, clustering2K
    import utils_test
    # for reproducibility
    disabling_blas()
    # seed_value = 42*i
    # utils_test.reproducibility(seed_value)
    dims_a = [x_train.shape[-1]] + dims[0]
    dims_b = [x_train.shape[-1]] + dims[1]
    save_dir = '../results_ae'
    ae1_wt = glob(save_dir + '/ae1_' + identifier +'_*')[i]
    ae3_wt = glob(save_dir + '/ae3_' + identifier +'_*')[i]
    n_clusters = n_clusters_list[i]
    ae1 = utils_test.autoencoder_(dims_a, uniqueID='0')
    ae3 = utils_test.autoencoder_(dims_b, uniqueID='2')
    # opt = tf.keras.optimizers.Nadam(learning_rate=1E-5)
    # opt = tf.keras.optimizers.SGD(learning_rate=0.1, nesterov=True)
    opt = tf.keras.optimizers.Adagrad(learning_rate=0.1)
    ae1.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae3.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae1.load_weights(ae1_wt)
    ae3.load_weights(ae3_wt)
    merged_hidden = concatenate([ae1.get_layer(name='encoder_' + '03').output,
                                 ae3.get_layer(name='encoder_' + '23').output])
    encoder = Model(inputs=[ae1.input, ae3.input], outputs=merged_hidden)
    clustering_layer = utils_test.ClusteringLayer(n_clusters, name='clustering')(merged_hidden)
    megaAE = Model(inputs=[ae1.input, ae3.input],
                   outputs=[clustering_layer, ae1.output, ae3.output])
    megaAE.compile(loss={'clustering': 'kld', 
                        'decoder_' + '00': 'mse',
                        'decoder_' + '20': 'mse'},
                   loss_weights=[0.5, 1/4, 1/4],
                   optimizer='Adam')#tf.keras.optimizers.Adam(learning_rate=0.0001))
    cl = clustering2K(model=megaAE, encoder=encoder, x=x_train, n_clusters=n_clusters, tol=0.03, batch_size=2**10, update_interval=140*5)
    megaAE.save(save_dir + '/megaAE152_' + identifier + '_' + str(i) + '.h5')
    return cl #, sample.to_list()
#best rn batch^10, intervalx1, tol0.03

def automated_cluster(x_train, identifier, dims):
    # set reproducibility
    import tensorflow as tf
    import csv
    import numpy as np
    import random
    import pandas as pd
    import flowkit as fk
    from glob import glob
    from rpy2 import robjects
    import rpy2.robjects.packages as rpackages
    from tensorflow.keras.callbacks import EarlyStopping
    from tensorflow.keras.initializers import glorot_normal, glorot_uniform, he_normal, lecun_normal
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
    from utils_test import clustering2, clustering3, get_cluster_num
    # for reproducibility
    disabling_blas()
    n_clusters = []
    x_train = x_train.to_numpy()
    dims_a = [x_train.shape[-1]] + dims[0]
    dims_b = [x_train.shape[-1]] + dims[1]
    save_dir = '../results_ae'
    for i in range(len(glob(save_dir + '/ae1_' + identifier +'_*'))):
        print('Working on best cluster number for rep {}'.format(i))
        ae1_wt = glob(save_dir + '/ae1_' + identifier +'_*')[i]
        ae3_wt = glob(save_dir + '/ae3_' + identifier +'_*')[i]
        ae1 = utils_test.autoencoder_(dims_a, uniqueID='0')
        ae3 = utils_test.autoencoder_(dims_b, uniqueID='2')
        opt = tf.keras.optimizers.Adagrad(learning_rate=0.1)
        ae1.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
        ae3.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
        ae1.load_weights(ae1_wt)
        ae3.load_weights(ae3_wt)
        merged_hidden = concatenate([ae1.get_layer(name='encoder_' + '03').output,
                                    ae3.get_layer(name='encoder_' + '23').output])
        encoder = Model(inputs=[ae1.input, ae3.input], outputs=merged_hidden)
        # get concat layer name
        layer_names = [layer.name for layer in encoder.layers]
        concat_ind = np.max(np.where(['concatenate' in layer for layer in layer_names]))
        concat_layer = layer_names[concat_ind]
        get_output = K.function([encoder.input], [encoder.get_layer(concat_layer).output])
        h = get_output([x_train, x_train, x_train])[0]
        n_clusters.append(int(get_cluster_num(h, maxK=40, subsampling_n=50000, 
                           plot_dir='rss_plots/distortions_' + identifier + '_rep' + str(i) + '.png')))
    return n_clusters



import os
import numpy as np
import pandas as pd
import flowkit as fk
from glob import glob
import multiprocessing
from joblib import Parallel, delayed
from sklearn.preprocessing import StandardScaler, QuantileTransformer
from itertools import combinations 
from sklearn.preprocessing import StandardScaler, MinMaxScaler, normalize, QuantileTransformer


num_cores = multiprocessing.cpu_count()
reps = 10
sess = 1
fcs_path = '../raw_data/max_events/fcs/'
files = np.sort(glob(fcs_path + '*_LowNo*.fcs'))
files = np.array([x for x in files if (('HF14-017' not in x) & ('HF14-083' not in x) & ('HF14-025' not in x))])
# omit 2 samples for testing
file_options = ['HF13-117', 'HF14-008', 'HF14-051', 'HF14-053', 'HF14-057', 'HF14-076']
pairs = list(combinations(file_options, 1))


# for p in pairs:
for p in pairs:
    pair = p[0]
    files_train = np.array([x for x in files if (pair not in x)]) # & (pair[1] not in x))])
    identifier = 'allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_' + ','.join(pair)
    # identifier = 'real10' #'allLowNoPresynaptic_105_SGDwithVal_lr_batch210'
    dims = [[512, 256, 128, 10], [512, 256, 128, 5]]
    # load files
    fcs_list = []
    sample = []

    for file in files_train:
        ff = fk.Sample(file)
        events = ff.get_orig_events()
        sample = sample + ([file.split('_')[-1]] * events.shape[0])
        fcs_list.append(events)

    df = pd.DataFrame(np.vstack(fcs_list), columns=ff.pnn_labels)
    sample = pd.Series(sample)
    excludedPro = ['b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin']
                #    'PARKIN', 'TMEM230_C20orf30', 'DJ-1_PARK7', 'GBA1'] #possible
    df = df.loc[:, ~df.columns.isin(excludedPro)]
    df = df.drop(['NET'], axis=1)
    ## %% --------------------------------------------------------------------------------
    x_train = np.array(df)#.sample(frac=0.02, random_state=0))

    n_clusters_list = [15]*10
    res_ = Parallel(n_jobs=reps)(delayed(fit_predict)(pd.DataFrame(x_train), identifier, dims, n_clusters_list, i) for i in range(reps))

    # predict and export to R ---------------------------------------------------------------------
    def get_predict(files, identifier_pred, reps, post=False):
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
        excludedPro = ['b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                        'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin']
                    #    'PARKIN', 'TMEM230_C20orf30', 'DJ-1_PARK7', 'GBA1'] #possible
        df = df.loc[:, ~df.columns.isin(excludedPro)]
        df = df.drop(['NET'], axis=1)
        ## %% --------------------------------------------------------------------------------
        x_train = np.array(df)#.sample(frac=0.02, random_state=0).reset_index(drop=True))
        res = Parallel(n_jobs=reps)(delayed(predict)(identifier, pd.DataFrame(x_train), i) for i in range(reps))
        cl_pred = [res[i] for i in range(len(res))]
        cl_pred = pd.DataFrame(np.column_stack(cl_pred))
        to_R = pd.concat([cl_pred, pd.DataFrame(sample_pred).reset_index(drop=True).rename(columns={0:'sample'})], axis=1)
        return to_R

    # get prediction of presynaptic in different groups
    fcs_path = '../raw_data/max_events/fcs/'

    identifier_pred = 'predLowNo' + '_maxK40_' + identifier
    to_R = get_predict(files_train, identifier_pred, reps)
    to_R.to_csv('R_py_exchange/presynTOF_AdamMegaAE152' + identifier_pred + '_sess_' + str(sess) + '_no_' + ','.join(pair) + '.csv')

    files_test = np.array([x for x in files if ((pair in x))])
    identifier_pred = 'predLowNo' + '_maxK40_' + identifier
    to_R = get_predict(files_test, identifier_pred, reps)
    to_R.to_csv('R_py_exchange/presynTOF_AdamMegaAE152' + identifier_pred.replace(',', '') + '_sess_' + str(sess) + '_for_' + pair + '.csv')
