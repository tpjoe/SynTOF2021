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
import os
os.chdir('/home/tpjoe/tpjoe@stanford.edu/project_SynTOF/scripts_synTOF/')
import sys
sys.path.insert(1, '/home/tpjoe/tpjoe@stanford.edu/project_SynTOF/scripts_synTOF/')


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
    from keras_adabound import AdaBound
    # for reproducibility
    disabling_blas()
    seed_value = 42*i
    # cb = EarlyStopping(monitor='loss', min_delta=0, patience=1, \
        # verbose=0, mode='auto', baseline=None, restore_best_weights=False) #test
    # cb = EarlyStopping(monitor='r_square', min_delta=0.0005, patience=1, \
        # verbose=0, mode='max', baseline=None, restore_best_weights=False)    #test2, real(ada) 0.57, real2(sgd) 0.62
    cb = EarlyStopping(monitor='r_square', min_delta=0.0025, patience=1, \
        verbose=0, mode='max', baseline=None, restore_best_weights=False)   #real4(ada) 0.71, real5(sgd) 0.62, real8(ada, n=6) , real9(lr=0.01) 0.66, real10(adaBound)
    # cb = EarlyStopping(monitor='r_square', min_delta=0.01, patience=1, \
        # verbose=0, mode='max', baseline=None, restore_best_weights=False)    #real6(ada) 0.59 , rea7(sgd) 0.69
    from utils_test import clustering2, clustering3, clustering3_u
    import utils_test
    utils_test.reproducibility(seed_value)
    init = glorot_normal(seed=seed_value)
    dims_a = [x_train.shape[-1]] + dims[0]
    dims_b = [x_train.shape[-1]] + dims[1]
    ae1 = utils_test.autoencoder_(dims_a, init=init, uniqueID='0')
    ae3 = utils_test.autoencoder_(dims_b, uniqueID='2', init=init)
    # opt = tf.keras.optimizers.Nadam(learning_rate=1E-5)
    # lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
    #     initial_learning_rate=0.1,
    #     decay_steps=1000,
    #     decay_rate=0.85)
    # opt = AdaBound(lr=1e-3, final_lr=0.1)#tf.keras.optimizers.SGD(nesterov=True, learning_rate=0.1, momentum=0.9) #learning_rate=0.1, 
    opt = tf.keras.optimizers.Adagrad(learning_rate=0.1)
    # opt = AdaBound(lr=1e-02, final_lr=0.1)
    ae1.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae3.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae1.fit(x=x_train, y=x_train, batch_size=2**15, epochs=20000, callbacks=[cb], shuffle=True)
    ae3.fit(x=x_train, y=x_train, batch_size=2**15, epochs=20000, callbacks=[cb], shuffle=True)
    # ae1.compile(optimizer='Nadam', loss='mse', metrics=[utils_test.r_square])
    # ae3.compile(optimizer='Nadam', loss='mse', metrics=[utils_test.r_square])
    # ae1.fit(x=x_train, y=x_train, batch_size=2**10, epochs=1, callbacks=[cb], shuffle=True)
    # ae3.fit(x=x_train, y=x_train, batch_size=2**10, epochs=1, callbacks=[cb], shuffle=True)
    save_dir = '../results_ae'
    ae1.save_weights(save_dir + '/ae1_' + identifier + '_' + str(i) + '.h5')
    ae3.save_weights(save_dir + '/ae3_' + identifier + '_' + str(i) + '.h5')


def pretrain_withval(x_train, identifier, dims, i):
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
    from keras_adabound import AdaBound
    # for reproducibility
    disabling_blas()
    seed_value = 42*i
    # cb = EarlyStopping(monitor='val_loss', min_delta=0, patience=3, \
        # verbose=0, mode='auto', baseline=None, restore_best_weights=False) #test
    # cb = EarlyStopping(monitor='r_square', min_delta=0.0005, patience=0, \
        # verbose=0, mode='max', baseline=None, restore_best_weights=False)    #test2 
    lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
        initial_learning_rate=0.1,
        decay_steps=2000*50,
        decay_rate=0.90)
    from utils_test import clustering2, clustering3, clustering3_u
    import utils_test
    x_trainVal = x_train.sample(frac=0.1, random_state=i)
    x_train = x_train.drop(x_trainVal.index, axis=0)
    utils_test.reproducibility(seed_value)
    init = glorot_normal(seed=seed_value)
    dims_a = [x_train.shape[-1]] + dims[0]
    dims_b = [x_train.shape[-1]] + dims[1]
    ae1 = utils_test.autoencoder_(dims_a, init=init, uniqueID='0')
    ae3 = utils_test.autoencoder_(dims_b, uniqueID='2', init=init)
    opt = tf.keras.optimizers.Adagrad(learning_rate=lr_schedule)
    ae1.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae3.compile(optimizer=opt, loss='mse', metrics=[utils_test.r_square])
    ae1.fit(x=x_train, y=x_train, batch_size=2**10, epochs=20000, callbacks=[cb], validation_data=(x_trainVal, x_trainVal), shuffle=True)
    ae3.fit(x=x_train, y=x_train, batch_size=2**10, epochs=20000, callbacks=[cb], validation_data=(x_trainVal, x_trainVal), shuffle=True)
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
    megaAE = tf.keras.models.load_model(save_dir + '/megaAE_' + identifier + '_' + str(i) + '.h5', 
             custom_objects={'ClusteringLayer': ClusteringLayer})
    q, _, _ = megaAE.predict([x_train, x_train])
    cl_pred = q.argmax(1)
    return cl_pred #, sample


# import tensorflow as tf
# import csv
# import numpy as np
# import random
# import pandas as pd
# import flowkit as fk
# from glob import glob
# from rpy2 import robjects
# import rpy2.robjects.packages as rpackages
# from tensorflow.keras.callbacks import EarlyStopping
# from tensorflow.keras.initializers import glorot_normal, glorot_uniform, he_normal, lecun_normal
# from tensorflow.keras.layers import concatenate
# from tensorflow.keras.models import Model
# from tensorflow.keras import backend as K
# gpus = tf.config.experimental.list_physical_devices('GPU')
# if gpus:
#     try:
#         for gpu in gpus:
#             tf.config.experimental.set_memory_growth(gpu, True)
#     except RuntimeError as e:
#         print(e)
# # import other scripts
# from importlib import reload
# from utils_test import clustering2, ClusteringLayer, get_cluster_num, clustering2K
# import utils_test
# y_pred = []
# x = np.array(df.sample(frac=0.02, random_state=0))
# for i in range(6):
#     dims_a = [x_train.shape[-1]] + dims[0]
#     dims_b = [x_train.shape[-1]] + dims[1]
#     save_dir = '../results_ae'
#     ae1_wt = glob(save_dir + '/ae1_' + identifier +'_*')[i]
#     ae3_wt = glob(save_dir + '/ae3_' + identifier +'_*')[i]
#     n_clusters = 16
#     ae1 = utils_test.autoencoder_(dims_a, uniqueID='0')
#     ae3 = utils_test.autoencoder_(dims_b, uniqueID='2')
#     ae1.compile(optimizer='Nadam', loss='mse', metrics=[utils_test.r_square])
#     ae3.compile(optimizer='Nadam', loss='mse', metrics=[utils_test.r_square])
#     ae1.load_weights(ae1_wt)
#     ae3.load_weights(ae3_wt)
#     merged_hidden = concatenate([ae1.get_layer(name='encoder_' + '03').output,
#                                     ae3.get_layer(name='encoder_' + '23').output])
#     # encoder = Model(inputs=ae1.input, outputs=ae1.output)
#     encoder = Model(inputs=[ae1.input, ae3.input], outputs=merged_hidden)
#     from sklearn.cluster import KMeans, MiniBatchKMeans
#     k_seed = None
#     kmeans = KMeans(n_clusters=n_clusters, random_state=k_seed, n_init=5)
#     y_pred += [kmeans.fit_predict(encoder.predict([x, x]))]


# from sklearn.metrics import normalized_mutual_info_score as NMI
# NMI(y_pred[3], y_pred[5])
# # pd.concat([pd.Series(y_pred), pd.Series(y_pred1)], axis=1).head(30)


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
                   optimizer=opt)
    cl = clustering2K(model=megaAE, encoder=encoder, x=x_train, n_clusters=n_clusters, tol=0.03, batch_size=2**10, update_interval=140*5)
    megaAE.save(save_dir + '/megaAE_' + identifier + '_' + str(i) + '.h5')
    return cl #, sample.to_list()


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
        ae1.compile(optimizer='Nadam', loss='mse', metrics=[utils_test.r_square])
        ae3.compile(optimizer='Nadam', loss='mse', metrics=[utils_test.r_square])
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


def get_hidden(x_train, identifier, i):
    # set reproducibility
    import tensorflow as tf
    import csv
    import numpy as np
    import random
    import pandas as pd
    import flowkit as fk
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
    
    from utils_test import clustering2, clustering3, ClusteringLayer, clustering2K
    import utils_test
    disabling_blas()
    save_dir = '../results_ae/'
    megaAE = tf.keras.models.load_model(save_dir + '/megaAE_' + identifier + '_' + str(i) + '.h5', 
             custom_objects={'ClusteringLayer': ClusteringLayer})
    layer_names = [layer.name for layer in megaAE.layers]
    concat_ind = np.max(np.where(['concatenate' in layer for layer in layer_names]))
    concat_layer = layer_names[concat_ind]
    encoder = Model(inputs=megaAE.input, outputs=megaAE.get_layer(name=concat_layer).output)
    hidden = encoder.predict([x_train, x_train])
    return hidden #, sample


import os
import numpy as np
import pandas as pd
import flowkit as fk
from glob import glob
import multiprocessing
from joblib import Parallel, delayed
from sklearn.preprocessing import StandardScaler, QuantileTransformer



num_cores = multiprocessing.cpu_count()
reps = 10
sess = 1
fcs_path = '../raw_data/max_events/fcs/'
files = np.sort(glob(fcs_path + '*_LowNo*.fcs'))
identifier = 'allLowNoPresynaptic_LowNo_08062020Batch210_105_Adagradlr01_noStd'
# identifier = 'real10' #'allLowNoPresynaptic_105_SGDwithVal_lr_batch210'
dims = [[512, 256, 128, 10], [512, 256, 128, 5]]


# load files
fcs_list = []
sample = []
for file in files:
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
# x_train = np.vstack([x_train, x_train_post])

# pretrain_withval(pd.DataFrame(x_train), identifier, dims, i)
 # adagrad0.1 0.0222/0.0236, 0.0424/0.0428 
 # sgd0.1 0.0212/0.0210, 0.0460/0.0475
 # adagrad
 # sgd 0.0196/0.0197

# automatic runs
Parallel(n_jobs=num_cores)(delayed(pretrain)(pd.DataFrame(x_train), identifier, dims, i) for i in range(reps))
n_clusters_list = automated_cluster(pd.DataFrame(x_train), identifier, dims)
# n_clusters_list = [16]*5
res_ = Parallel(n_jobs=reps)(delayed(fit_predict)(pd.DataFrame(x_train), identifier, dims, n_clusters_list, i) for i in range(reps))



# predict and export to R ---------------------------------------------------------------------
import pandas as pd
import flowkit as fk
from sklearn.preprocessing import StandardScaler, MinMaxScaler, normalize, QuantileTransformer


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
    # scaler = StandardScaler(with_mean=False)
    # x_train = scaler.fit_transform(x_train)
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
# scaler = StandardScaler(with_mean=False)
# x_train = scaler.fit_transform(x_train)

res = Parallel(n_jobs=reps)(delayed(get_hidden)(pd.DataFrame(x_train), identifier, i) for i in range(reps))
hidden = [res[i] for i in range(len(res))]
hidden_ = pd.DataFrame(np.column_stack(hidden))

to_R = pd.concat([hidden_, pd.DataFrame(sample_pred).rename(columns={0:'sample'})], axis=1)
to_R.to_csv('R_py_exchange/hidden_' + identifier + '_sess_' + str(sess) + '.csv')
