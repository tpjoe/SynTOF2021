# libraries
import pandas as pd
from dask import dataframe as dd
import matplotlib.pyplot as plt

# import data
mibi = dd.read_csv('../raw_data/tifData.csv').compute()

# screen for rough typr of excitory/non-excitory synaposes
#| (mibi.CD56>0.004) (mibi.TotalTau>0.01) (mibi.VGLUT2<0.002)
excite_condition = ((mibi.ApoE>0.004) | (mibi.VGLUT1>0.004) | (mibi.TotalTau>0.001) | (mibi.CD56>0.004))
excite_syn = mibi.loc[excite_condition& ((mibi.Synaptophysin>0)), :]
nonExcite_syn = mibi.loc[~excite_condition & ((mibi.Synaptophysin>0)), :]

excite_syn_ad = excite_syn.loc[excite_syn.Point.isin([1, 2, 3, 4, 5, 6]), :]
excite_syn_hc = excite_syn.loc[excite_syn.Point.isin([7, 8, 9, 10, 11, 12]), :]
nonExcite_syn_ad = nonExcite_syn.loc[nonExcite_syn.Point.isin([1, 2, 3, 4, 5, 6]), :]
nonExcite_syn_hc = nonExcite_syn.loc[nonExcite_syn.Point.isin([7, 8, 9, 10, 11, 12]), :]

# frequency
excite_syn_ad.groupby(['Point']).count().iloc[:, 0]
excite_syn_hc.groupby(['Point']).count().iloc[:, 0]
nonExcite_syn_ad.groupby(['Point']).count().iloc[:, 0]
nonExcite_syn_hc.groupby(['Point']).count().iloc[:, 0]


marker = 'CD47'
marker = 'PHF1Tau'

excite_syn_ad.groupby(['Point']).mean().loc[:, marker].median()
excite_syn_hc.groupby(['Point']).mean().loc[:, marker].median()
nonExcite_syn_ad.groupby(['Point']).mean().loc[:, marker].mean()
nonExcite_syn_hc.groupby(['Point']).mean().loc[:, marker].mean()

excite_syn_ad.groupby(['Point']).mean().loc[:, marker].std()
excite_syn_hc.groupby(['Point']).mean().loc[:, marker].std()
nonExcite_syn_ad.groupby(['Point']).mean().loc[:, marker].std()
nonExcite_syn_hc.groupby(['Point']).mean().loc[:, marker].std()


ad_tau_low = nonExcite_syn_ad.groupby(['Point']).mean().loc[:, 'PHF1Tau']
hc_tau_low = nonExcite_syn_hc.groupby(['Point']).mean().loc[:, 'PHF1Tau']
ad_tau_high = excite_syn_ad.groupby(['Point']).mean().loc[:, 'PHF1Tau']
hc_tau_high = excite_syn_ad.groupby(['Point']).mean().loc[:, 'PHF1Tau']

ad_cd47_low = nonExcite_syn_ad.groupby(['Point']).mean().loc[:, 'CD47']
hc_cd47_low = nonExcite_syn_ad.groupby(['Point']).mean().loc[:, 'CD47']
ad_cd47_high = excite_syn_ad.groupby(['Point']).mean().loc[:, 'CD47']
hc_cd47_high = excite_syn_hc.groupby(['Point']).mean().loc[:, 'CD47']


syn = mibi.loc[((mibi.Synaptophysin>0)), :]




plt.figure()
mibi.loc[(mibi.Synaptophysin>0) & mibi.Point.isin([1, 2, 3, 4, 5, 6]) & (mibi.VGLUT1>0.0), 'VGLUT1'].hist(bins=100)
# plt.xlim(-0.001, 2)
plt.savefig('test.png')



plt.figure()
mibi.loc[(mibi.Synaptophysin>0) & mibi.Point.isin([7, 8, 9, 10, 11, 12]) & (mibi.VGLUT1>0.0), 'VGLUT1'].hist(bins=100)
# plt.xlim(-0.001, 2)
plt.savefig('test2.png')







