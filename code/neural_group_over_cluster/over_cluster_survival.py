import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style="ticks", font_scale=1)
plt.rcParams["font.family"] = "Inter"

from sklearn.model_selection import train_test_split
from lifelines.statistics import logrank_test
from lifelines.statistics import pairwise_logrank_test
from lifelines import KaplanMeierFitter

figdir = "figures"
if not os.path.exists(figdir):
	os.mkdir(figdir)
df_cl = pd.read_excel("clinical_cohort/metadata.xlsx", index_col=0, decimal=",")
df_cl = df_cl.set_index("IDAT-file")

df_cl["Neural Score"] = np.round(df_cl["Neural Score"], 2)
print(df_cl["Neural Score"].mean())

df_cl = df_cl.loc[~df_cl.index.duplicated(keep="first")]
print(df_cl.Source.value_counts())

ks = list(range(2,8))

df = pd.DataFrame(index=ks, columns=list(range(1)))
df1 = pd.DataFrame(index=ks, columns=list(range(1)))
plt.figure(figsize=(12,7))
cv=0
dict_pvals = {}
dict_results = {}
dict_breaks = {}
for i in ks:
    breaks = jenkspy.jenks_breaks(tmp['Neural Score'], n_classes=i)
    if i==2:
        breaks[1] = 0.41
    tmp_cl[f'groups_{i}'] = pd.cut(df_cl['Neural Score'], bins=breaks, labels=list(range(1, i+1)), 
                                  include_lowest=True, right=False, )
    groups = tmp_cl[f"groups_{i}"].copy()
    durations = tmp_cl["Overall survival [months]"]
    event_observed = tmp_cl["Dead at last follow-up [0=no, 1=yes, 2=n/a]"]
    tmp1 = pd.DataFrame({'durations': durations,
                           'events': event_observed,
                           'groups': groups})
    if i>2:
        results=pairwise_logrank_test(tmp1['durations'], tmp1['groups'], tmp1['events'])
    else:
        results=logrank_test(tmp1['durations'], tmp1['groups'], tmp1['events'])
        
    
    counts = {}
    labels = []
    ax = plt.subplot(2,3,i-1)
    for group in sorted(groups.unique()):
        kmf1 = KaplanMeierFitter()
        kmf1.fit(durations.loc[tmp_cl[f"groups_{i}"]==group], event_observed.loc[tmp_cl[f"groups_{i}"]==group], label=group)
        a1 = kmf1.plot(ci_show=True)
        if k==1:
            counts[1] = 215
            counts[2] = 150
        else:
            counts[group] = groups.value_counts()[group]
        labels.append(f"{group} (n={counts[group]})")
        
            
    if i-1 in [2,3,5,6]:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel("Probability of Survival [%]")
    if i-1 in [1,2,3]:
        ax.set_xticklabels([])
        ax.set_xlabel("")
    else:
        ax.set_xlabel("Overall Survivial [OS]")
    # kmf1.plot(ax=a1, ci_show=False)
    plt.title(f"Number of groups: {i}")
    handles, _ = ax.get_legend_handles_labels()
    plt.legend(handles, labels, title="Groups")
    ax.set_xlim(0,84)
    
    if i>2:
        sig_n = 0
        for j in results.p_value:
            if j < 0.05:
                sig_n+=1
        df1.loc[i, cv] = sig_n
        if sig_n==len(results.p_value):
            df.loc[i, cv] = "sig"
        else:
            df.loc[i, cv] = "NonSig"
    if i==2:
        if results.p_value < 0.05:
            df.loc[i, cv] = "sig"
            df1.loc[i, cv] = 2
        else:
            df.loc[i, cv] = "Nonsig"
            df1.loc[i, cv] = 0
    dict_pvals[i] = results.p_value
    dict_results[i] = results
    dict_breaks[i] = breaks
    
plt.suptitle(f"Clinical cohort (n={tmp_cl.shape[0]})")
plt.tight_layout()
plt.savefig("figures/jenks_breaks.pdf", bbox_inches="tight")
