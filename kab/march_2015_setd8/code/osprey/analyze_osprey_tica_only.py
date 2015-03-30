import osprey.config, osprey.trials
import pandas as pd
import seaborn as sns

config = osprey.config.Config("./config.yaml")

session = config.trials()
items = [cursor.to_dict() for cursor in session.query(osprey.trials.Trial).all()]
df = pd.DataFrame(items).set_index('id')

for key in df.iloc[0].parameters.keys():
    df[key] = df.parameters.map(lambda x: x[key])


df[["tica__gamma", "mean_test_score"]]

plot(df.tica__gamma, df.mean_test_score, 'o')
xscale('log')
xlabel("gamma")
ylabel("GMRQ")
title("tICA / GMRQ cross validation (dihedrals)")
savefig("../figures/tica_cross_validation.png", bbox_inches="tight")
