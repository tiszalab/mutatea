import pandas as pd
df = pd.read_csv("/data/tisza/analyses/crm/mutatea/test_input_data/clinical_input_data_H1N1/sequences (2).csv", low_memory=False)
mask = pd.to_datetime(df["Collection_Date"], errors="coerce").isna()
print(df.loc[mask, "Collection_Date"].value_counts())