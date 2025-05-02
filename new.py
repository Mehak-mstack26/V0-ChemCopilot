import pandas as pd

# Read the Parquet file
df = pd.read_parquet('reaction_classification_checkpoint.parquet')

# Display the first 5 rows
pd.set_option('display.max_rows', None)

# Display all rows
# print(df)
print(df.columns)

