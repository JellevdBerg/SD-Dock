import pandas as pd
import os

# set working directory
os.chdir("./Docking")

# subset the filtered data to the first 50 rows and only the columns we need
data = pd.read_csv("filtered_data.csv")
subset = data.iloc[:20000, [2]]
subset.to_csv("subset.csv", index=False)