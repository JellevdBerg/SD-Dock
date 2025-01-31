import pandas as pd
import os
from qsprpred.models import QSPRModel
from qsprpred.utils.stopwatch import StopWatch
from tqdm import tqdm

output_folder = "./data/Predictions"
model = QSPRModel.fromFile("./data/Classification model/Q09472 (EP300) XGboost Model_meta.json")
df = pd.read_csv(
    "./data/2024.07_Enamine_REAL_DB_9.6M.cxsmiles.bz2",
    compression="bz2",
    chunksize=100_000,
    sep="\t"
)

for chunk_id, df_chunk in enumerate(tqdm(df, desc="Processing chunks")):
    print(f"Running chunk: {chunk_id}")
    watch = StopWatch()
    preds = model.predictMols(df_chunk["smiles"].to_list(), use_probas=True)
    df_preds = pd.DataFrame({
    "Prediction (active)": [p[1] for p in preds[0]],
    "ID": df_chunk["id"],
    "SMILES": df_chunk["smiles"]
})
    output_file = os.path.join(output_folder, "Preds_" + f"{chunk_id}".zfill(5) + ".csv")
    df_preds.to_csv(output_file, index=False, header=True)
    watch.stop()