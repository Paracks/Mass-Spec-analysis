import pandas as pd

# Load file
df = pd.read_csv("Final_with_TAIR_best_hits_localblast.csv")

# Identify rows where any field contains the word "mitochond"
mask = df.apply(lambda row: row.astype(str).str.contains("mitochond", case=False, na=False)).any(axis=1)

# Subset mitochondrial proteins
mito_df = df[mask]

# Save output
mito_df.to_csv("mitochondrial_interactors.csv", index=False)

print("Mitochondrial proteins saved to mitochondrial_interactors.csv")
print(mito_df.head())
