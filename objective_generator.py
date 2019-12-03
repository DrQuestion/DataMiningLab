import pandas as pd
import numpy as np

clinical = pd.read_csv ("clinical.tsv", header=0, sep="\t")
rows = []

for i, row in clinical.iterrows():
  #code for filtering Invidumed samples(with barcodes AA and AG)
  barcode = row[1].split("-")
  if barcode[1] != "AA" and barcode[1] != "AG":
  # Filter for dead and alive > 2years
    if row[8] == "Dead" and row[9] != "--":
      rows.append(i)
    else:
      if row[28] != "--" and int(row[28]) > 2*365: #calculated in days and not in years (days_to_last_followup)
        rows.append(i)

clinical = clinical.iloc[rows]
over2years = []
for i, row in clinical.iterrows():
  if row[8] == "Dead" and int(row[9]) <= 2*365:
    over2years.append(0)
  else:
    over2years.append(1)

objective = pd.DataFrame( list(zip(clinical.iloc[:,1], over2years)), columns=["barcode", "over2years"])
objective.to_csv("objective.tsv", sep="\t")
