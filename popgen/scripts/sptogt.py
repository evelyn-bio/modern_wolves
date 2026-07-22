import pandas as pd

df = pd.read_csv(
    "/projects/psg/people/pkb156/MW_sub2/globetrotter/test/testpaint_reftarget_chunklength.txt.gz",
    sep=r"\s+",
    index_col="indnames"
)

# Move row names into a column
df.reset_index(inplace=True)

# Remove the _0/_1 suffix
df["Recipient"] = df["indnames"].str.rsplit("_", n=1).str[0]

# Sum the two haplotypes for each individual
df = df.groupby("Recipient", as_index=False).sum(numeric_only=True)

# Rename the first column back to Recipient
df.rename(columns={"Recipient": "Recipient"}, inplace=True)

df.to_csv(
    "/projects/psg/people/pkb156/MW_sub2/globetrotter/test/testpaint_reftarget_chunklength_combined.txt",
    sep=" ",
    index=False
)