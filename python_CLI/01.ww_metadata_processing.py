# load modules
import pandas as pd
import glob
import os

# load in metadata
metadata_folder="/Users/camillemazurek2025/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/data2/metadata"
metadata_files=glob.glob(os.path.join(metadata_folder,"*.xlsx"))

md_list=[pd.read_excel(file) for file in metadata_files]

# merge metadata files into one large metadata file
### crm: not sure if I should ignore the index here
metadata=pd.concat(md_list, ignore_index=True)

# create a dictionary of expected cities and their public health regions
city_region = {
    "Houston, TX": "6_5S",
    "El Paso, TX": "9_10",
    "Lubbock, TX": "1",
    "Brownsville, TX": "11",
    "Wichita Falls, TX": "2_3",
    "Baytown, TX": "6_5S",
    "Humble, TX": "6_5S",
    "Missouri City, TX": "6_5S",
    "Austin, TX": "7",
    "Laredo, TX": "11",
    "Waco, TX": "7",
    "Fort Worth, TX": "2_3",
    "Palestine, TX": "4_5N",
    "Athens, TX": "4_5N",
    "Dallas, TX": "2_3",
    "Katy, TX": "6_5S"
}

# use dictionary to add region column to metadata
metadata["Region"] = metadata["City"].map(city_region)

# add warning for any unexpected cities, would need to update city_region dictionary
print("Unknown cities:", metadata.loc[metadata["Region"].isna(), "City"].unique())

# add a column for month_year to the metadata
metadata["Month_Year"] = metadata["Date"].dt.strftime("%m.%Y")

# organize the metadata in this specific way
metadata = metadata[
    [
        "City",
        "Sample_ID",
        "Site",
        "Date",
        "Flow",
        "PoolID",
        "SiteCode",
        "Region",
        "Month_Year",
    ]
]

# export metadata as tsv
metadata.to_csv("/Users/camillemazurek2025/Downloads/metadata_combined.csv", sep=",", index=False)

# print the time range of the wastewater samples
earliest_date = metadata["Date"].min()
latest_date = metadata["Date"].max()
print(f"The wastewater samples range from {earliest_date.strftime("%m/%d/%Y")} to {latest_date.strftime("%m/%d/%Y")}, you should use clinical data that matches this time range")

query = (
    '"Influenza A virus"[Organism] '
    'AND "Homo sapiens"[Host] '
    'AND Texas[Location] '
    f'AND {earliest_date.strftime("%m/%d/%Y")}[Collection Date] : {latest_date.strftime("%m/%d/%Y")}[Collection Date] '
    'AND "complete genome"'
)

print(f"Use this NCBI Virus query: {query}")