# script should take an inputted flu subtype and download the clinical sequences for that subtype
# collection dates should match the time range of wastewater samples, be from a human host, in Texas, with complete assembly

from Bio import Entrez, SeqIO
import pandas as pd
from datetime import datetime

# ---- USER SETTINGS ----
Entrez.email = "u255582@bcm.edu"  # required by NCBI

# From wastewater script (YYYY-MM-DD format is best here)
wastewater_start = "2022-05-04"
wastewater_end = "2025-11-05"

virus = "Influenza A virus"
host = "Homo sapiens"
location = "USA: TX" 
subtype = input("Enter flu subtype (e.g. H1N1, H3N2, H5N1): ").strip() 

# Convert wastewater dates to the YYYY/MM/DD format used in the shell script
start_cd = wastewater_start.replace("-", "/")
end_cd = wastewater_end.replace("-", "/")

# ---- Build query for clinical sequences ----
# We use the simpler combination that we know returns hits in nuccore:
#   "Influenza A virus"[Organism] AND "Homo sapiens"[Host]
#   AND "USA: Texas"[All Fields] AND "<subtype>"[All Fields]
query = (
    '"Influenza A virus"[Organism] '
    'AND "Homo sapiens"[Host] '
    'AND "USA: Texas"[All Fields] '
    f'AND "{subtype}"[All Fields]'
)

print("NCBI query:")
print(query)

# ---- Search NCBI with the full query ----
print("\nSearching NCBI with full query...")
search_handle = Entrez.esearch(
    db="nuccore",
    term=query,
    retmax=10000,  # adjust if needed
)
search_result = Entrez.read(search_handle)
search_handle.close()

# Show total matches reported by NCBI
total_count = int(search_result.get("Count", 0))
print(f"NCBI reports {total_count} total matches for this query")

ids = search_result["IdList"]
print(f"Retrieving {len(ids)} sequences")

if not ids:
    exit(0)

# ---- Fetch sequences and metadata as GenBank ----
print("Fetching GenBank records for metadata and sequence filtering...")
fetch_handle = Entrez.efetch(
    db="nuccore",
    id=",".join(ids),
    rettype="gb",
    retmode="text",
)
records = list(SeqIO.parse(fetch_handle, "gb"))
fetch_handle.close()

# ---- Build metadata table ----
metadata_rows = []
for rec in records:
    accession = rec.id
    collection_date = None
    country = None
    host_qual = None
    subtype_qual = None

    for feat in rec.features:
        if feat.type == "source":
            q = feat.qualifiers
            collection_date = q.get("collection_date", [None])[0]
            country = q.get("country", [None])[0]
            host_qual = q.get("host", [None])[0]
            subtype_qual = q.get("serotype", [None])[0] if "serotype" in q else q.get("subtype", [None])[0]
            break

    metadata_rows.append(
        {
            "accession": accession,
            "collection_date": collection_date,
            "country": country,
            "host": host_qual,
            "subtype": subtype_qual,
        }
    )

md_df = pd.DataFrame(metadata_rows)

# ---- Filter metadata by collection date range ----
start_dt = datetime.fromisoformat(wastewater_start)
end_dt = datetime.fromisoformat(wastewater_end)

def parse_collection_date(cd: str):
    if not cd:
        return None
    for fmt in ("%Y-%m-%d", "%Y-%m", "%Y", "%d-%b-%Y"):
        try:
            return datetime.strptime(cd, fmt)
        except ValueError:
            continue
    return None

md_df["collection_dt"] = md_df["collection_date"].apply(parse_collection_date)

filtered_md = md_df[(md_df["collection_dt"].notna()) & (md_df["collection_dt"] >= start_dt) & (md_df["collection_dt"] <= end_dt)]

# ---- Write metadata CSV ----
md_csv_path = f"{subtype}_clinical_md.csv"
filtered_md.drop(columns=["collection_dt"], inplace=True)
filtered_md.to_csv(md_csv_path, index=False)
print(f"Filtered metadata written to {md_csv_path}")

# ---- Write FASTA for sequences within date range ----
keep_accessions = set(filtered_md["accession"].tolist())
filtered_records = [rec for rec in records if rec.id in keep_accessions]

output_fasta = f"{subtype}_clinical.fasta"
print(f"Writing filtered sequences to {output_fasta} ...")
with open(output_fasta, "w") as out_f:
    SeqIO.write(filtered_records, out_f, "fasta")

print("Done.")