import sys, os
import pandas as pd

# Get parameters
if len(sys.argv) < 5:
    print("Usage: python _8.py <UUID_DIR> <ASO_COUNT> <PRIORITY> <HOMOLOGOUS_SPECIES>")
    sys.exit(1)

uid = sys.argv[1]
aso_count = int(sys.argv[2])
priority = sys.argv[3]
homologous_species = sys.argv[4]

# Get gene name from the directory or filename
# Try to find the gene name from ASO_AllCandidates_*.txt files
gene_name = ""
for file in os.listdir(uid):
    if file.startswith("ASO_AllCandidates_") and file.endswith(".txt"):
        gene_name = file.split("ASO_AllCandidates_")[1].split(".txt")[0]
        break

if not gene_name:
    print("[ERROR] Could not find gene name from ASO_AllCandidates_*.txt files")
    sys.exit(1)

# Input and output files
input_file = f"{uid}/ASO_AllCandidates_{gene_name}.txt"
output_file = f"{uid}/ASO_FilteredCandidates_{gene_name}.xlsx"

print(f"[INFO] Running _8.py with parameters:")
print(f"[INFO] UUID directory: {uid}")
print(f"[INFO] ASO count: {aso_count}")
print(f"[INFO] Priority: {priority}")
print(f"[INFO] Homologous species: {homologous_species}")
print(f"[INFO] Gene name: {gene_name}")
print(f"[INFO] Input file: {input_file}")
print(f"[INFO] Output file: {output_file}")

# Read the input file
try:
    # Read the input file with proper column names
    df = pd.read_csv(input_file, sep="\t")
    
    print(f"[INFO] Successfully read input file with {len(df)} rows")
    print(f"[INFO] Columns: {list(df.columns)}")
    
    # Process based on priority
    if priority == "Homologous":
        print("[INFO] Filtering based on Homologous priority")
        
        # Map species name to column name
        species_map = {
            "crab-eating_macaque": "Macaca fascicularis (crab-eating macaque) homology",
            "mouse": "Mus musculus (mouse) homology",
            "rat": "Rattus norvegicus (rat) homology",
            "pig": "Sus scrofa (pig) homology",
            "rabbit": "Oryctolagus cuniculus (rabbit) homology",
            "guinea_pig": "Cavia porcellus (Guinea pig) homology"
        }
        
        if homologous_species == "none":
            print("[INFO] No homology filtering required")
            # Use all rows without filtering
            df_filtered = df.copy()
            print(f"[INFO] Using all {len(df_filtered)} rows")
        elif homologous_species.lower() == "all":
            print("[INFO] Filtering for rows with at least one species homology = 1.00")
            # Get all species columns
            species_columns = list(species_map.values())
            # Filter columns that exist in the DataFrame
            existing_species_columns = [col for col in species_columns if col in df.columns]
            
            if not existing_species_columns:
                print("[WARNING] No species homology columns found in input file, using all rows")
                df_filtered = df.copy()
                print(f"[INFO] Using all {len(df_filtered)} rows")
            else:
                # Create a filter that checks if any species has homology = 1.00
                filter_condition = False
                for col in existing_species_columns:
                    filter_condition |= (df[col] == 1.00)
                
                df_filtered = df[filter_condition]
                print(f"[INFO] Filtered to {len(df_filtered)} rows with at least one species homology = 1.00")
        elif homologous_species in species_map:
            print(f"[INFO] Filtering for rows with {homologous_species} homology = 1.00")
            species_column = species_map[homologous_species]
            
            if species_column in df.columns:
                # Filter rows where the specific species homology is 1.00
                df_filtered = df[df[species_column] == 1.00]
                print(f"[INFO] Filtered to {len(df_filtered)} rows with {homologous_species} homology = 1.00")
            else:
                print(f"[WARNING] Column for {homologous_species} not found in input file, using all rows")
                df_filtered = df.copy()
                print(f"[INFO] Using all {len(df_filtered)} rows")
        else:
            print(f"[WARNING] Invalid homologous species: {homologous_species}, using all rows")
            df_filtered = df.copy()
            print(f"[INFO] Using all {len(df_filtered)} rows")
    
    elif priority == "efficiency":
        print("[INFO] Filtering based on efficiency priority")
        
        # Define the columns to check
        mismatch_columns = [
            "0 mismatch genes",
            "1 mismatch genes",
            "2 mismatch genes",
            "3 mismatch genes",
            "4 mismatch genes"
        ]
        
        # Check if all mismatch columns exist
        for col in mismatch_columns:
            if col not in df.columns:
                print(f"[ERROR] Column '{col}' not found in input file")
                sys.exit(1)
        
        # Apply filters
        # 1. All mismatch counts <= 2
        mismatch_filter = (df[mismatch_columns] <= 2).all(axis=1)
        
        # 2. GC content between 0.25 and 0.6
        if "GC content" in df.columns:
            gc_filter = (df["GC content"] >= 0.25) & (df["GC content"] <= 0.6)
        else:
            print("[ERROR] Column 'GC content' not found in input file")
            sys.exit(1)
        
        # 3. ASO MFE >= -1
        if "ASO MFE" in df.columns:
            mfe_filter = df["ASO MFE"] >= -1
        else:
            print("[ERROR] Column 'ASO MFE' not found in input file")
            sys.exit(1)
        
        # Combine all filters
        df_filtered = df[mismatch_filter & gc_filter & mfe_filter]
        print(f"[INFO] Filtered to {len(df_filtered)} rows based on efficiency criteria")
        
        # Sort by RNase H score (ascending)
        if "RNase H score" in df.columns:
            df_filtered = df_filtered.sort_values(by="RNase H score", ascending=True)
            print("[INFO] Sorted by RNase H score (ascending)")
        else:
            print("[ERROR] Column 'RNase H score' not found in input file")
            sys.exit(1)
    
    else:
        print(f"[ERROR] Invalid priority: {priority}")
        sys.exit(1)
    
    # Limit to requested number of ASOs
    df_final = df_filtered.head(aso_count)
    print(f"[INFO] Final selection: {len(df_final)} rows")
    
    # Write to Excel file
    df_final.to_excel(output_file, index=False)
    print(f"[INFO] Final selection written to {output_file}")
    
    # Also create a text version for debugging
    text_output = f"{uid}/ASO_FilteredCandidates_{gene_name}.txt"
    df_final.to_csv(text_output, sep="\t", index=False)
    print(f"[INFO] Text version written to {text_output}")
    
except Exception as e:
    print(f"[ERROR] Error in _8.py: {str(e)}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
