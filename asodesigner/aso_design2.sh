#!/bin/bash

# Function to monitor system resources
monitor_resources() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] System Resource Usage:"
    echo "CPU Usage: $(top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1}')%"
    echo "Memory Usage: $(free -m | awk '/Mem/{printf "%.2f%%", $3/$2 * 100}')"
    echo "Disk Usage: $(df -h / | awk 'NR==2{print $5}')"
}

# Display current step and timestamp
print_step() {
    local STEP=$1
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Executing: $STEP"
}

# Error handling function
handle_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: Step $1 failed!"
    exit 1
}

# Main program
{
    # Receive parameters
    GENE_NAME="$1"
    ASO_LEN="$2"
    ASO_COUNT="$3"  # New ASO count parameter
    PRIORITY="$4"   # New priority parameter
    HOMOLOGOUS_SPECIES="$5"  # New homologous species parameter
    UUID="$6"       # New UUID parameter
    [ -z "$GENE_NAME" ] && { echo "Please provide gene name parameter"; exit 1; }
    [ -z "$ASO_LEN" ] && { echo "Please provide ASO length parameter"; exit 1; }
    [ -z "$ASO_COUNT" ] && { echo "Please provide ASO count parameter"; exit 1; }
    [ -z "$PRIORITY" ] && { echo "Please provide priority parameter"; exit 1; }
    [ -z "$HOMOLOGOUS_SPECIES" ] && { echo "Please provide homologous species parameter"; exit 1; }
    [ -z "$UUID" ] && { echo "Please provide UUID parameter"; exit 1; }  # Validate UUID parameter
    
    print_step "$UUID $ASO_LEN $ASO_COUNT $PRIORITY $HOMOLOGOUS_SPECIES $GENE_NAME"
    echo "[DEBUG] ASO_COUNT: $ASO_COUNT"
    echo "[DEBUG] PRIORITY: $PRIORITY"
    echo "[DEBUG] HOMOLOGOUS_SPECIES: $HOMOLOGOUS_SPECIES"
    
    # Create UUID-specific output directory
    UUID_DIR=${UUID}
    
    # Ensure UUID directory exists
    mkdir -p "$UUID_DIR"
    
    # Detect Python interpreter, prioritize miniconda3 python3, then system python3, then python
    PYTHON_CMD=$(which /root/miniconda3/bin/python3 2>/dev/null || which python3 2>/dev/null || which python 2>/dev/null)
    if [ -z "$PYTHON_CMD" ]; then
        echo "[ERROR] Python interpreter not found!"
        handle_error "Python detection"
    fi
    
    # Add miniconda3 bin directory to PATH to ensure commands are found
    export PATH="/root/miniconda3/bin:$PATH"
    
    echo "[DEBUG] Python path: $PYTHON_CMD"
    echo "[DEBUG] Python version: $($PYTHON_CMD --version)"
    echo "[DEBUG] Full PATH: $PATH"
    
    # Detect Python modules
    echo "[DEBUG] Checking Python modules..."
    
    # Check pandas
    $PYTHON_CMD -c "import pandas" 2>/dev/null
    if [ $? -ne 0 ]; then
        echo "[ERROR] Python module pandas is not installed!"
        handle_error "Python module detection"
    fi
    echo "[DEBUG] ✅ pandas module is installed"
    
    # Check other potentially required modules
    MODULES=("Bio" "numpy" "requests")
    for MODULE in "${MODULES[@]}"; do
        $PYTHON_CMD -c "import $MODULE" 2>/dev/null
        if [ $? -ne 0 ]; then
            echo "[WARNING] Python module $MODULE is not installed"
        else
            echo "[DEBUG] ✅ $MODULE module is installed"
        fi
    done
    
    # ---------------------------------------------------------------
    # Step 1: Extract gene GTF data
    print_step "Extracting gene GTF data (1/11)"
    GTF_FILE="$UUID_DIR/gene.gtf"  # Modified path
    cat /asodesigner/reference/human/genomic.gtf | grep "$GENE_NAME" > "$GTF_FILE" || handle_error "GTF extraction"
    
    # Find matching line
    MATCH_LINE=$(grep -E "gene_id[[:space:]]+\"$GENE_NAME\"" "$GTF_FILE" | awk -F'\t' '$3 == "gene"' | head -n 1)
    
    if [ -z "$MATCH_LINE" ]; then
        echo "Gene '$GENE_NAME' not found"
        exit 1
    fi
    
    # ---------------------------------------------------------------
    print_step "Parsing genome position (2/11)"
    # Extract key information
    CHR=$(echo "$MATCH_LINE" | awk -F'\t' '{print $1}')
    START=$(echo "$MATCH_LINE" | awk -F'\t' '{print $4}')
    END=$(echo "$MATCH_LINE" | awk -F'\t' '{print $5}')
    DIR=$(echo "$MATCH_LINE" | awk -F'\t' '{print $7}')
    GENE_ID=$(echo "$MATCH_LINE" | grep -o 'GeneID:[0-9]\+' | sed 's/GeneID://')
    echo $CHR $START $END $DIR $GENE_ID
    
    # Chromosome mapping
    declare -A CHR_MAP=(
        ["NC_000001"]="chr1"  ["NC_000002"]="chr2"  ["NC_000003"]="chr3"
        ["NC_000004"]="chr4"  ["NC_000005"]="chr5"  ["NC_000006"]="chr6"
        ["NC_000007"]="chr7"  ["NC_000008"]="chr8"  ["NC_000009"]="chr9"
        ["NC_000010"]="chr10" ["NC_000011"]="chr11" ["NC_000012"]="chr12"
        ["NC_000013"]="chr13" ["NC_000014"]="chr14" ["NC_000015"]="chr15"
        ["NC_000016"]="chr16" ["NC_000017"]="chr17" ["NC_000018"]="chr18"
        ["NC_000019"]="chr19" ["NC_000020"]="chr20" ["NC_000021"]="chr21"
        ["NC_000022"]="chr22" ["NC_000023"]="chrX"  ["NC_000024"]="chrY"
    )
    
    CHR_PREFIX=${CHR%%.*} # Extract NC number
    HUMAN_CHR=${CHR_MAP["$CHR_PREFIX"]}
    
    [ -z "$HUMAN_CHR" ] && {
        echo "Unknown chromosome: $CHR_PREFIX"
        exit 1
    }
    
    # Extract common prefix of start and end positions
    COMMON=""
    for ((i=0; i<${#START} && i<${#END}; i++)); do
        [ "${START:$i:1}" == "${END:$i:1}" ] && COMMON+="${START:$i:1}" || break
    done
    
    # ---------------------------------------------------------------
    print_step "Obtaining homologous genes (3/11)"
    echo $GENE_ID $UUID_DIR
    OUTPUT_FILE="$UUID_DIR/gene.gtf"  # Modified path
    echo -e "${CHR}\t${START}\t${END}" > "$OUTPUT_FILE"
    
    # Enhanced debugging information
    echo "[DEBUG] Homologous gene query details:"
    echo "[DEBUG] Current working directory: $(pwd)"
    echo "[DEBUG] GENE_ID: $GENE_ID"
    echo "[DEBUG] Does TaxId.csv file exist?: $(if [ -f /asodesigner/text/TaxId.csv ]; then echo "Yes"; else echo "No"; fi)"
    echo "[DEBUG] Does gene_orthologs file exist?: $(if [ -f /asodesigner/text/gene_orthologs ]; then echo "Yes"; else echo "No"; fi)"
    echo "[DEBUG] Does output directory exist?: $(if [ -d "$UUID_DIR/" ]; then echo "Yes"; else echo "No"; fi)"
    echo "[DEBUG] Does appdesignASO.py exist?: $(if [ -f /asodesigner/scr/appdesignASO.py ]; then echo "Yes"; else echo "No"; fi)"
    echo "[DEBUG] Is appdesignASO.py executable?: $(if [ -x /asodesigner/scr/appdesignASO.py ]; then echo "Yes"; else echo "No"; fi)"
    
    # Check file permissions
    echo "[DEBUG] File permissions:"
    ls -la /asodesigner/text/TaxId.csv 2>&1 || echo "[ERROR] Cannot view TaxId.csv permissions"
    ls -la /asodesigner/text/gene_orthologs 2>&1 || echo "[ERROR] Cannot view gene_orthologs permissions"
    ls -la /asodesigner/scr/appdesignASO.py 2>&1 || echo "[ERROR] Cannot view appdesignASO.py permissions"
    
    # Capture complete output
    echo "[DEBUG] Executing appdesignASO.py..."
    $PYTHON_CMD /asodesigner/scr/appdesignASO.py "$GENE_ID" /asodesigner/text/TaxId.csv /asodesigner/text/gene_orthologs "$UUID_DIR/" > "$UUID_DIR/appdesign_output.log" 2>&1
    
    # Check return code
    APPDESIGN_RESULT=$?
    echo "[DEBUG] appdesignASO.py return code: $APPDESIGN_RESULT"
    
    # If failed, output complete log
    if [ $APPDESIGN_RESULT -ne 0 ]; then
        echo "[ERROR] Homologous gene query failed, complete output:"
        cat "$UUID_DIR/appdesign_output.log"
        handle_error "Homologous gene query"
    else
        echo "[DEBUG] Homologous gene query successful, output:"
        cat "$UUID_DIR/appdesign_output.log"
    fi
    
    # ---------------------------------------------------------------
    print_step "Extracting SNP data (4/11)"
    
    # Enhanced debugging information
    echo "[DEBUG] SNP extraction details:"
    echo "[DEBUG] Does bigBedToBed command exist?: $(if [ -x $(which bigBedToBed 2>/dev/null) ]; then echo "Yes"; else echo "No"; fi)"
    echo "[DEBUG] bigBedToBed path: $(which bigBedToBed 2>/dev/null)"
    echo "[DEBUG] Does dbSnp155Common.bb file exist?: $(if [ -f /asodesigner/reference/human/dbSnp155Common.bb ]; then echo "Yes"; else echo "No"; fi)"
    echo "[DEBUG] dbSnp155Common.bb file size: $(if [ -f /asodesigner/reference/human/dbSnp155Common.bb ]; then ls -lh /asodesigner/reference/human/dbSnp155Common.bb | awk '{print $5}'; else echo "N/A"; fi)"
    echo "[DEBUG] HUMAN_CHR: $HUMAN_CHR"
    echo "[DEBUG] START: $START"
    echo "[DEBUG] END: $END"
    echo "[DEBUG] Output file: $UUID_DIR/snp.gtf"
    
    # Execute command and capture complete output
    /root/miniconda3/bin/bigBedToBed /asodesigner/reference/human/dbSnp155Common.bb \
        -chrom="$HUMAN_CHR" -start="$START" -end="$END" "$UUID_DIR/snp.gtf" > "$UUID_DIR/bigBedToBed_output.log" 2>&1
    
    # Check return code
    if [ $? -ne 0 ]; then
        echo "[ERROR] SNP extraction failed, detailed output:"
        cat "$UUID_DIR/bigBedToBed_output.log"
        handle_error "SNP extraction"
    else
        echo "[DEBUG] SNP extraction successful, output file size: $(ls -lh "$UUID_DIR/snp.gtf" | awk '{print $5}')"
    fi
    
    # ---------------------------------------------------------------
    print_step "Obtaining gene sequence (5/11)"
    /usr/bin/bedtools getfasta -fi /asodesigner/reference/human/genomic.fna \
        -bed "$UUID_DIR/gene.gtf" > "$UUID_DIR/gene.fa" || handle_error "Gene sequence extraction"  # Modified path
    
    # ---------------------------------------------------------------
    print_step "ASO design (6/11)"
    # All Python scripts get UUID directory parameter
    $PYTHON_CMD /asodesigner/scr/_1.py "$UUID_DIR" || handle_error "ASO step 1"
    $PYTHON_CMD /asodesigner/scr/_2.py "$ASO_LEN" "$DIR" "$CHR" "$UUID_DIR" || handle_error "ASO step 2"
    $PYTHON_CMD /asodesigner/scr/_4.py "$UUID_DIR"  || handle_error "ASO step 4"
    

    # ---------------------------------------------------------------
    print_step "RNA secondary structure prediction (8/11)"
    /root/miniconda3/bin/RNAfold \
        -d2 --noPS --noconv --noLP --infile "$UUID_DIR/forG.txt" > "$UUID_DIR/Gresult.txt" || handle_error "RNAfold"  # Modified path
    
    # ---------------------------------------------------------------
    print_step "Off-target analysis (9/11)"
    # All Python scripts get UUID directory parameter
    monitor_resources
    $PYTHON_CMD /asodesigner/scr/_5.py "$UUID_DIR" || handle_error "Off-target analysis"
    # Split query sequences
    split -l 50 "$UUID_DIR/foroff.fa" "$UUID_DIR/foroff_part_"
    
    # Run BLAST for each part
    for part in "$UUID_DIR/foroff_part_"*; do
        /root/miniconda3/bin/blastn -query "$part" \
            -db /asodesigner/reference/human/blast/humangene.blastdb \
            -task blastn-short -word_size 16 -evalue 1 \
            -out "$UUID_DIR/human_blstn_result_$(basename $part).txt" -outfmt 6
    done
    
    # Merge results
    cat "$UUID_DIR/human_blstn_result_"* > "$UUID_DIR/human_blstn_result.txt"
    
    # ---------------------------------------------------------------
    print_step "Exon analysis (10/11)"
    $PYTHON_CMD /asodesigner/scr/_6.py "$GENE_NAME" "$ASO_LEN" "$UUID_DIR" || handle_error "Off-target result analysis"
    $PYTHON_CMD /asodesigner/scr/_7.py "$UUID_DIR" || handle_error "Off-target result analysis 2"
    
    # ---------------------------------------------------------------
    print_step "Exon analysis (10/11)"
    grep -E "gene_id \"$GENE_NAME\".*\"GeneID:$GENE_ID\"|\"GeneID:$GENE_ID\".*gene_id \"$GENE_NAME\"" /asodesigner/reference/human/genomic.gtf > "$UUID_DIR/gene2.gtf"
    # All Python scripts get UUID directory parameter
    $PYTHON_CMD /asodesigner/scr/_71.py "$UUID_DIR" || handle_error "Exon analysis 1"
    /usr/bin/bedtools getfasta -fi /asodesigner/reference/human/genomic.fna \
        -bed "$UUID_DIR/exon.txt" > "$UUID_DIR/exonfasta.fa" || handle_error "Exon sequence extraction"
    $PYTHON_CMD /asodesigner/scr/_72.py "$DIR" "$UUID_DIR" || handle_error "Exon analysis 2"
    
    # Exon analysis 3 - detailed logging
    echo "[DEBUG] Exon analysis 3 parameters:"
    echo "[DEBUG] PYTHON_CMD: $PYTHON_CMD"
    echo "[DEBUG] ASO_LEN: $ASO_LEN"
    echo "[DEBUG] GENE_NAME: $GENE_NAME"
    echo "[DEBUG] UUID_DIR: $UUID_DIR"
    echo "[DEBUG] Does _73.py script exist?: $(if [ -f /asodesigner/scr/_73.py ]; then echo "Yes"; else echo "No"; fi)"
    if [ -f /asodesigner/scr/_73.py ]; then
        echo "[DEBUG] _73.py script permissions: $(ls -la /asodesigner/scr/_73.py)"
    fi
    echo "[DEBUG] Does UUID directory exist?: $(if [ -d "$UUID_DIR" ]; then echo "Yes"; else echo "No"; fi)"
    if [ -d "$UUID_DIR" ]; then
        echo "[DEBUG] UUID directory contents: $(ls -la "$UUID_DIR")"
    fi
    
    # Execute and capture complete output
    echo "[DEBUG] Starting exon analysis 3 script..."
    $PYTHON_CMD /asodesigner/scr/_73.py "$ASO_LEN" "$GENE_NAME" "$UUID_DIR" > "$UUID_DIR/_73_output.log" 2>&1
    EXIT_CODE=$?
    # ---------------------------------------------------------------
    print_step "ASO final selection (7/11)"
    $PYTHON_CMD /asodesigner/scr/_8.py "$UUID_DIR" "$ASO_COUNT" "$PRIORITY" "$HOMOLOGOUS_SPECIES" || handle_error "ASO final selection"
    
    echo "[DEBUG] Exon analysis 3 script completed, exit code: $EXIT_CODE"
    echo "[DEBUG] Script output:"
    cat "$UUID_DIR/_73_output.log"
    
    if [ $EXIT_CODE -ne 0 ]; then
        handle_error "Exon analysis 3"
    fi
    echo "[DEBUG] Exon analysis 3 executed successfully!"
    
    # ---------------------------------------------------------------
    print_step "ASO design completed! (12/12)"
    echo "==================== All steps completed ===================="
    date
    
    # Move final results to UUID directory
    mv /asodesigner/outfile/ASO_AllCandidates_*.xlsx "$UUID_DIR/" 2>/dev/null || true
    mv /asodesigner/outfile/*.png "$UUID_DIR/" 2>/dev/null || true
    mv /asodesigner/outfile/*.pdf "$UUID_DIR/" 2>/dev/null || true
} | tee "$UUID_DIR/aso_design_$(date +%Y%m%d_%H%M%S).log"  # Log file with timestamp to avoid overwriting

exit 0