import time
import sys
import cobra
from Bio import SeqIO
from Bio.Seq import Seq
import re

GENOME_FASTA = "data/ncbi_dataset/data/GCA_000027325.1/cds_from_genomic.fna"
MODEL_FILE = "src/iPS189.xml"

def typing_effect(text, delay=0.005):
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)
    print()

def visualize_central_dogma(record):
    """
    Real-time visualization of the central dogma processes: DNA -> RNA -> Protein
    """

    # Parse the header to find the locus_tag
    # Looking for something like "locus_tag=MG_001" in the header to extract a clean gene ID
    header = record.description
    match = re.search(r'locus_tag=(MG_?\d+)', header)

    if not match:
        return None
    
    raw_id = match.group(1) # e.g. MG_001
    clean_id = raw_id.replace("_", "") # e.g. MG001
    gene_name = 'Unknown'

    # Try to find the common name [gene=dnaN]
    name_match = re.search(r'gene=([^\]]+)', header)
    if name_match:
        gene_name = name_match.group(1)

    # Transcription
    dna_seq = record.seq
    mrna_seq = dna_seq.transcribe()

    # Translation
    # Use table=4 for Mycoplasma, but fallback to 11 if it fails (some genes might not follow the standard code)
    protein_seq = mrna_seq.translate(table=4, to_stop=True)

    # Visualization
    display_dna = str(dna_seq)[:20] + "..." + str(dna_seq)[-5:]
    display_prot = str(protein_seq)[:15] + "..."
    
    sys.stdout.write(f"\r🧬 READING: {clean_id:<6} ({gene_name:<5}) | DNA: {display_dna:<30}")
    time.sleep(0.01)
    sys.stdout.write(f"\r🏗️ BUILDING: {clean_id:<6} ({gene_name:<5}) | PEP: {display_prot:<30}\n")
    time.sleep(0.01)
    
    return clean_id

def run_genesis():
    print("\n")
    print("="*60)
    print("   GENESIS SIMULATOR: FROM RAW CODE TO LIFE")
    print("="*60)

    # --- PHASE 1: SEQUENCING & SYNTHESIS ---
    print(f"\n[PHASE 1] READING RAW GENOME ({GENOME_FASTA})...")
    time.sleep(1)
    
    synthesized_proteome = set()
    
    try:
        # Load the file
        fasta_sequences = list(SeqIO.parse(GENOME_FASTA, "fasta"))
        print(f"   > Loaded {len(fasta_sequences)} sequences.")
        print("   > Initiating Central Dogma (Transcription/Translation)...\n")
        time.sleep(1)
        
        # Process first 10 visually, then speed up
        for i, record in enumerate(fasta_sequences):
            gene_id = visualize_central_dogma(record)
            if gene_id:
                synthesized_proteome.add(gene_id)
            
            if i == 15:
                print(f"\n   ... [Accelerating Synthesis for {len(fasta_sequences)-15} remaining genes] ...")
        
        # Process the rest silently
        for record in fasta_sequences[16:]:
            match = re.search(r'locus_tag=(MG_?\d+)', record.description)
            if match:
                clean_id = match.group(1).replace("_", "")
                synthesized_proteome.add(clean_id)

    except FileNotFoundError:
        print("❌ Error: sequence.fasta not found.")
        return

    print(f"\n   ✅ PROTEOME ASSEMBLED: {len(synthesized_proteome)} proteins active.")
    time.sleep(1)

    # --- PHASE 2: METABOLIC DOCKING ---
    print("\n[PHASE 2] MAPPING PROTEINS TO METABOLIC NETWORK...")
    model = cobra.io.read_sbml_model(MODEL_FILE)
    
    # Check which model genes were actually found in your FASTA
    model_genes = {g.id for g in model.genes}
    found_genes = model_genes.intersection(synthesized_proteome)
    missing_genes = model_genes - synthesized_proteome
    
    print(f"   > Metabolic Enzymes Required: {len(model_genes)}")
    print(f"   > Enzymes Found in Genome:    {len(found_genes)}")
    print(f"   > Enzymes Missing:            {len(missing_genes)}")
    
    # CRITICAL STEP: Knockout the missing genes
    # If your FASTA is missing a gene, the simulated organism shouldn't have it!
    if missing_genes:
        print("   ⚠️  Removing missing enzymes from the simulation...")
        cobra.manipulation.delete_model_genes(model, list(missing_genes))
    
    # --- PHASE 3: IGNITION ---
    print("\n[PHASE 3] IGNITING METABOLISM...")
    model.objective = "Biomass"
    
    # Open the diet (Omnivore) to give it a chance
    for rxn in model.exchanges:
        rxn.lower_bound = -10.0
        
    try:
        solution = model.optimize()
        rate = solution.objective_value
        
        if rate > 0.000001:
            print("\n" + "="*60)
            print(f"   🚀 IT'S ALIVE! ")
            print("="*60)
            print(f"   Growth Rate: {rate:.5f} (Doublings per hour)")
            print(f"   Phenotype:   Viable")
            
            # Simple visual of the cell breathing
            print("\n   [Visualizing Metabolic Flux]")
            print("   Glucose  --> [ Glycolysis ] --> ATP")
            print("   No. Rxns --> [    " + str(len(model.reactions)) + "     ] --> Biomass")
        else:
            print("\n" + "="*60)
            print("   ☠️  ORGANISM FAILED TO GROW")
            print("="*60)
            print("   Reason: Critical metabolic genes were missing from the FASTA.")
            print("   (Check if your FASTA is the complete genome or just a fragment)")

    except Exception as e:
        print(f"   Error: {e}")

if __name__ == "__main__":
    run_genesis()