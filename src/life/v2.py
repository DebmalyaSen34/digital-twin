import time
import sys
import random
import re
import cobra
from Bio import SeqIO
from Bio.Seq import Seq

# --- CONFIGURATION ---
GENOME_FASTA = "data/ncbi_dataset/data/GCA_000027325.1/cds_from_genomic.fna"
MODEL_FILE = "src/iPS189.xml"

def typing_effect(text, delay=0.005):
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)
    print()

def mutate_dna(dna_seq):
    """Introduces a mutation."""
    bases = ['A', 'T', 'C', 'G']
    seq_list = list(dna_seq)
    changes = []
    
    # Force a mutation
    pos = random.randint(0, len(seq_list)-1)
    original = seq_list[pos]
    new_base = random.choice([b for b in bases if b != original])
    seq_list[pos] = new_base
    changes.append(f"Position {pos}: {original} -> {new_base}")
    
    return Seq("".join(seq_list)), changes

def analyze_mutation_impact(original_prot, mutated_prot):
    if original_prot == mutated_prot:
        return "SILENT"
    elif "*" in mutated_prot[:-1]:
        return "NONSENSE" # Knockout
    else:
        return "MISSENSE" # Slow Enzyme

def run_strict_evolution():
    print("\n" + "="*60)
    print("   🧬 STRICT EVOLUTIONARY SIMULATOR: FORCED BOTTLENECK")
    print("="*60)
    
    print("[PHASE 1] LOADING METABOLIC MODEL...")
    model = cobra.io.read_sbml_model(MODEL_FILE)
    model.objective = "Biomass"
    # Ensure Omnivore diet to avoid instant death
    for rxn in model.exchanges: rxn.lower_bound = -10.0
    
    # 1. BASELINE CHECK (Find Active Genes)
    print("   > Calculating Healthy Flux...")
    solution = model.optimize()
    base_growth = solution.objective_value
    print(f"   ✅ Healthy Growth Rate: {base_growth:.5f}")
    
    # Get fluxes to see which genes are ACTUALLY working
    fluxes = solution.fluxes.abs()
    
    # Find candidates: Genes whose reactions have > 0 flux
    active_genes = []
    for gene in model.genes:
        # Check if ANY reaction controlled by this gene is active
        is_active = False
        for rxn in gene.reactions:
            if fluxes[rxn.id] > 0.0001:
                is_active = True
                break
        if is_active:
            cleaned_id = gene.id[:2] + "_" + gene.id[2:] if not gene.id.startswith("MG_") else gene.id
            active_genes.append(cleaned_id)
            
    print(f"   > Active Genes Detected: {len(active_genes)} (Skipping lazy ones)")
    
    # 2. PICK A TARGET (Only from active list)
    target_id = random.choice(active_genes)
    print(f"   > Target Selected: {target_id} (Carries significant flux)")
    
    # 3. LOAD GENOME DATA
    print(f"\n[PHASE 2] EXPOSING {target_id} TO MUTAGEN...")
    original_records = list(SeqIO.parse(GENOME_FASTA, "fasta"))
    candidate_record = None
    
    # Find the FASTA record matching our chosen ID
    for rec in original_records:
        if target_id in rec.description: # Simple substring match
            candidate_record = rec
            break
            
    if not candidate_record:
        print(f"   ❌ Error: {target_id} not found in FASTA. Using generic mutation.")
        # Proceed with simulation only
    else:
        # Apply Mutation
        original_dna = candidate_record.seq
        mutated_dna, logs = mutate_dna(original_dna)
        print(f"   ⚠️  MUTATION: {logs[0]}")
        
        # Check Impact
        orig_prot = original_dna.translate(table=4, to_stop=True)
        mut_prot = mutated_dna.translate(table=4, to_stop=True)
        impact = analyze_mutation_impact(orig_prot, mut_prot)
        
        print(f"   Mutation Type: {impact}")

        # 4. SIMULATE PENALTY
        print("\n[PHASE 3] APPLYING PHYSIOLOGICAL PENALTY...")
        target_gene = model.genes.get_by_id(target_id.replace("_", ""))
        
        if impact == "SILENT":
            print("   > Effect: None.")
        
        elif impact == "MISSENSE":
            # Here is the FIX: Instead of 50% capacity, limit to 10% of NEEDED flux.
            print("   > Effect: Enzyme damaged. Limiting to 10% capacity.")
            
            # Find the max flux this gene was carrying before
            max_flux_needed = 0
            affected_rxns = []
            for rxn in target_gene.reactions:
                current_flux = fluxes[rxn.id]
                if current_flux > max_flux_needed:
                    max_flux_needed = current_flux
                affected_rxns.append(rxn)
            
            # Set the limit TIGHTER than what is needed
            new_limit = max_flux_needed * 0.1
            
            for rxn in affected_rxns:
                # Constrain both directions
                if rxn.upper_bound > 0: rxn.upper_bound = new_limit
                if rxn.lower_bound < 0: rxn.lower_bound = -new_limit
                
        elif impact == "NONSENSE":
            print("   > Effect: Knockout (Zero capacity).")
            target_gene.knock_out()

        # 5. FINAL RESULT
        new_growth = model.optimize().objective_value
        
        print("\n" + "="*60)
        print("   EVOLUTIONARY OUTCOME")
        print("="*60)
        print(f"   Gene:       {target_id}")
        print(f"   Mutation:   {impact}")
        print(f"   Old Rate:   {base_growth:.5f}")
        print(f"   New Rate:   {new_growth:.5f}")
        
        loss = (base_growth - new_growth) / base_growth * 100
        if loss < 1:
            print("   Status:     ROBUST (Redundant pathways found)")
        elif loss > 99:
            print("   Status:     DEAD (Lethal Mutation)")
        else:
            print(f"   Status:     SICK ({loss:.1f}% Fitness Loss)")
        print("="*60 + "\n")

if __name__ == "__main__":
    run_strict_evolution()