import time
import sys
import random
import re
import cobra
from Bio import SeqIO
from Bio.Seq import Seq

GENOME_FASTA = "data/ncbi_dataset/data/GCA_000027325.1/cds_from_genomic.fna"
MODEL_FILE = "src/iPS189.xml"

def typing_effect(text, delay=0.005):
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)
    print()

def mutate_dna(dna_seq):
    """
    Introduces random point mutations into the DNA sequence to stimulate evolutionary processes.
    Returns: (Mutated Sequence, List of Changes)
    """
    bases = ['A', 'T', 'C', 'G']
    seq_list = list(dna_seq)
    changes = []

    stop_codons = ['TAA', 'TAG', 'TGA']

    if random.random() < 0.4 and len(seq_list) >= 6:
        num_codons = len(seq_list) // 3
        if num_codons > 2:
            target_codon_idx = random.randint(1, int(num_codons*0.8))
            codon_start = target_codon_idx * 3
            orignal_codon = seq_list[codon_start:codon_start+3]

            random.shuffle(stop_codons)
            mutated = False
            for stop in stop_codons:
                diffs = [(i, orignal_codon[i], stop[i]) for i in range(3) if orignal_codon[i] != stop[i]]
                if len(diffs) == 1:
                    i, orig_base, new_base = diffs[0]
                    pos = codon_start + i
                    seq_list[pos] = new_base
                    changes.append(f"Position {pos}: {orig_base} -> {new_base} (Introduced STOP codon)")
                    mutated = True
                    break
            
            if not mutated:
                stop = random.choice(stop_codons)
                i = random.randint(0, 2)
                pos = codon_start + i
                original_base = seq_list[pos]
                new_base = stop[i]
                if new_base == original_base:
                    # Change a different position in the codon
                    i = (i + 1) % 3
                    pos = codon_start + i
                    original_base = seq_list[pos]
                    new_base = stop[i]
                seq_list[pos] = new_base
                changes.append(f"Position {pos}: {original_base} -> {new_base}")
        else:
            pos = random.randint(0, len(seq_list) - 1)
            original_base = seq_list[pos]
            new_base = random.choice([b for b in bases if b != original_base])
            seq_list[pos] = new_base
            changes.append(f"Position {pos}: {original_base} -> {new_base}")
    else:
        pos = random.randint(0, len(seq_list) - 1)
        original_base = seq_list[pos]
        new_base = random.choice([b for b in bases if b != original_base])
        seq_list[pos] = new_base
        changes.append(f"Position {pos}: {original_base} -> {new_base}")

    # # Force at least one mutation for now
    # pos = random.randint(0, len(seq_list)-1)
    # original_base = seq_list[pos]
    # new_base = random.choice([b for b in bases if b != original_base])
    # seq_list[pos] = new_base
    # changes.append(f"Position {pos}: {original_base} -> {new_base}")

    return Seq("".join(seq_list)), changes

def analyze_mutation_impact(original_prot, mutated_prot):
    """
    Compares the healthy protein vs the mutated one.
    """

    if original_prot == mutated_prot:
        return "SILENT" # No change in amino acid sequence
    elif "*" in mutated_prot[:-1]:
        return "NONSENSE" # Premature stop codon introduced - Deadly
    else:
        return "MISSENSE" # Amino acid changed - maybe tolerable, maybe not
    
def run_evolution_simulator():
    print("\n" + "="*60)
    print("   🧬 GENOMIC MUTATION LAB: EVOLUTION SIMULATOR")
    print("="*60)
    
    # 1. LOAD BASELINE
    print("[PHASE 1] LOADING HEALTHY GENOME...")
    original_records = list(SeqIO.parse(GENOME_FASTA, "fasta"))
    print(f"   > Loaded {len(original_records)} genes.")
    
    # Pick a random gene to mutate (that actually exists in the metabolic model)
    # We load the model just to check which genes matter
    model = cobra.io.read_sbml_model(MODEL_FILE)
    model_gene_ids = [g.id for g in model.genes]
    
    # Find a candidate gene from your FASTA that maps to the model
    candidate = None
    candidate_record = None
    
    while not candidate:
        rec = random.choice(original_records)
        match = re.search(r'locus_tag=(MG_?\d+)', rec.description)
        print(match)
        if match:
            clean_id = match.group(1).replace("_", "")
            if clean_id in model_gene_ids:
                candidate = clean_id
                print(f"   > Found candidate gene: {candidate} for mutation.")
                candidate_record = rec


    print(f"   > Target Selected: {candidate} (Essential Enzyme)")
    
    # 2. APPLY RADIATION (MUTATION)
    print(f"\n[PHASE 2] EXPOSING {candidate} TO MUTAGENIC RADIATION...")
    time.sleep(1)
    
    original_dna = candidate_record.seq
    mutated_dna, logs = mutate_dna(original_dna)
    
    print(f"   ⚠️  MUTATION DETECTED: {logs[0]}")
    
    # 3. TRANSLATION ANALYSIS
    print("\n[PHASE 3] ANALYZING PROTEIN STRUCTURE...")
    # Translate both
    orig_protein = original_dna.translate(table=4, to_stop=True)
    mut_protein = mutated_dna.translate(table=4, to_stop=True)
    
    impact = analyze_mutation_impact(orig_protein, mut_protein)
    
    print(f"   Original Protein: {str(orig_protein)[:20]}...")
    print(f"   Mutated Protein:  {str(mut_protein)[:20]}...")
    print(f"   Mutation Type:    {impact}")

    # 4. METABOLIC SIMULATION
    print("\n[PHASE 4] SIMULATING PHYSIOLOGICAL EFFECT...")
    model.objective = "Biomass"
    # Open diet
    for rxn in model.exchanges: rxn.lower_bound = -10.0
    
    # Baseline Growth
    base_growth = model.optimize().objective_value
    print(f"   ✅ Healthy Growth Rate: {base_growth:.5f}")
    
    # Apply Damage
    if impact == "SILENT":
        print("   > Effect: None. The protein sequence is unchanged.")
        final_growth = base_growth
        
    elif impact == "MISSENSE":
        print("   > Effect: Enzyme Efficiency Reduced (50%).")
        # Find the reactions this gene controls and slow them down
        target_gene = model.genes.get_by_id(candidate)
        for rxn in target_gene.reactions:
            rxn.lower_bound = rxn.lower_bound * 0.5
            rxn.upper_bound = rxn.upper_bound * 0.5
        final_growth = model.optimize().objective_value
        
    elif impact == "NONSENSE":
        print("   > Effect: PROTEIN DESTROYED (Knockout).")
        target_gene = model.genes.get_by_id(candidate)
        # Use 'with model' context to temporarily kill it
        target_gene.knock_out()
        final_growth = model.optimize().objective_value
    
    # 5. FINAL REPORT
    print("\n" + "="*60)
    print("   EVOLUTIONARY OUTCOME REPORT")
    print("="*60)
    print(f"   Gene:           {candidate}")
    print(f"   Mutation:       {logs[0]}")
    print(f"   Classification: {impact}")
    print(f"   Survival:       {'YES' if final_growth > 0.000001 else 'NO'}")
    print(f"   New Growth:     {final_growth:.5f}")
    
    if final_growth < base_growth and final_growth > 0.000001:
        print("   Status:         SICK (Metabolic function impaired)")
    elif final_growth < 0.000001:
        print("   Status:         DEAD (Lethal Mutation)")
    else:
        print("   Status:         HEALTHY (Silent Mutation)")
    print("="*60 + "\n")

if __name__ == "__main__":
    run_evolution_simulator()