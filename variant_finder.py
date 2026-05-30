from Bio.Seq import Seq

def variant_calling_pipeline(reference_seq, patient_seq):
    """
    Automated Variant Calling Engine (v2.0.0)
    Designed to align and detect Single Nucleotide Polymorphisms (SNPs)
    for clinical diagnostics assessment.
    """
    # 1. STANDARDIZATION & AUDIT
    ref = reference_seq.strip().upper()
    pat = patient_seq.strip().upper()
    
    if len(ref) != len(pat):
        raise ValueError("Sequence Alignment Error: Reference and Patient sequence lengths must match.")
    
    mutations_found = []
    total_bases = len(ref)
    
    # 2. MUTATION DETECTION ENGINE
    for position in range(total_bases):
        ref_base = ref[position]
        pat_base = pat[position]
        
        if ref_base != pat_base:
            mutation_entry = {
                "position": position + 1,  # Converting to biological 1-index
                "change": f"{ref_base} -> {pat_base}"
            }
            mutations_found.append(mutation_entry)
            
    # 3. CLINICAL REPORTING OUTPUT
    print(f"{'='*40}")
    print(f"🧬 CLINICAL VARIANT DETECTION REPORT 🧬")
    print(f"{'='*40}")
    print(f"Total Genomic Track Length: {total_bases} bases")
    print(f"Mutations Detected: {len(mutations_found)}")
    print(f"{'-'*40}")
    
    if mutations_found:
        for mut in mutations_found:
            print(f"🚨 ALERT: Variant at Position {mut['position']}: {mut['change']}")
    else:
        print("✅ ANALYSIS COMPLETE: 100% Match with Reference.")
        
    print(f"{'='*40}")
    
    return mutations_found

# --- SIMULATING CLINICAL RUN (TEST CASE) ---
healthy_ref = "ATGCGATCGATCGATC"
patient_dna = "ATGCGATCGCTCGATC"  # Mismatch at position 10 (A -> C)

detected_variants = variant_calling_pipeline(healthy_ref, patient_dna)
