MEDIUM_PRESETS = {
    "Rich (Default)": {
        "desc": "SP4-like rich medium. All nutrients fully available.",
        "color": "#22c55e",
        "overrides": {},  # Use model defaults (all open at -10)
    },
    "Glucose-Rich": {
        "desc": "Excess glucose, normal amino acids. Simulates high-sugar environment.",
        "color": "#3b82f6",
        "overrides": {
            "EX_glc_D_b": -50.0,   # Flood glucose
            "EX_fru_b": -50.0,     # Fructose also available
            "EX_glyc_b": -50.0,    # Glycerol excess
        },
    },
    "Minimal (Starvation)": {
        "desc": "Only essential carbon/nitrogen. Amino acids severely limited.",
        "color": "#eab308",
        "overrides": {
            # Cut most amino acid dipeptides and amino acids to trace levels
            "EX_arg_L_b": -0.5,
            "EX_asn_L_b": -0.5,
            "EX_glu_L_b": -0.5,
            "EX_his_L_b": -0.5,
            "EX_ile_L_b": -0.5,
            "EX_leu_L_b": -0.5,
            "EX_lys_L_b": -0.5,
            "EX_met_L_b": -0.5,
            "EX_phe_L_b": -0.5,
            "EX_pro_L_b": -0.5,
            "EX_ser_L_b": -0.5,
            "EX_trp_L_b": -0.5,
            "EX_tyr_L_b": -0.5,
            "EX_val_L_b": -0.5,
            # Dipeptides severely limited
            "EX_ala_L_asp_L_b": -0.5,
            "EX_ala_L_gln_L_b": -0.5,
            "EX_ala_L_glu_L_b": -0.5,
            "EX_ala_L_his_L_b": -0.5,
            "EX_ala_L_leu_L_b": -0.5,
            "EX_ala_L_thr_L_b": -0.5,
            "EX_gly_asn_L_b": -0.5,
            "EX_gly_asp_L_b": -0.5,
            "EX_gly_gln_L_b": -0.5,
            "EX_gly_glu_L_b": -0.5,
            "EX_gly_met_L_b": -0.5,
            "EX_gly_pro_L_b": -0.5,
            "EX_met_L_ala_L_b": -0.5,
            # Nucleotide precursors limited
            "EX_datp_b": -1.0,
            "EX_dgtp_b": -1.0,
            "EX_gtp_b": -1.0,
            # Keep minimal carbon
            "EX_glc_D_b": -2.0,
            "EX_glyc_b": -2.0,
            "EX_glyc3p_b": -2.0,
        },
    },
    "No Nucleotides": {
        "desc": "Nucleotide precursors removed. Tests de novo synthesis capacity.",
        "color": "#a855f7",
        "overrides": {
            "EX_datp_b": 0.0,
            "EX_dgtp_b": 0.0,
            "EX_gtp_b": 0.0,
            "EX_thymd_b": 0.0,
            "EX_csn_b": 0.0,
            "EX_ura_b": 0.0,
            "EX_ade_b": 0.0,
        },
    },
    "No Lipid Precursors": {
        "desc": "Glycerol and glycerol-3-phosphate removed. Membrane synthesis stressed.",
        "color": "#f97316",
        "overrides": {
            "EX_glyc_b": 0.0,
            "EX_glyc3p_b": 0.0,
        },
    },
    "Iron/Metal Depleted": {
        "desc": "All metal ions removed. Metalloenzymes disabled.",
        "color": "#ef4444",
        "overrides": {
            "EX_fe2_b": 0.0,
            "EX_fe3_b": 0.0,
            "EX_mn2_b": 0.0,
            "EX_zn2_b": 0.0,
            "EX_cu2_b": 0.0,
            "EX_cobalt2_b": 0.0,
            "EX_ni2_b": 0.0,
            "EX_mg2_b": 0.0,
            "EX_mobd_b": 0.0,
            "EX_ca2_b": 0.0,
        },
    },
    "Vitamin Depleted": {
        "desc": "Cofactor/vitamin sources removed. Tests auxotrophy.",
        "color": "#ec4899",
        "overrides": {
            "EX_thm_b": 0.0,     # Thiamine
            "EX_ribflv_b": 0.0,  # Riboflavin
            "EX_nac_b": 0.0,     # Nicotinic acid
            "EX_pydx_b": 0.0,    # Pyridoxal
            "EX_fol_b": 0.0,     # Folate
            "EX_coa_b": 0.0,     # Coenzyme A
            "EX_mqn7_b": 0.0,    # Menaquinone
        },
    },
}