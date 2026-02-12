DRUG_DB = {
    "Control (No Treatment)": {
        "targets": [], "desc": "Baseline metabolism.", "mic_uM": 0, "typical_dose_uM": 0,
    },
    "Trimethoprim": {
        "targets": ["MG228"],
        "desc": "Inhibits DHFR (Folate synthesis). Competitive inhibitor with Ki ≈ 5 nM.",
        "mic_uM": 2.0, "typical_dose_uM": 10.0,
    },
    "Methotrexate": {
        "targets": ["MG228", "MG006"],
        "desc": "Blocks dTMP production. Tight-binding DHFR inhibitor (Ki ≈ 1 nM) + thymidylate synthase.",
        "mic_uM": 0.5, "typical_dose_uM": 5.0,
    },
    "Fosmidomycin": {
        "targets": ["MG066"],
        "desc": "Inhibits DXR in isoprenoid synthesis (MEP pathway). Slow tight-binding inhibitor.",
        "mic_uM": 10.0, "typical_dose_uM": 50.0,
    },
    "Cerulenin": {
        "targets": ["MG212", "MG114"],
        "desc": "Irreversible (covalent) inhibitor of fatty acid synthase. Time-dependent killing.",
        "mic_uM": 5.0, "typical_dose_uM": 20.0,
    },
    "Mupirocin": {
        "targets": ["MG345"],
        "desc": "Inhibits isoleucyl-tRNA synthetase. Competitive with isoleucine (Ki ≈ 20 nM).",
        "mic_uM": 0.1, "typical_dose_uM": 1.0,
    },
    "Generic Glycolysis Inhibitor": {
        "targets": ["MG041", "MG429"],
        "desc": "Blocks glucose import/phosphorylation. Cooperative inhibition.",
        "mic_uM": 20.0, "typical_dose_uM": 100.0,
    },
}