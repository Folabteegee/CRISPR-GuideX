# CRISPR-GuideX
An advanced, R-native bioinformatics platform for dynamic CRISPR gRNA design, featuring algorithmic support for next-gen nucleases (hfCas12Max), dual-guide deletion architectures, and bi-directional off-target prediction.

# 🧬 CRISPR GuideX (v3.1)

**CRISPR GuideX** is a production-level, local-first bioinformatics platform designed to accelerate pre-clinical gRNA design and therapeutic genome editing workflows. 

While legacy tools often restrict researchers to standard SpCas9 parameters and single-target workflows, GuideX is engineered for complex clinical architectures—including dual-guide hotspot deletions and reading-frame restoration surgeries.

### 🚀 Key Features

* **Multi-Nuclease Algorithmic Support:** Native design capabilities for SpCas9 (NGG), SaCas9 (NNGRRT), Cas12a (TTTV), hfCas12Max (TN/TNN), and fully custom regex PAM patterns.
* **Dual-Guide Deletion Architecture:** Batch-process multiple target sequences simultaneously to instantly generate the safest, highest-efficiency "Best Pair" upstream and downstream anchors.
* **Bi-Directional Off-Target Engines:** Weighted seed-region mismatch scoring that aggressively scans both the forward (+) and reverse (-) strands.
* **Dynamic Cut Geometry Simulation:** Mathematically predicts cut sites and outputs strand-aware visual maps of flanking DNA to assist with PCR primer design.
* **Absolute Patient Privacy:** Powered entirely by R and Shiny WebSockets with zero database retention. Patient clinical variant data processed via the NCBI fetch tool remains strictly in local browser memory.

### 🔬 Clinical Applications
GuideX v3.1 has been heavily optimized to support pre-clinical research for Duchenne Muscular Dystrophy (DMD), seamlessly handling deep intronic targeting and mutant allele sequence mapping.

### 💻 Built With
* R / Shiny
* `httr` & `jsonlite` (NCBI eUtils REST API Integration)
* `plotly` & `DT` (Interactive Data Visualization)

---
*Developed by Taiwo Afolabi for educational and pre-clinical research use.*
