# SV Caller Benchmark: Deletions & Insertions

This project benchmarks four structural variant (SV) callers—**Manta**, **DELLY**, **LUMPY**, and **Breakdancer**—on synthetic human genome data simulated with **SURVIVOR**. We evaluate how their algorithmic approaches affect performance and usability in detecting **deletions** and **insertions**.

## 🔍 Objective

To understand how different SV calling algorithms influence accuracy, efficiency, and usability in structural variant detection.

## 🧪 Methods

- **Data**: 10 simulated chr21 of GRCh38 with deletions and insertions (20–10,000 bp), generated using SURVIVOR and ART.
- **Tools**: Manta, DELLY, LUMPY, Breakdancer.
- **Metrics**: Precision, recall, F1-score, accuracy, and usability.

## 📁 Structure

├── VCF/ # Result vcf files
├── scripts/ # SV calling and evaluation scripts
├── fig/ # Plots and summaries
