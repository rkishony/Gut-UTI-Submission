# Gut-pathobiome profiling predicts antibiotic resistance of urinary tract infections

Custom analysis code for the manuscript *"Gut-pathobiome profiling predicts antibiotic resistance of urinary tract infections"*.

The code is organised into three independent modules — **Phenotyping**, **Genotyping**, and **Clinic** — each described below after the shared setup.

## 1. System requirements

**Operating systems**
- Developed and tested on Windows 10/11 (MATLAB) and Linux/WSL2 (Ubuntu) for the Python pipeline.
- No non-standard hardware is required. A normal desktop/laptop is sufficient to run the Clinic demo and the Phenotyping statistics/figures. The full metagenomic pipeline (`Genotyping`) is compute-intensive and was run on a Linux server.

**MATLAB (`Clinic`, `Phenotyping`)**
- MATLAB **R2023a** (tested).
- **Statistics and Machine Learning Toolbox** (uses `fishertest`, `kstest2`, `prctile`, etc.).

**Python (`Genotyping`)**
- Python **3.11**.
- PyPI packages: see [`Genotyping/requirements.txt`](Genotyping/requirements.txt) (`numpy`, `pandas`, `matplotlib`, `scipy`, `biopython`).
- External command-line tools (invoked via `Genotyping/tools/`): cutadapt, bowtie2, seqtk, mash, MetaPhlAn, StrainPhlAn, ResFinder, PointFinder, MetaMLST, mlst, unicycler — plus their reference databases (MetaPhlAn/StrainPhlAn marker DB, ResFinder and PointFinder databases, MetaMLST bowtie2 index). Exact versions are listed in the manuscript Methods.

## 2. Installation guide

```bash
git clone https://github.com/rkishony/Gut-UTI-Submission.git
cd Gut-UTI-Submission
```

- **MATLAB modules** need no compilation. Each module adds its own paths via its `add_required_paths.m` (or, for the Clinic demo, `run_demo.m` does this automatically). Typical setup time: < 1 minute.
- **Python module**: create an environment and install dependencies:

  ```bash
  python -m venv venv && source venv/bin/activate   # or conda
  pip install -r Genotyping/requirements.txt
  ```

  Installing the external bioinformatics tools and their databases is the time-consuming step (tens of minutes to hours, mostly database downloads).

## 3. Phenotyping (MATLAB)

[`Phenotyping/`](Phenotyping/) starts from the analysed images of the individual micro-plating wells: the annotated PDFs and `results.csv` count tables, both provided in the manuscript **supplementary materials**. Layout: `loading_scripts/` (combine / MIC conversion), `analyses/` (isolate and community analyses), `metadata/`.

Download the community and isolate phenotyping packages and unpack them into:

```text
Phenotyping/source_data/Phenotyping/
  CommunityPhenotyping/
  CommunityPhenotypingB/
  IsolateResA/
  IsolateResB/
  IsolateResC/
```

Each experiment folder contains `results.csv` (colony counts / growth areas), annotated plate PDFs for verification, and layout metadata (`PlateLayout.xlsx`, `strips.xlsx`, `readme.xlsx`).

```matlab
cd Phenotyping
generate_phen_figures
```

This combines isolate experiments A/B/C, computes phenotypic-distance statistics, and writes figures under `Phenotyping/figures/`, `IsolateResData.csv`, and `phen_values.csv`.

## 4. Genotyping (Python + MATLAB)

[`Genotyping/`](Genotyping/) contains the metagenomic processing pipeline and figure scripts. Layout: `steps/` (pipeline stages: trim, downsample, assemble, mash, metaphlan, strainphlan, resfinder, metamlst, mlst, align, …), `tools/` (CLI wrappers), `utils/`, `scripts/` (Python figures), `matlab/` (MATLAB genotyping figures).

1. Download FASTQs from [`PRJNA1150262`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1150262) (faecal metagenomes) and [`PRJNA1470189`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1470189) (urine-isolate genomes). Use Supplementary Table 1 to match SRA records to study sample labels.
2. Install the external tools and databases listed under [System requirements](#1-system-requirements).
3. Create a local config:

   ```bash
   cd Genotyping
   cp config.example.py config.py
   # edit: FASTQ folders, TOOLS, DB paths, OUTPUT_DIR, sample sheets, ...
   ```

4. Install Python deps and run stages, then figures:

   ```bash
   pip install -r requirements.txt
   # each stage under steps/ via its run_cli, and/or:
   python scripts/create_all_figures.py
   ```

**Outputs** under `OUTPUT_DIR` include trimmed reads, taxonomic/strain profiles, ResFinder/PointFinder/MLST calls, and summaries. Two summaries are consumed by `Clinic`: `summary/trim_stats.csv` and `summary/alignment_stats_coalesced.csv`.


## 5. Clinic (MATLAB)

[`Clinic/`](Clinic/) contains the clinical and epidemiological analysis linking faecal (FIT) pathobiome resistance to urine-culture (UTI) resistance.
Clinical data are available from Maccabi Healthcare Services but restrictions apply to the availability of these data, which were used under license for the current study, and so are not publicly available. Access to the data is however available upon reasonable request and signing an MTA agreement with Maccabi Healthcare Services.
Entry point for the full analysis: `Clinic/generateFigures.m`.

Without access to the full clinical data, run the same analysis functions on synthetic data shipped in `Clinic/demo/` (structure matches the real inputs; numbers have no biological meaning):

```matlab
cd Clinic/demo
run_demo
```

`run_demo` creates `demo_data.mat` on first run via `make_synthetic_data.m` (`rng(1)`). To regenerate explicitly: `make_synthetic_data`.

- **Expected output**: `Clinic/demo/output/` with `values.xls` and 15 SVG figures (genotypic and phenotypic CDFs, risk bar charts, odds-ratio-vs-time, a resistance-gene coverage heatmap, and population statistics).
- **Expected run time**: under ~2 minutes on a normal desktop (demo uses 2,000 permutations/bootstraps for the odds-ratio figure; the full analysis uses 10,000).


## License

Released under the MIT License - see [`LICENSE`](LICENSE).
