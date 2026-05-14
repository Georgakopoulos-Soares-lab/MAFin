#!/usr/bin/env python3

import io
import json
import multiprocessing
import os
import shutil
import subprocess
import tempfile
import textwrap
import zipfile
import sys
import time
import psutil
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(page_title="MAFin", layout="wide")



CURRENT_DIR = Path(__file__).resolve().parent
TEST_DIR = CURRENT_DIR / "test"


def count_lines_fast(filepath):
    try:
        with open(filepath, 'rb') as f:
            return sum(1 for _ in f)
    except Exception:
        return None


def render_input_stats(label, content):
    if content is None:
        return

    if hasattr(content, "getvalue"):
        raw_bytes = content.getvalue()
        text_value = raw_bytes.decode("utf-8", errors="ignore")
    elif isinstance(content, bytes):
        raw_bytes = content
        text_value = raw_bytes.decode("utf-8", errors="ignore")
    else:
        text_value = str(content)
        raw_bytes = text_value.encode("utf-8")

    lines = text_value.splitlines()
    nonempty_lines = [line.strip() for line in lines if line.strip()]
    first_line = nonempty_lines[0] if nonempty_lines else "n/a"

    stat_cols = st.columns(4)

    if "mer" in label.lower():
        total_kmers = len(nonempty_lines)
        k_len = len(nonempty_lines[0]) if nonempty_lines else 0
        
        # Alphabet check (DNA: A, C, G, T, N)
        dna_chars = set("ACGTNacgtn")
        all_dna = all(set(line).issubset(dna_chars) for line in nonempty_lines) if nonempty_lines else True
        
        stat_cols[0].metric("Total K-mers", f"{total_kmers:,}", help="The total number of unique sequences to search.")
        stat_cols[1].metric("Length (K)", f"{k_len}", help="The number of nucleotides in each sequence.")
        stat_cols[2].metric("Alphabet", "DNA (ACGTN)" if all_dna else "Mixed/Invalid", help="Checks if the input contains only standard DNA bases.")
        stat_cols[3].metric("First Line", textwrap.shorten(first_line, width=28, placeholder="..."), help="A preview of the first valid line in the input.")

    elif "Regex" in label:
        import re
        total_patterns = len(nonempty_lines)
        
        # Validate regexes
        invalid_count = 0
        for p in nonempty_lines:
            try:
                re.compile(p)
            except:
                invalid_count += 1
        
        stat_cols[0].metric("Total Patterns", f"{total_patterns:,}", help="The total number of flexible patterns to evaluate.")
        stat_cols[1].metric("Input Size", f"{len(text_value):,} chars", help="Total number of characters in the regex input.")
        stat_cols[2].metric("Validation", "All Valid" if invalid_count == 0 else f"{invalid_count} Invalid", help="Checks if all patterns are mathematically valid regular expressions.")
        stat_cols[3].metric("First Line", textwrap.shorten(first_line, width=28, placeholder="..."), help="A preview of the first valid line in the input.")

    else:
        # Fallback for JASPAR or other files
        stat_cols[0].metric(f"{label} Size", f"{len(raw_bytes):,} bytes", help="The total file or pattern size in bytes.")
        stat_cols[1].metric("Lines", f"{len(lines):,}", help="The total number of lines in the input.")
        stat_cols[2].metric("Non-empty", f"{len(nonempty_lines):,}", help="The number of lines containing actual text/data.")
        stat_cols[3].metric("First Line", textwrap.shorten(first_line, width=28, placeholder="..."), help="A preview of the first valid line in the input.")


def execute_mafin_analysis(tmpdir, maf_file_obj, maf_name, search_type_param, search_data, search_in, rev_comp, procs, pval, bg_freqs, detailed, genome_ids_data=None):
    tmpdir_path = Path(tmpdir)
    maf_path = tmpdir_path / maf_name

    with open(maf_path, "wb") as file_handle:
        shutil.copyfileobj(maf_file_obj, file_handle)

    cmd = [sys.executable, str(CURRENT_DIR / "mafin.py"), str(maf_path)]

    if search_type_param == "kmers":
        kmers_path = os.path.join(tmpdir, "kmers.txt")
        with open(kmers_path, "wb") as file_handle:
            file_handle.write(search_data)
        cmd.extend(["--kmers", kmers_path])
    elif search_type_param == "regexes":
        patterns = [pattern.strip() for pattern in search_data.split("\n") if pattern.strip()]
        cmd.extend(["--regexes"] + patterns)
    elif search_type_param == "jaspar_file":
        jaspar_bytes, jaspar_name = search_data
        jaspar_path = os.path.join(tmpdir, jaspar_name)
        with open(jaspar_path, "wb") as file_handle:
            file_handle.write(jaspar_bytes)
        cmd.extend(["--jaspar_file", jaspar_path])

    cmd.extend([
        "--search_in",
        search_in,
        "--reverse_complement",
        rev_comp,
        "--processes",
        str(procs),
    ])

    if search_type_param == "jaspar_file":
        cmd.extend([
            "--pvalue_threshold",
            str(pval),
            "--background_frequencies",
        ] + bg_freqs.split())

    if detailed:
        cmd.append("--detailed_report")

    if genome_ids_data:
        genome_ids_path = os.path.join(tmpdir, "genome_ids.txt")
        with open(genome_ids_path, "wb") as f:
            f.write(genome_ids_data)
        cmd.extend(["--genome_ids", genome_ids_path])

    class AnalysisResult:
        def __init__(self, returncode, stdout, stderr, elapsed_time, max_memory_mb):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr
            self.elapsed_time = elapsed_time
            self.max_memory_mb = max_memory_mb

    start_time = time.time()
    process = subprocess.Popen(cmd, cwd=tmpdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    max_memory_mb = 0.0
    try:
        p = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                mem = p.memory_info().rss
                for child in p.children(recursive=True):
                    mem += child.memory_info().rss
                mem_mb = mem / (1024 * 1024)
                if mem_mb > max_memory_mb:
                    max_memory_mb = mem_mb
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            time.sleep(0.1)
    except psutil.NoSuchProcess:
        pass
        
    stdout, stderr = process.communicate()
    end_time = time.time()
    
    result = AnalysisResult(
        returncode=process.returncode,
        stdout=stdout,
        stderr=stderr,
        elapsed_time=end_time - start_time,
        max_memory_mb=max_memory_mb
    )
    return result, cmd


def render_analysis_results(workspace_dir, elapsed_time=None, max_memory=None):
    workspace = Path(workspace_dir)
    bed_files = list(workspace.glob("*_motif_hits.bed"))
    json_files = list(workspace.glob("*_motif_hits.json"))
    csv_files = list(workspace.glob("*_motif_hits.csv"))

    st.markdown("### Analysis Results & Execution Statistics")

    total_matches = sum((count_lines_fast(bed_file) or 0) for bed_file in bed_files)
    
    unique_chromosomes = set()
    for bed_file in bed_files:
        if bed_file.stat().st_size == 0:
            continue
        try:
            preview_df = pd.read_csv(bed_file, sep="\t", header=None, nrows=1000)
            if preview_df.shape[1] >= 1:
                unique_chromosomes.update(preview_df.iloc[:, 0].dropna().astype(str).tolist())
        except Exception:
            continue

    col_t1, col_t2, col_t3, col_t4 = st.columns(4)
    with col_t1:
        st.metric(label="Total Sequence Matches Found", value=f"{total_matches:,}", help="The total number of matches found across all searched genomes.")
    with col_t2:
        st.metric("Chromosomes Seen", f"{len(unique_chromosomes)}", help="The number of unique chromosomes or sequence identifiers that contained at least one motif match.")
    with col_t3:
        if elapsed_time is not None:
            st.metric("Time Elapsed", f"{elapsed_time:.2f} sec", help="The total time MAFin spent analyzing the data.")
    with col_t4:
        if max_memory is not None:
            st.metric("Peak Memory Usage", f"{max_memory:.2f} MB", help="The maximum amount of RAM consumed by the analysis engine.")

    st.markdown("<br>", unsafe_allow_html=True)
    summary_cols = st.columns(3)
    
    with summary_cols[0]:
        with st.container(border=True):
            st.metric("BED Files", f"{len(bed_files)}", help="BED (Browser Extensible Data) files contain the genomic coordinates (chromosome, start, end) of where the motif was found. They are standard in bioinformatics and can be directly loaded into genome browsers like UCSC or IGV.")
            if bed_files:
                if len(bed_files) == 1:
                    with open(bed_files[0], "rb") as f:
                        st.download_button("📥 Download", f, file_name=bed_files[0].name, mime="text/plain", key="dl_bed_top", use_container_width=True)
                else:
                    zip_buf = io.BytesIO()
                    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                        for bf in bed_files:
                            zf.write(bf, arcname=bf.name)
                    st.download_button("📦 Download All", zip_buf.getvalue(), file_name="bed_results.zip", mime="application/zip", key="dl_beds_top", use_container_width=True)

    with summary_cols[1]:
        with st.container(border=True):
            st.metric("CSV Reports", f"{len(csv_files)}", help="CSV (Comma-Separated Values) files provide a tabular, spreadsheet-friendly report of all hits. This includes detailed scoring, exact matches, and similarity vectors. Best for importing into Excel or statistical software for downstream analysis.")
            if csv_files:
                if len(csv_files) == 1:
                    with open(csv_files[0], "rb") as f:
                        st.download_button("📥 Download", f, file_name=csv_files[0].name, mime="text/csv", key="dl_csv_top", use_container_width=True)
                else:
                    zip_buf = io.BytesIO()
                    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                        for cf in csv_files:
                            zf.write(cf, arcname=cf.name)
                    st.download_button("📦 Download All", zip_buf.getvalue(), file_name="csv_results.zip", mime="application/zip", key="dl_csvs_top", use_container_width=True)

    with summary_cols[2]:
        with st.container(border=True):
            st.metric("JSON Output", f"{len(json_files)}", help="JSON (JavaScript Object Notation) files store the raw, hierarchical data from the run. This contains the most comprehensive structured output, perfect for programmatic access, parsing with scripts, or integrating into automated data pipelines.")
            if json_files:
                if len(json_files) == 1:
                    with open(json_files[0], "rb") as f:
                        st.download_button("📥 Download", f, file_name=json_files[0].name, mime="application/json", key="dl_json_top", use_container_width=True)
                else:
                    zip_buf = io.BytesIO()
                    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                        for jf in json_files:
                            zf.write(jf, arcname=jf.name)
                    st.download_button("📦 Download All", zip_buf.getvalue(), file_name="json_results.zip", mime="application/zip", key="dl_jsons_top", use_container_width=True)



    result_tabs = []
    if bed_files:
        result_tabs.append("BED Coordinates")
    if csv_files:
        result_tabs.append("CSV Report")
    if json_files:
        result_tabs.append("Raw JSON")

    if not result_tabs:
        st.info("Run completed, but no standard output files were detected.")
        return

    tabs = st.tabs(result_tabs)

    if bed_files:
        with tabs[result_tabs.index("BED Coordinates")]:
            for bed_file in bed_files:
                st.write(f"{bed_file.name}")
                if bed_file.stat().st_size == 0:
                    st.info(f"{bed_file.name} is empty (no matches found).")
                else:
                    bed_df = pd.read_csv(bed_file, sep="\t", header=None, nrows=1000)
                    if bed_df.shape[1] >= 6:
                        bed_df.columns = ["Chromosome", "Start", "End", "Motif", "Score", "Strand"]

                    bed_stats_cols = st.columns(4)
                    bed_stats_cols[0].metric("Preview Rows", f"{len(bed_df):,}", help="The number of data rows currently loaded in this visual preview.")
                    bed_stats_cols[1].metric("Columns", f"{bed_df.shape[1]}", help="The total number of columns in the generated BED file (usually 6: Chromosome, Start, End, Motif, Score, Strand).")
                    bed_stats_cols[2].metric("Unique Chromosomes", f"{bed_df['Chromosome'].nunique() if 'Chromosome' in bed_df.columns else 0}", help="The distinct number of chromosomes/sequence IDs found in this preview.")
                    bed_stats_cols[3].metric("Score Mean", f"{bed_df['Score'].astype(float).mean():.2f}" if 'Score' in bed_df.columns and not bed_df.empty else "n/a", help="The average conservation percentage across the motif matches shown here.")


                    bed_col_config = {
                        "Chromosome": st.column_config.TextColumn("Chromosome", help="The specific chromosome or sequence identifier where the motif was found."),
                        "Start": st.column_config.NumberColumn("Start", help="The exact 0-indexed starting position of the motif match on the chromosome.", format="%d"),
                        "End": st.column_config.NumberColumn("End", help="The exact ending position of the motif match (non-inclusive).", format="%d"),
                        "Motif": st.column_config.TextColumn("Motif", help="The exact sequence or motif name that was matched at this location."),
                        "Strand": st.column_config.TextColumn("Strand", help="Indicates if the match was on the forward (+) or reverse (-) DNA strand.")
                    }
                    if 'Score' in bed_df.columns:
                        bed_col_config["Score"] = st.column_config.ProgressColumn(
                            "Score",
                            help="Conservation percentage. 100% means the motif is perfectly conserved across all aligned genomes.",
                            format="%.2f%%",
                            min_value=0.0,
                            max_value=100.0,
                        )
                    st.dataframe(bed_df, use_container_width=True, column_config=bed_col_config)

                    if not bed_df.empty and "Chromosome" in bed_df.columns:
                        with st.expander("📊 View Match Distribution", expanded=True):
                            st.markdown("**Matches per Chromosome/Sequence**")
                            chrom_counts = bed_df["Chromosome"].value_counts()
                            st.bar_chart(chrom_counts)

                    total_rows = count_lines_fast(bed_file)
                    if total_rows is not None:
                        st.caption(f"Previewing first 1000 rows. Total rows: {total_rows}")
                    else:
                        st.caption("Previewing first 1000 rows.")


    if csv_files:
        with tabs[result_tabs.index("CSV Report")]:
            for csv_file in csv_files:
                st.write(f"{csv_file.name}")
                if csv_file.stat().st_size == 0:
                    st.info(f"{csv_file.name} is empty (no matches found).")
                else:
                    csv_preview = pd.read_csv(csv_file, nrows=1000)

                    csv_stats_cols = st.columns(3)
                    csv_stats_cols[0].metric("Preview Rows", f"{len(csv_preview):,}", help="The number of data rows currently loaded in this visual preview.")
                    csv_stats_cols[1].metric("Columns", f"{csv_preview.shape[1]}", help="The total number of columns present in the detailed CSV report.")
                    csv_stats_cols[2].metric("Preview Columns", ", ".join(csv_preview.columns[:3]) if len(csv_preview.columns) else "n/a", help="A quick reference to the first three column headers in the CSV file.")

                    csv_col_config = {
                        "Chromosome": st.column_config.TextColumn("Chromosome", help="The specific chromosome or sequence identifier where the motif was found."),
                        "Start": st.column_config.NumberColumn("Start", help="The exact 0-indexed starting position of the motif match on the chromosome.", format="%d"),
                        "End": st.column_config.NumberColumn("End", help="The exact ending position of the motif match (non-inclusive).", format="%d"),
                        "Motif": st.column_config.TextColumn("Motif", help="The exact sequence or motif name that was matched at this location."),
                        "Strand": st.column_config.TextColumn("Strand", help="Indicates if the match was on the forward (+) or reverse (-) DNA strand."),
                        "Similarity_Vector": st.column_config.TextColumn("Similarity_Vector", help="A visual string showing exactly which bases matched (1), mismatched (0), or were gaps (-) across the alignment block.")
                    }
                    if 'Score' in csv_preview.columns:
                        csv_col_config["Score"] = st.column_config.ProgressColumn(
                            "Score",
                            help="Conservation percentage. 100% means the motif is perfectly conserved across all aligned genomes.",
                            format="%.2f%%",
                            min_value=0.0,
                            max_value=100.0,
                        )
                    st.dataframe(csv_preview, use_container_width=True, column_config=csv_col_config)
                    total_rows = count_lines_fast(csv_file)
                    if total_rows is not None:
                        st.caption(f"Previewing first 1000 rows. Total rows: {total_rows}")
                    else:
                        st.caption("Previewing first 1000 rows.")


    if json_files:
        with tabs[result_tabs.index("Raw JSON")]:
            for json_file in json_files:
                st.write(f"{json_file.name}")
                size = json_file.stat().st_size
                if size == 0:
                    st.info(f"{json_file.name} is empty (no matches found).")
                elif size < 2_000_000:
                    with open(json_file, "r") as file_handle:
                        json_data = json.load(file_handle)
                        preview = json_data[:10] if isinstance(json_data, list) else json_data
                        json_stats_cols = st.columns(2)
                        if isinstance(json_data, list):
                            json_stats_cols[0].metric("Entries", f"{len(json_data):,}", help="The total number of data objects (entries) within the JSON array.")
                            json_stats_cols[1].metric("Preview Count", f"{len(preview):,}", help="The number of JSON objects currently displayed in this visual preview.")
                        else:
                            json_stats_cols[0].metric("Top-Level Type", type(json_data).__name__, help="The root structural type of the JSON data (e.g., dictionary or list).")
                            json_stats_cols[1].metric("Keys", f"{len(json_data.keys())}" if isinstance(json_data, dict) else "n/a", help="The total number of top-level keys in the JSON object.")
                        st.json(preview)
                else:
                    st.caption(f"JSON file too large to preview ({size} bytes). Download to inspect.")





with st.sidebar:
    st.image("https://camo.githubusercontent.com/58572d328397031d6eef8d5622024a6865c050331519a919a154a4fec7c33894/68747470733a2f2f692e706f7374696d672e63632f5635516e784d57702f6c6f676f2d6e6f2d6261636b67726f756e642e706e67", use_container_width=True)
    st.divider()
    nav_selection = st.radio(
        "Navigation",
        ["🚀 Run Motif Search", "⚡ Quick Start", "💡 Guides", "📖 Architecture", "📚 Glossary"],
        label_visibility="collapsed",
        help="Use this sidebar to navigate through the different sections of the MAFin application."
    )
    st.divider()
    st.markdown("**License**")
    st.caption("This project is licensed under the GNU GPL v3.")
    st.markdown("**Contact**")
    st.caption("For any questions or support, please contact:<br>izg5139@psu.edu<br>mpp5977@psu.edu<br>kap6605@psu.edu<br>ioannis.mouratidis@psu.edu", unsafe_allow_html=True)

if nav_selection == "⚡ Quick Start":
    st.title("⚡ Welcome to MAFin")
    
    st.markdown("### The Story Behind MAFin")
    st.markdown(
        """
        Imagine trying to find a specific phrase in a book, but the book is three billion letters long, 
        written in a four-letter alphabet (A, C, G, T), and you need to cross-reference it with the same 
        "book" from dozens of different species. 
        
        This is the daily reality of genomics. As scientists sequence more and more genomes, they align 
        them together to trace evolutionary history. These alignments are stored in massive files called 
        **Multiple Alignment Format (MAF)** files — chosen because MAF stores not just the DNA letters, but also the 
        precise coordinate mapping between each species' genome. Embedded within these files are the secrets of evolution—conserved 
        genetic sequences that might explain how diseases work, how species adapted, or how genes are regulated.
        
        However, exploring these massive files typically required complex command-line scripting and deep technical 
        expertise. The barrier to entry was high, locking away valuable biological insights from many researchers 
        and students.
        """
    )
    
    st.info("**MAFin (Multiple Alignment Format Motif Finder) was built to bridge this gap.**")
    
    st.markdown("### What MAFin Does")
    st.markdown(
        """
        MAFin is a high-performance search engine designed specifically for MAF files, wrapped in an easy-to-use 
        visual interface. It acts as your personal genomic detective.
        
        Instead of writing code, you simply tell MAFin what you are looking for:
        - **Exact Strings (K-mers):** A specific DNA sequence you know exists (the "K" means the sequence length).
        - **Flexible Patterns (Regex):** A sequence where some letters might vary.
        - **Biological Motifs (PWMs):** Complex binding profiles, such as those from the JASPAR database.
        
        MAFin instantly scans through the evolutionary alignments, pinpoints exactly where your sequence occurs 
        in the primary reference genome, and calculates how well that sequence has been conserved across all the 
        other species in the alignment block.
        """
    )
    
    st.markdown("### Our Purpose")
    st.markdown(
        """
        Our goal is to democratize genomic analysis. By combining a lightning-fast backend engine with an intuitive frontend, MAFin ensures that anyone—from seasoned bioinformaticians 
        to undergraduate students—can seamlessly explore the evolutionary history hidden within our DNA.
        """
    )
    
    st.divider()
    
    st.markdown("### Ready to explore?")
    col1, col2 = st.columns(2)
    with col1:
        with st.container(border=True):
            st.markdown("#### 1️⃣ Learn the Terms")
            st.caption("Understanding the vocabulary helps you interpret results correctly. Check the **📚 Glossary** for quick definitions.")
    with col2:
        with st.container(border=True):
            st.markdown("#### 2️⃣ Start Searching")
            st.caption("Go to **🚀 Run Motif Search** to upload your own `.maf` files.")
            
    st.markdown("<br>", unsafe_allow_html=True)
    st.success("🎉 **System Ready:** MAFin is successfully installed and fully operational on your system.")

elif nav_selection == "💡 Guides":
    st.title("🧬 Understanding MAFin: Complete Workflow")
    
    # Create tabs for detailed scenarios
    tab_intro, tab_maf, tab1, tab2, tab3, tab_faq = st.tabs(["📖 Introduction & Setup", "📂 Reading MAF Files", "🔍 Scenario 1: Exact K-mers", "🧬 Scenario 2: Flexible Patterns", "📊 Scenario 3: Motif Profiles", "❓ Common Questions"])
    
    with tab_intro:
        st.subheader("📖 Introduction & Setup Guide")
        st.markdown("""
        **MAFin** (Motif Detection in Multiple Alignment Files) is a powerful bioinformatics tool that searches for biological motifs in **multiple genome alignments**. 
        Here's what it does:
        
        ### The Big Picture
        1. **Input**: You provide a MAF file (Multiple Alignment Format) containing DNA sequences aligned across multiple species/strains
        2. **Process**: MAFin searches for specific motifs (DNA patterns) you're interested in. Motifs mark functional regions — gene switches, binding sites, mutation signatures — so finding them reveals which parts of DNA actually *do* something.
        3. **Output**: It tells you exactly WHERE each motif was found, including the genomic coordinates, strand, and conservation scores
        
        Think of it like a "Find" function (Ctrl+F) in a text editor, but for DNA sequences across multiple aligned genomes!
        """)
    
        st.markdown("""
        ### 🚀 Step-by-Step Guide: How to Use MAFin
        
        If you are new to MAFin, follow these simple steps to run your first analysis:
        """)

        with st.expander("Step 1: Upload your Data", expanded=True):
            st.markdown("""
            - Go to the **🚀 Run Motif Search** tab.
            - Click **"Browse files"** to upload your `.maf` file. This file contains the aligned genomes you want to search. If you don't have one, you can download the provided *Sample Data*.
            """)

        with st.expander("Step 2: Choose your Search Method", expanded=False):
            st.markdown("""
            - **Exact K-mers:** Type exact DNA sequences (e.g., `ATGC`). Best for finding specific, known strings.
            - **Regex:** Type flexible patterns (e.g., `A[GT]C`). Best for motifs with known variations across species.
            - **PWM (JASPAR):** Upload a `.jaspar` matrix file. JASPAR is a public database of experimentally-determined binding profiles. A matrix is used because real biology isn't black-and-white — at each position, some letters are preferred but others are tolerated.
            """)

        with st.expander("Step 3: Adjust Settings (Optional)", expanded=False):
            st.markdown("""
            - **Search Target:** "Reference Sequence" (faster — sufficient if you expect the motif in the primary genome) or "All Aligned Sequences" (use when the motif might only exist in another species).
            - **Strands to Search:** DNA is double-stranded. Choose "Both Strands (+/-)" if you're unsure which strand your motif sits on — this ensures you don't miss matches on the opposite strand.
            """)

        with st.expander("Step 4: Run Analysis & Download Results", expanded=False):
            st.markdown("""
            - Click the **"🚀 Run MAFin Analysis"** button.
            - Behind the scenes, MAFin scans the massive genome blocks, matches your patterns, and calculates how conserved they are across different species.
            - Once finished, you can view interactive tables of your results and download them as **BED** (coordinates), **CSV** (detailed report), or **JSON** files.
            """)
    
        st.markdown("""
        ### 📊 Comparing the Three Methods
    
        | Feature | K-mers | Regex | PWM |
        |---------|--------|-------|-----|
        | **Flexibility** | None (exact only) | Medium (pattern rules) | High (probability-based) |
        | **Speed** | ⚡ Fastest | ⚡⚡ Medium | ⚡⚡⚡ Slower |
        | **Complexity** | Simple | Medium | Complex |
        | **Best For** | Known exact sequences | Motifs with variants | Biological motif discovery |
        | **Example** | Primer sequences | Degenerate binding sites | TFBS from JASPAR |
        """)
    
    def highlight_maf(snippet, highlights):
        html = snippet
        for h in highlights:
            html = html.replace(h, f"<span style='background-color: #fde047; color: #0f172a; font-weight: 800; padding: 0 2px; border-radius: 2px;'>{h}</span>")
        html = html.replace('\n', '<br>')
        return f'<div style="background-color: #0e1117; color: #c9d1d9; padding: 1rem; border-radius: 0.5rem; font-family: monospace; font-size: 14px; white-space: nowrap; overflow-x: auto; border: 1px solid #30363d; margin-bottom: 1rem;">{html}</div>'

    maf_snippet = textwrap.dedent("""\
    ##maf version=1 scoring=blastz
    a score=19253.6
    s ref_genome.chr29        25671 1950 + 1000000   AAAATTTTGGGGCCCCACCTGGCAGGGCAGTCCGAATGGGCCAGCAAGTGGAGATGACT-AGCCCGGCT-GCTCGACCCGACCGCATCAAATAATCGGGGCGAAG-TTCGATCAGTCAATCGG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATCCTATATGACGCGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_3.chr2     25710 1953 - 27402155  AAAATTTTGGGGCCCCACCTGGCAGGGCAGTCCGAATGGGCCAGCAAGTGGAGATGACT-AGCCCGGCT-GCTCGACCCGACCGCATCAAATAATCGGGGCGAAG-TTCGATCAGTCAATCGG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATCCTATATGACGCGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_5.chr2     25685 1951 + 41169554  AAAATTTTGGGGCCCCACCTGGCAGGGCAGTCCGAATGGGCCAGCAAGTGGAGATGACT-AGCCCGGCT-GCTCGACCCGACCGCATCAAATAATCGGGGCGAAG-TTCGATCAGTCAATCGG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATCCTATATGACGCGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_7.chr2     25715 1950 + 64160004  AAAATTTTGGGGCCCCACCTGGCAGGGCAGTCCGAATGGGCCAGCAAGTGGAGATGACT-AGCCCGGCT-GCTCGACCCGACCGCATCAAATAATCGGGGCGAAG-TTCGATCAGTCAATCGG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATCCTATATGACGCGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT

    a score=18945.3
    s ref_genome.chr29        24940 1800 + 1000000   TCGATAAAATTTTGGGGCCCCACCTGGCAGGGCATTCCGTAATGGGCCAGCAAGTGGAGATGACT-AGCCCGGCT-GCTCGACCCGACCGCATCAAATAATCGGGGCGAAG-TTTTTTTTTTTTTTTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATCCCATATATAGCGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_3.chr2     24980 1803 - 27402155  TCGATAAAATTTTGGGGCCCCACCTGGCAGGGCATTCCGTAATGGGCCAGCAAGTGGAGATGACT-AGCCCGGCT-GCTCGACCCGACCGCATCAAATAATCGGGGCGAAG-TTTTTTTTTTTTTTTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATCCCATATATAGCGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_5.chr2     24955 1801 + 41169554  TCGATAAAATTTTGGGGCCCCACCTGGCAGGGCATTCCGTAATGGGCCAGCAAGTGGAGATGACT-AGCCCGGCT-GCTCGACCCGACCGCATCAAATAATCGGGGCGAAG-TTTTTTTTTTTTTTTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATCCCATATATAGCGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT

    a score=15022.1
    s ref_genome.chr24        19300 1500 + 1000000   TCGACCCGACCGCATCAATTTTTTTTTTTTTTTTTGGGGGCGAAG-AGTCCCCCGTTAATTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGCATATA-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_3.chr5     19350 1502 - 12402155  TCGACCCGACCGCATCAATTTTTTTTTTTTTTTTTGGGGGCGAAG-AGTCCCCCGTTAATTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGCATATA-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_5.chrX     10520 1500 + 31169554  TCGACCCGACCGCATCAATTTTTTTTTTTTTTTTTGGGGGCGAAG-AGTCCCCCGTTAATTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGCATATA-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT

    a score=13000.5
    s ref_genome.chr29        26400 1500 + 1000000   TCGACCCGACCGCATCAAAAATTTTTTTTTTTTTTTTGGGGGCGAAG-AGTCCCCCGTTAATTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_3.chr2     26440 1502 - 27402155  TCGACCCGACCGCATCAAAAATTTTTTTTTTTTTTTTGGGGGCGAAG-AGTCCCCCGTTAATTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_5.chr2     26415 1500 + 41169554  TCGACCCGACCGCATCAAAAATTTTTTTTTTTTTTTTGGGGGCGAAG-AGTCCCCCGTTAATTTG--ACCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT


    a score=9876.5
    s ref_genome.chr18        10200 1500 + 1000000   GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_3.chr4     10240 1502 - 12402155  GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT
    s fake_species_5.chr7     10215 1500 + 31169554  GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCCATCTAGAAACTCTCCGATGTTTCCGCCGTTCACATAACCCTACACTACATTAATTTACGACGAGGTTAAACAGATGCCCTGAACGAACTCAGTTACGACTGTACACGCAA-GTGTTACCTGTAGGCATCCCCCAACGGCCTC-ACAAGTTCCC-CGGGCACTACGGGGGTTCTAGAA-GGTATATG-ATAAGTTCTTTGATCTGTTCGTGTACTCCAGACTCAACATTCGCGAACTCGCCCATGAGTATGCATGTT""")
    
    with tab_maf:
        st.subheader("📂 How to Read a MAF File")
        st.markdown("""
        Before running any search, it helps to understand the file you're searching through.
        A **MAF file** is simply a text file that stores the results of aligning DNA from multiple species.
        
        💡 **Where can I find MAF files?**
        The most common source is the **[UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html)**. Navigate to your species (e.g., Human), then the "Downloads" page, and look for the **"Multiple Alignments"** section (often labeled as `multiz` or `n-way` alignments).
        """)

        with st.expander("1️⃣ The Header Line: `##maf`", expanded=False):
            st.markdown("""
            Every valid MAF file begins with a header line that looks like this:
            """)
            st.code("##maf version=1 scoring=blastz", language="text")
            st.markdown("""
            **What it is:** The first line of the file, always starting with `##maf`. It specifies the format version and optionally the scoring scheme used to generate the alignments.
            
            **Why it exists:** Bioinformatics tools need to verify they are reading the correct file type. Without this header, a tool might try to read a completely different file as a MAF file, causing crashes. It also provides metadata (like `scoring=blastz`) to tell researchers *how* the alignment was created.
            
            **How it is used:** MAFin checks this exact line to confirm the file is valid before it begins searching. If this line is missing, MAFin will reject the file to prevent processing invalid data.
            """)

        with st.expander("2️⃣ What is an Alignment Block?", expanded=True):
            st.markdown("""
            Imagine you have the DNA of a human, a mouse, and a dog. Scientists have compared their genomes
            and found regions where the DNA is similar — these similar regions are **aligned** together.

            Each aligned region is called an **alignment block**. A MAF file is just a long list of these blocks,
            one after another. Each block starts with a line beginning with `a` (for "alignment") and contains
            one `s` line (for "sequence") per species:
            """)
            st.code("""a score=19253.6
s ref_genome.chr29   25671 1950 + 1000000   AAAATTTTGGGGCCCC...
s fake_species_3.chr2  25710 1953 - 27402155  AAAATTTTGGGGCCCC...
s fake_species_5.chr2  25685 1951 + 41169554  AAAATTTTGGGGCCCC...""", language="text")
            st.markdown("""
            - The `a` line marks the start of a block and carries a quality score (higher = better alignment quality, meaning more matching letters and fewer gaps).
            - Each `s` line is **one species' DNA** in that aligned region.
            - The first `s` line is always the **reference genome** — the primary species everything else is compared to.
            """)

        with st.expander("3️⃣ Reading an `s` Line — Field by Field", expanded=False):
            st.markdown("""
            Every `s` line has 7 fields. Instead of just listing them, let's understand **why each one exists**
            — what question does it answer, and what would be missing without it?
            """)
            st.code("s ref_genome.chr24   17744 2457 + 1000000   GGCTACGACGA...", language="text")

            st.markdown("---")
            st.markdown("##### Field 1: `s` — Line Type")
            st.markdown("""
            **What it is:** The letter `s` at the start of the line.

            **Why it exists:** MAF files can contain other line types (`a` for alignment headers, `#` for comments,
            `q` for quality scores, etc.). The `s` tells any program reading the file: *"this line contains
            actual DNA sequence data."*
            """)

            st.markdown("---")
            st.markdown("##### Field 2: `ref_genome.chr24` — Source")
            st.markdown("""
            **What it is:** The identity of this sequence — which species and which chromosome.

            **Why it exists:** A MAF file can contain dozens of species, each with many chromosomes.
            Without this field, you'd have no idea *whose DNA* you're looking at. The format is always
            `species.chromosome`, separated by a dot:
            - `ref_genome` → the species/genome name
            - `chr24` → the specific chromosome within that genome

            **Why the dot matters:** It lets MAFin (and other tools) split the field programmatically.
            `ref_genome.chr24` and `fake_species_3.chr24` both say "chr24" but belong to completely
            different organisms. The species prefix before the dot is what keeps them apart.
            """)

            st.markdown("---")
            st.markdown("##### Field 3: `17744` — Start Position")
            st.markdown("""
            **What it is:** The position on the chromosome where this aligned region begins.

            **Why it exists:** Chromosomes are millions of letters long, but each alignment block only covers
            a tiny slice. The start position tells you *where* on the chromosome this slice sits.
            Without it, you'd know the DNA letters but have no idea where they came from on the chromosome.

            **Important detail:** This number is **0-indexed** (counting starts at 0, not 1). This is a standard
            convention in computer science and genomics. The very first letter on a chromosome is at position 0.
            """)

            st.markdown("---")
            st.markdown("##### Field 4: `2457` — Size (Number of Actual Letters)")
            st.markdown("""
            **What it is:** How many real DNA letters (A, C, G, T) are in this sequence — **not counting gaps**.

            **Why it exists:** The sequence field (field 7) contains alignment gaps (`-`) that were inserted
            to keep species lined up. So if the sequence looks like `ACGT--ACG`, the *size* is 7 (only counting
            actual letters), even though the displayed sequence is 9 characters long. You need this field to know
            how many real nucleotides this block covers.

            **How it connects to Start:** Together, Start + Size tells you the exact span on the chromosome.
            Start = 17744, Size = 2457 means this block covers positions **17744 through 20200** on the chromosome.
            """)

            st.markdown("---")
            st.markdown("##### Field 5: `+` — Strand")
            st.markdown("""
            **What it is:** Which direction the sequence is read — forward (`+`) or reverse complement (`-`).

            **Why it exists:** DNA is double-stranded. A gene can be encoded on either strand. When a sequence
            is on the `-` strand, the letters in the MAF file are the **reverse complement** of what sits on
            the chromosome. This field tells tools how to correctly map positions back to the genome.

            **What happens with `-` strand:** The coordinates are counted from the *end* of the chromosome
            (relative to the Source Size). This is why the Source Size field is essential — you need it to
            convert `-` strand coordinates back to forward-strand positions.
            """)

            st.markdown("---")
            st.markdown("##### Field 6: `1000000` — Source Size")
            st.markdown("""
            **What it is:** The total length of the entire chromosome this sequence comes from.

            **Why it exists — and why Start + Size alone isn't enough:**

            For **forward strand** (`+`) sequences, Start + Size is indeed sufficient to locate the region.
            But for **reverse strand** (`-`) sequences, the Start position is counted from the *end* of the
            chromosome. To convert that back to a normal forward-strand coordinate, you need to know the
            total chromosome length.
            """)
            st.info("""
            **Example:** If Source Size = 1,000,000 and a `-` strand sequence has Start = 25710 and Size = 1953,
            the forward-strand coordinates are: 1,000,000 − 25710 − 1953 = **972,337** to **974,290**.

            Without the Source Size field, this calculation would be impossible.
            """)

            st.markdown("---")
            st.markdown("##### Field 7: `GGCTACGACGA...` — The Aligned Sequence")
            st.markdown("""
            **What it is:** The actual DNA letters for this species in this aligned region.

            **Why it exists:** This is the heart of the data — the actual genetic information.

            **About the gaps (`-`):** You'll see `-` characters scattered through the sequence. These are
            **alignment gaps** — positions where the alignment algorithm determined that one species has DNA
            that another species doesn't (an insertion or deletion event in evolution). The gaps make all
            species' sequences the same length within a block, so you can compare them position by position.
            """)

        with st.expander("4️⃣ The Species.Chromosome Naming Convention", expanded=False):
            st.markdown("""
            The **Source** field (Field 2) always follows the format **`species.chromosome`**.
            This is how every tool — including MAFin — identifies whose DNA is whose.
            """)
            naming_df = pd.DataFrame([
                ["ref_genome.chr24", "ref_genome", "chr24", "Chromosome 24 of the reference genome"],
                ["fake_species_3.chr24", "fake_species_3", "chr24", "Chromosome 24 of a different species"],
                ["fake_species_5.chrX", "fake_species_5", "chrX", "Chromosome X of another species"],
            ], columns=["Full Identifier", "Species (before dot)", "Chromosome (after dot)", "Interpretation"])
            st.dataframe(naming_df, use_container_width=True, hide_index=True)
            st.warning("""
            **Key insight:** `ref_genome.chr24` and `fake_species_3.chr24` both say `chr24`, but they belong
            to **completely different genomes**. They are never confused — the species prefix always disambiguates them.
            """)

        with st.expander("5️⃣ Multiple Blocks for the Same Chromosome", expanded=False):
            st.markdown("""
            A single chromosome often appears in **many different alignment blocks**. This is completely normal!

            **Why?** Chromosomes are millions of letters long, but only *some* regions are similar enough between
            species to be aligned. Each aligned region becomes its own block, covering a different genomic window.
            """)
            st.code("""# Block A — covers positions 17744–20201 on chr24
a score=21475.0
s ref_genome.chr24   17744 2457 + 1000000   GGCTACGACGA...

# Block B — covers positions 18225–20813 on chr24
a score=18945.3
s ref_genome.chr24   18225 2588 + 1000000   TCTGGTAGACC...

# Block C — covers positions 18716–19129 on chr24
a score=15022.1
s ref_genome.chr24   18716  413 + 1000000   AAAGCTTTTGG...""", language="text")
            st.markdown("""
            MAFin processes **each block independently**. If your motif is found in Block A at position 19318
            and also in Block B at position 19540, both hits appear as separate rows in the output — each
            with their own conservation scores calculated from the species aligned in *that specific block*.
            """)

        with st.expander("6️⃣ Why Output Only Shows Reference Chromosomes", expanded=False):
            st.markdown("""
            When you look at MAFin's output, you'll see chromosome names like `chr24` and `chr29` —
            but **never** `fake_species_3.chr2` or `fake_species_5.chrX`. Why?

            **MAFin reports all hits relative to the reference genome** (the first `s` line in each block).
            The other species are not ignored — their sequences are used to calculate the **Conservation Score**,
            which tells you how well the motif is preserved across evolution.
            """)
            ref_output_df = pd.DataFrame([
                ["chr24 at position 19318", "The motif was found on the reference genome's chromosome 24"],
                ["Score = 100%", "The motif is perfectly conserved in all other aligned species at that position"],
                ["Score = 85%", "85% of letters matched across the other species — some variation exists"],
            ], columns=["What you see in output", "What it means"])
            st.dataframe(ref_output_df, use_container_width=True, hide_index=True)
            st.markdown("This design keeps the output clean: **one row per reference hit**, with conservation summarizing all species.")

    with tab1:
        st.subheader("🔍 Scenario 1: Find Exact DNA Sequences (K-mers)")
        
        st.markdown("""
        ### What is a K-mer Search?
        A **K-mer** is a sequence of exactly K nucleotides (the "K" is the length — a 16-letter sequence is a "16-mer"). This search method finds **exact matches** with no variations allowed.
        
        **Use this when:**
        - You have a specific DNA sequence you know and want to find it exactly
        - You're looking for primer binding sites (short sequences used in lab experiments to copy specific DNA regions)
        - You're searching for confirmed regulatory elements (DNA sequences that control when genes turn on or off)
        - Speed is important (fastest search method)
        
        *See the **📚 Glossary** for detailed definitions of these terms.*
        """)
        
        st.info("💡 **Real-world example:** A cancer researcher discovers a 16-letter-long mutation in a patient's tumor (e.g., `ATGCATGCATGCATGC`). They search a MAF file with human + chimp + mouse genomes to find where this exact sequence appears. If it appears identically in all species at the same location, it's likely an important conserved region. If only in humans, it's a recent human-specific mutation.")
        
        st.markdown("### 🔄 Step-by-Step Execution Workflow")
        st.info("Follow along as we see exactly how the inputs transform into outputs.")

        with st.expander("1️⃣ STEP 1: Provide Your Search Sequences", expanded=True):
            st.markdown("""
            **What is the input data?**
            A list of exact DNA sequences (K-mers), provided with one sequence per line. In this example, these are 16-mers (sequences of exactly 16 nucleotides).
            
            **Why did we prepare it this way?**
            We use exact sequences when looking for a specific, known DNA string with absolutely no variations allowed (like a known mutation or a confirmed regulatory element).
            
            **How does the application use it?**
            MAFin reads each sequence from this list and slides through every alignment block, checking if any sequence window matches your query exactly, letter for letter.
            """)
            st.code("""AAAATTTTGGGGCCCC
TTTTTTTTTTTTTTTT
GGGGGGGGGGGGGGGG""", language="text")

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)
        
        with st.expander("2️⃣ STEP 2: Provide the MAF File (Aligned Genomes)", expanded=False):
            st.markdown("Next, provide the alignment file that MAFin will search through. This contains sequences from multiple species aligned together (you can find a MAF file at the **[UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html)**)")
            st.markdown(highlight_maf(maf_snippet, ["AAAATTTTGGGGCCCC", "TTTTTTTTTTTTTTTT"]), unsafe_allow_html=True)
            st.caption("Contains aligned sequences from reference genome + multiple other species. *(Note: This snippet is just a tiny window into the chromosome. A real MAF file has thousands of these blocks!)*")
            
            st.info("💡 **New to MAF files?** Check out the **📂 Reading MAF Files** tab for a full breakdown of what each field means and how species identifiers work.")

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)
        
        with st.expander("3️⃣ STEP 3: MAFin Finds the Matches", expanded=False):
            st.markdown("""
            Behind the scenes, MAFin maps your input sequences to the genomic coordinates:
            1. It reads each query sequence from your K-mer list (e.g., `AAAATTTTGGGGCCCC` — a sequence 16 letters long).
            2. It slides through every alignment block, checking if any 16-character window matches your query exactly, letter for letter.
            
            *As MAFin scans through all the thousands of blocks in the full MAF file, it finds multiple hits across different coordinates and chromosomes:*
            """)
            
            explain_df = pd.DataFrame([
                ["AAAATTTTGGGGCCCC", "chr29:25671-25687", "+", "Exact forward-strand hit"],
                ["AAAATTTTGGGGCCCC", "chr29:24945-24961", "+", "Second exact forward-strand hit"],
                ["TTTTTTTTTTTTTTTT", "chr24:19318-19334", "+", "Exact forward-strand hit"],
                ["TTTTTTTTTTTTTTTT", "chr29:26421-26437", "+", "Second exact forward-strand hit"]
            ], columns=["K-mer", "BED Interval", "Strand", "Interpretation"])
            st.dataframe(explain_df, use_container_width=True, hide_index=True)

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)

        with st.expander("4️⃣ STEP 4: Output BED Format Results", expanded=False):
            st.markdown("Finally, every exact hit is reported as one row in a BED file with its coordinates and conservation score.")
            
            df1 = pd.DataFrame([
                ["chr29", 25671, 25687, "AAAATTTTGGGGCCCC", 100.0, "Forward (+)"],
                ["chr29", 24945, 24961, "AAAATTTTGGGGCCCC", 100.0, "Forward (+)"],
                ["chr24", 19318, 19334, "TTTTTTTTTTTTTTTT", 90.18, "Forward (+)"],
                ["chr29", 26421, 26437, "TTTTTTTTTTTTTTTT", 85.71, "Forward (+)"]
            ], columns=["Chromosome", "Start", "End", "Motif", "Score", "Strand"])
            st.table(df1)
            
            with st.expander("📚 What each column means:", expanded=False):
                output_data = {
                    "Column": ["Chromosome", "Start", "End", "Motif", "Score", "Strand"],
                    "Meaning": [
                        "Which chromosome/sequence was hit",
                        "Position where match begins (0-based, meaning the first letter on a chromosome is position 0)",
                        "Position where match ends (non-inclusive: End = Start + sequence length, pointing to the position *after* the last matched letter)",
                        "The exact sequence that was matched",
                        "Conservation % (100% = found in all genomes)",
                        "+ (forward) or - (reverse complement)"
                    ]
                }
                st.dataframe(pd.DataFrame(output_data), use_container_width=True, hide_index=True)
                
            with st.expander("🔬 Exact rows from the demo files", expanded=False):
                st.markdown("**Exact BED rows (directly from `test_16mer_motif_hits.bed`):**")
                st.code(
                    """chr24\t19318\t19334\tTTTTTTTTTTTTTTTT\t90.18\t+
chr29\t25671\t25687\tAAAATTTTGGGGCCCC\t100.0\t+
chr29\t26421\t26437\tTTTTTTTTTTTTTTTT\t85.71\t+
chr29\t24945\t24961\tAAAATTTTGGGGCCCC\t100.0\t+""",
                    language="text",
                )
                st.caption("Each query sequence is 16 characters long, so End − Start = 16 in every row. The End coordinate is non-inclusive (it points to the position *after* the last matched letter). If a query K-mer has no exact match anywhere in the reference (for example `GGGGGGGGGGGGGGGG` in the bundled demo), it simply does not appear in the output.")
        
        st.markdown("---")
        st.markdown("### Understanding the Output")
        st.markdown("""
        - **Score = 100%**: Pattern found identically in all aligned genomes → highly conserved (likely functionally important — nature preserved it across evolution)
        - **Score < 100%**: Pattern found with variations in some genomes → partially conserved (some evolutionary drift occurred)
        - **Strand**: Indicates which DNA strand the hit is on — use the genome browser with these coordinates to visualize
        """)
        
        st.info("💡 **Why only reference chromosomes?** Output lists `chr29` and `chr24` (reference genome) — not `fake_species_3.chr2`. Other species are folded into the **Score**. See the **📂 Reading MAF Files** and **❓ Common Questions** tabs for details.")
        

    
    with tab2:
        st.subheader("🧬 Scenario 2: Find Patterns with Variations (Regular Expressions)")
        
        st.markdown("""
        ### What is a Regular Expression (Regex)?
        A **regex** is a pattern with **allowed flexibility** at certain positions. This lets you capture biological variants — because evolution introduces small changes, the same functional motif might differ slightly between species.
        
        **Use this when:**
        - You want to allow certain positions to vary (e.g., using `.` for any nucleotide)
        - You're looking for motifs with known variable positions
        - You want to be more flexible than exact matching but don't have a PWM
        - You're searching for degenerate motifs (variants of the same functional sequence)
        
        ### Common Regex Patterns
        """)
        
        regex_help = pd.DataFrame([
            ["[GT]", "Matches either G or T"],
            ["[ACG]", "Matches any of A, C, or G"],
            [".", "Matches any single nucleotide"],
            [".{3}", "Matches exactly 3 of any nucleotide"],
            [".+", "Matches 1 or more of any nucleotide"],
            ["A", "Matches the literal letter 'A'"],
            ["A[GT]T", "Matches: AGT or ATT"],
            ["CC.{3}AT", "Matches: CC + any 3 letters + AT"],
        ], columns=["Pattern", "What it Matches"])
        st.dataframe(regex_help, use_container_width=True, hide_index=True)
        st.caption("💡 **When to use which:** Use `[XY]` when you know *which* specific letters are allowed. Use `.` when *any* letter is acceptable. Use `.{N}` for variable-length spacers between conserved elements.")
        
        st.info("💡 **Real-world example:** An immunologist knows immune regulatory elements have a conserved `TTCC` core but flexible flanking regions. They use regex pattern `[AC]TTCC[GT]` to find all variants (`ATTCCG`, `ATTCCT`, `CTTCCG`, etc.) that would be missed by exact k-mer searching.")
        
        st.markdown("### 🔄 Step-by-Step Execution Workflow")
        st.info("Follow along as we see exactly how the inputs transform into outputs.")

        with st.expander("1️⃣ STEP 1: Define Your Flexible Patterns", expanded=True):
            st.markdown("""
            **What is the input data?**
            A list of regular expression (regex) patterns, where special characters define flexibility.
            
            **Why did we prepare it this way?**
            Biology is messy, and the same functional motif might differ slightly between species. Instead of listing every possible variant, regex allows us to capture them all by specifying which positions can vary (e.g., `[GT]` means the position can be G or T).
            
            **How does the application use it?**
            MAFin evaluates every sequence chunk in the MAF file against these rules to see if the sequence fits the allowed flexibility of your pattern.
            """)
            st.code("""A[GT]TCCG[AT]A
CC.{3}AT.{3}GC
GG[TC]ATAT[GA]""", language="text")
            st.caption("💡 **Tip:** Characters inside `[ ]` are alternatives (one character matches one DNA letter). `.{3}` means 'any 3 letters'. So `A[GT]TCCG[AT]A` matches sequences that are exactly **8 letters long** — the brackets don't add length, they add *flexibility* at that position.")

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)
        
        with st.expander("2️⃣ STEP 2: Provide the MAF File", expanded=False):
            st.markdown("Next, provide the same alignment file:")
            st.markdown(highlight_maf(maf_snippet, ["AGTCCGAA", "ATTCCGTA", "CCTATATGACGC", "CCCATATATAGC", "GGTATATG", "GGCATATA"]), unsafe_allow_html=True)
            st.caption("Contains aligned sequences from reference genome + multiple other species")

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)
        
        with st.expander("3️⃣ STEP 3: MAFin Evaluates Windows Against Rules", expanded=False):
            st.markdown("""
            Behind the scenes, MAFin:
            1. Slides across the reference sequence in the MAF file.
            2. Evaluates every sequence chunk to see if it matches the flexibility allowed by your pattern.
            
            **How pattern length maps to match length:** A pattern like `A[GT]TCCG[AT]A` has 14 visible characters in the pattern text, but `[GT]` and `[AT]` each represent *one* position in the DNA. So the pattern matches sequences that are exactly **8 letters long**.
            """)
            
            pattern_match_df = pd.DataFrame([
                ["A[GT]TCCG[AT]A", "AGTCCGAA", "8 letters", "G matches [GT] and A matches [AT]"],
                ["A[GT]TCCG[AT]A", "ATTCCGTA", "8 letters", "T matches [GT] and T matches [AT]"],
                ["CC.{3}AT.{3}GC", "CCTATATGACGC", "12 letters", ".{3} matches 'TAT' and 'GAC'"],
                ["CC.{3}AT.{3}GC", "CCCATATATAGC", "12 letters", ".{3} matches 'CAT' and 'ATA'"],
                ["GG[TC]ATAT[GA]", "GGTATATG", "8 letters", "T matches [TC] and G matches [GA]"],
                ["GG[TC]ATAT[GA]", "GGCATATA", "8 letters", "C matches [TC] and A matches [GA]"],
            ], columns=["Pattern", "Matched Sequence", "Match Length", "Why it matches"])
            st.dataframe(pattern_match_df, use_container_width=True, hide_index=True)

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)

        with st.expander("4️⃣ STEP 4: Output All Matching Patterns", expanded=False):
            st.markdown("Every sequence that fits the rules is captured and reported in the BED output. These coordinates come directly from the sample MAF data shown in Step 2:")
            
            df2 = pd.DataFrame([
                ["chr29", 25699, 25707, "AGTCCGAA", 100.0, "Forward (+)"],
                ["chr29", 25870, 25882, "CCTATATGACGC", 100.0, "Forward (+)"],
                ["chr29", 25976, 25984, "GGTATATG", 100.0, "Forward (+)"],
                ["chr29", 24973, 24981, "ATTCCGTA", 100.0, "Forward (+)"],
                ["chr29", 25144, 25156, "CCCATATATAGC", 100.0, "Forward (+)"],
                ["chr24", 19540, 19548, "GGCATATA", 90.18, "Forward (+)"],
            ], columns=["Chromosome", "Start", "End", "Motif Found", "Score", "Strand"])
            st.table(df2)
            st.caption("Each row maps directly back to a highlighted match in the MAF snippet above. Notice how different blocks on different chromosomes produce separate output rows.")
        
        st.markdown("---")
        st.markdown("### Key Differences from K-mers")
        st.markdown("""
        | Feature | K-mers (Exact) | Regex (Flexible) |
        |---------|---|---|
        | Matches | Only exact sequences | All sequences matching the pattern |
        | Speed | Fastest (simple letter comparison) | Slightly slower (rule evaluation at each position) |
        | Sensitivity | Low (might miss variations) | High (catches variants) |
        | Use Case | Known fixed sequences | Motifs with known variations |
        """)
        

    
    with tab3:
        st.subheader("📊 Scenario 3: Search Using Motif Profiles (PWM)")
        
        st.markdown("""
        ### What is a Position Weight Matrix (PWM)?
        A **PWM** (or Position Frequency Matrix) is a **probabilistic model** of a motif. Instead of specifying exact sequences 
        or patterns, you provide the probability of each nucleotide at each position.
        
        **Why use probabilities instead of exact letters or regex?** Real transcription factor binding sites don't follow 
        rigid rules — they have *preferences*. A PWM captures this biological reality: position 1 might strongly prefer G (40%) 
        but tolerate C (30%). Neither K-mers (too rigid) nor regex (only yes/no per position) can express 
        "I prefer G but will accept C with lower confidence."
        
        **Use this when:**
        - You're searching for a well-characterized biological motif (like transcription factor binding sites)
        - You have prior knowledge of which bases are preferred at each position
        - You want scoring based on information content (how confident the match is)
        - You're working with JASPAR or similar motif databases
        """)
        
        st.info("💡 **Real-world example:** A developmental biologist downloads a TP53 transcription factor PWM from JASPAR (a public database of experimentally-determined binding profiles). They search aligned mammal genomes to find all potential p53 binding sites and identify which are conserved across species (likely functional for tumor suppression).")
        
        st.markdown("### How PWM Scoring Works")
        st.markdown("**Why multiply probabilities?** For each candidate sequence, MAFin looks up the probability of each observed letter at each position in the PWM, then multiplies them together. This gives the overall likelihood that the motif model would produce this exact sequence. Higher likelihood = better match.")
        
        pwm_example_data = {
            "Position": ["Pos 1", "Pos 2", "Pos 3", "Pos 4"],
            "A": [0.2, 0.1, 0.4, 0.3],
            "C": [0.3, 0.5, 0.1, 0.2],
            "G": [0.4, 0.3, 0.4, 0.5],
            "T": [0.1, 0.1, 0.1, 0.0],
        }
        st.dataframe(pd.DataFrame(pwm_example_data), use_container_width=True, hide_index=True)
        
        st.markdown("""
        **Interpretation:**
        - Position 1: G is most likely (40%), C is second (30%)
        - Position 2: C is strongly preferred (50%)
        - Position 3: A and G are equally likely (40% each)
        - Position 4: G is most likely (50%), T is impossible (0%)
        
        **What this means:** A sequence like **GCAG** would score high, while **TCAT** would score very low.
        """)
        
        # STEP-BY-STEP WORKFLOW
        st.markdown("### 🔄 Step-by-Step Execution Workflow")
        st.info("Follow along as we see exactly how the inputs transform into outputs.")

        with st.expander("1️⃣ STEP 1: Provide PWM Profile (JASPAR Format)", expanded=True):
            st.markdown("""
            **What is the input data?**
            A mathematical probability matrix (Position Weight Matrix or PWM) in JASPAR format, showing the probability of each nucleotide (A, C, G, T) occurring at each position.
            
            **Why did we prepare it this way?**
            Real biological binding sites have *preferences*, not rigid rules. A PWM captures this reality, allowing us to express complex rules like "Position 1 strongly prefers G (40%) but will tolerate C (30%)."
            
            **How does the application use it?**
            MAFin scores every candidate sequence window by multiplying the probability of each observed base at its position. It takes the logarithm to produce a final Information Content Score. Matches scoring above your statistical p-value threshold are reported.
            """)
            st.code(""">TFBS_motif
A  [ 0.2  0.1  0.4  0.3 ]
C  [ 0.3  0.5  0.1  0.2 ]
G  [ 0.4  0.3  0.4  0.5 ]
T  [ 0.1  0.1  0.1  0.0 ]""", language="text")
            st.caption("JASPAR format: probabilities for each nucleotide at each position")

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)
        
        with st.expander("2️⃣ STEP 2: Provide the MAF File", expanded=False):
            st.markdown("Next, provide the same alignment file:")
            st.markdown(highlight_maf(maf_snippet, ["GCAG", "GCGA"]), unsafe_allow_html=True)

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)
        
        with st.expander("3️⃣ STEP 3: MAFin Scores All Sequence Windows", expanded=False):
            st.markdown("""
            Behind the scenes, MAFin slides across the alignment and scores every 4-letter candidate window against the PWM profile.
            
            **How scores are calculated:** For each window, MAFin multiplies the probability in each column of the PWM for the observed base, then takes the logarithm. Why logarithm? It converts tiny probability products (like 0.00032) into human-readable scores (like 8.45). Without the log, you'd be comparing unwieldy decimal numbers. Higher scores mean the sequence is a better match to the motif model.
            """)
            
            scoring_df = pd.DataFrame([
                ["GCAG", "G(0.4)×C(0.5)×A(0.4)×G(0.5)", "8.45", "✅ Above threshold → reported"],
                ["AAAT", "A(0.2)×A(0.1)×A(0.4)×T(0.0)", "−∞", "❌ T is impossible at Pos 4 → rejected"],
                ["GCGA", "G(0.4)×C(0.5)×G(0.4)×A(0.3)", "7.92", "✅ Above threshold → reported"],
                ["ATGA", "A(0.2)×T(0.1)×G(0.4)×A(0.3)", "3.10", "❌ Below threshold → ignored"],
            ], columns=["Window", "Per-Position Probabilities", "Score", "Verdict"])
            st.dataframe(scoring_df, use_container_width=True, hide_index=True)
            st.caption("Probabilities are read directly from the PWM table in Step 1. A 0% probability at any position makes the score −∞ (impossible match).")

        st.markdown("<h2 style='text-align: center; color: #64748b; margin-top: 10px; margin-bottom: 10px;'>⬇️</h2>", unsafe_allow_html=True)

        with st.expander("4️⃣ STEP 4: Output High-Scoring Matches", expanded=False):
            st.markdown("Only the windows that score above the p-value threshold are saved to the output BED file. Unlike K-mer and Regex results, the **Motif** column shows the JASPAR matrix ID (not the matched sequence), because PWM matches are probabilistic rather than exact.")
            
            df3 = pd.DataFrame([
                ["chr29", 25692, 25696, "MA0001.3", 8.45, "Forward (+)"],
                ["chr29", 25768, 25772, "MA0001.3", 7.92, "Forward (+)"],
                ["chr29", 24966, 24970, "MA0001.3", 8.45, "Forward (+)"],
                ["chr24", 19339, 19343, "MA0001.3", 7.92, "Forward (+)"],
            ], columns=["Chromosome", "Start", "End", "Motif", "Score", "Strand"])
            st.table(df3)
            st.caption("Score = information content score (higher = better match to PWM profile). Coordinates are illustrative using the same sample MAF data.")
        
        st.markdown("---")
        st.markdown("### Understanding PWM Output")
        st.markdown("""
        - **Score > 8.0**: Excellent match to the motif profile
        - **Score 6.0-8.0**: Good match
        - **Score < 6.0**: Weak match (you can filter these with p-value thresholds)
        - **P-value**: MAFin can assign statistical significance to each match
        
        **Note:** These score ranges are general guidelines. The actual meaningful threshold depends on your specific motif and genome. The **p-value filter** (set in Step 3 of the Run page) is the statistically rigorous way to decide what counts as a real match — it accounts for genome composition and motif length.
        """)
        

    
    with tab_faq:
        st.subheader("❓ Common Questions")
        st.markdown("Answers to the most frequently asked questions about MAFin, organized by topic.")

        st.markdown("### 📄 Understanding the Input")

        with st.expander("Why does the same chromosome appear in multiple alignment blocks?", expanded=False):
            st.markdown("""
            This is completely normal. A MAF file breaks the genome into many alignment blocks, and each block covers
            a **different region** of the chromosome. Think of it like chapters in a book — chromosome 24 is the book,
            and each block is a different chapter.

            MAFin scans **every block independently**, so a motif found at position 19318 in one block and position
            19540 in another block will both appear as separate rows in your output.
            """)

        with st.expander("What is the difference between 'Reference Sequence' and 'All Aligned Sequences'?", expanded=False):
            st.markdown("""
            This setting controls **where MAFin looks** for your motif:

            - **Reference Sequence** (default): MAFin only searches the **first sequence** in each alignment block.
              This is the primary genome (e.g., the Human genome). It is faster because it scans fewer sequences.
            - **All Aligned Sequences**: MAFin searches **every species** in the block. Use this if your motif might
              only exist in a non-reference species. This is slower but more exhaustive.

            In both cases, the conservation score is still calculated by comparing against all other species.
            """)

        with st.expander("What does searching 'Both Strands (+/-)' do?", expanded=False):
            st.markdown("""
            DNA is **double-stranded**. Each strand runs in the opposite direction, and every letter has a complement:
            A↔T and C↔G. Your motif might sit on either strand.

            - **Forward Strand Only (+)**: MAFin only searches the sequence as written in the MAF file.
            - **Both Strands (+/-)**: MAFin also generates the **reverse complement** of your motif and searches for
              that too. For example, searching for `ATGC` on both strands also searches for `GCAT`.

            If you're unsure which strand your motif is on, choose **Both Strands**.
            """)

        with st.expander("What happens with gaps (-) in the alignment?", expanded=False):
            st.markdown("""
            Gaps (`-`) are inserted by alignment tools to keep comparable positions lined up across species. MAFin
            handles them as follows:

            - **During search**: Gaps are **stripped out** before searching. MAFin scans the continuous DNA letters only.
            - **During conservation scoring**: Gaps matter. If the reference has a letter but another species has a gap
              at that position, it counts as a **mismatch** (lowers the conservation score). If both have a gap, that
              position is simply skipped.
            """)

        st.markdown("### 📊 Understanding the Output")

        with st.expander("Why does my output say 'chr24' when another species also has 'chr24'?", expanded=False):
            st.markdown("""
            MAFin reports hits relative to the **reference genome only**. When the output says `chr24`, it always means
            the reference genome's chromosome 24 — never another species' chromosome.

            Internally, MAFin uses the full `species.chromosome` identifier (e.g., `ref_genome.chr24` vs `fake_species_3.chr24`)
            to keep them apart. The species prefix is stripped from the output because BED files follow the convention of using
            the reference genome's coordinate system.
            """)

        with st.expander("What does a Score of 85% actually mean?", expanded=False):
            st.markdown("""
            The Score is the **average conservation percentage** across all other species aligned at that position.

            - **100%** → Every letter in the motif is identical across all aligned species. Highly conserved!
            - **85%** → On average, 85% of the letters matched across other species. Some variation exists.
            - **50%** → Only half the letters matched. The motif region is divergent across species.

            A high score suggests the region is under evolutionary pressure to stay the same — potentially
            because it has an important biological function.
            """)

        with st.expander("What is the Similarity Vector in the CSV output?", expanded=False):
            st.markdown("""
            The Similarity Vector is a character-by-character comparison between the reference genome and one other species.
            Each position in the vector tells you exactly what happened at that letter:

            - `1` = **Match** — both genomes have the same letter
            - `0` = **Mismatch** — the genomes have different letters
            - `-` = **Gap** — one genome has a gap at that position

            **Example:** A vector of `111101100` means the first four letters matched, the fifth mismatched,
            the sixth and seventh matched, and the eighth and ninth mismatched.
            """)

        with st.expander("What is the difference between BED, CSV, and JSON output files?", expanded=False):
            st.markdown("""
            MAFin produces three output formats. Each contains the same results but in a different structure:

            | Format | Best For | What's Inside |
            |--------|----------|---------------|
            | **BED** | Genome browsers (UCSC, IGV) | Coordinates only: chromosome, start, end, motif, score, strand |
            | **CSV** | Spreadsheets (Excel, Google Sheets) | Full detail: coordinates + similarity vectors for each species |
            | **JSON** | Scripts and data pipelines (Python, R) | Hierarchical structured data with all fields |

            If you just want to know *where* your motifs are, use **BED**. If you need per-species detail, use **CSV**.
            If you're writing code to process the results, use **JSON**.
            """)

        st.markdown("### 🔍 Troubleshooting")

        with st.expander("Why is my K-mer not appearing in the output?", expanded=False):
            st.markdown("""
            If a K-mer doesn't appear in the output, it means **no exact match was found** in any reference
            sequence block. Common reasons:

            - The sequence may exist in a non-reference species but not in the reference (and you searched in "Reference" mode)
            - The sequence may span a gap (`-`) in the alignment, which breaks the exact match
            - There may be a single-letter difference that prevents an exact match — try using **Regex** mode with flexible positions instead
            """)

        with st.expander("Why do I get different results when I search 'Both Strands' vs 'Forward Only'?", expanded=False):
            st.markdown("""
            When you search both strands, MAFin finds matches on the **reverse complement** strand that would
            be invisible to a forward-only search. These additional hits will have a strand value of `-` in the output.

            This can sometimes **double** (or more) the number of results — this is expected behavior, not an error.
            The motif genuinely exists on the opposite strand at those locations.
            """)

        with st.expander("My PWM search returns too many / too few results. How do I tune it?", expanded=False):
            st.markdown("""
            The **p-value threshold** controls how strict your PWM search is:

            - **Lower p-value** (e.g., `1e-5`) → Stricter. Only the very best matches are reported. Fewer results.
            - **Higher p-value** (e.g., `1e-3`) → More permissive. Weaker matches are also included. More results.

            If you're getting too many hits, lower the p-value. If you're getting too few, raise it.
            The **background frequencies** also matter — if your genome is AT-rich, adjust them from the default
            `0.25 0.25 0.25 0.25` to reflect the actual composition.
            """)

    st.markdown("---")

    st.markdown("""
    ## 🎯 Next Steps

    Ready to start?
    1. **Prepare your MAF file** (aligned sequences from multiple genomes)
    2. **Choose a search method** based on your biological question
    3. **Head to "🚀 Run Motif Search"** to execute your analysis
    """)

elif nav_selection == "📚 Glossary":
    st.title("📚 Glossary of Terms")
    st.markdown("Each term below includes a short explanation and a concrete example. **Click on any term to expand it.**")
    
    search_term = st.text_input("🔍 Search Glossary", placeholder="e.g. Genome, BED file, PWM...", help="Start typing to instantly filter the glossary terms.").strip().lower()
    
    glossary_terms = {
        "Genome ID": "- **What it is:** a short identifier representing the species or genome assembly in the MAF file.\n- **Why it matters:** It's how MAFin (and all bioinformatics tools) know which species a sequence belongs to. Without it, you'd have raw DNA with no way to tell human from mouse.\n- **Example:** In the sequence name `hg38.chr1`, `hg38` is the Genome ID (Human genome assembly 38) and `chr1` is the chromosome.",
        "Genome": "- **What it is:** the complete set of DNA instructions found in a cell. In MAFin, it's the entire sequence of DNA for a specific species (like Human or Mouse).\n- **Why it matters:** Genomes are the raw material MAFin searches through. Each species has its own genome, and comparing them reveals what's shared (conserved) vs. what's unique.\n- **Example:** The entire Human genome has over 3 billion letters.",
        "Genome Block (Alignment Block)": "- **What it is:** a specific, small section of the genome where the DNA of multiple different species has been aligned together because they share an evolutionary history.\n- **Why it matters:** MAFin processes each block independently. Your motif might be found in multiple blocks at different positions, each with its own conservation score.\n- **Example:** A MAF file is essentially just a large list of these individual genome blocks.",
        "Alignment": "- **What it is:** putting sequences on top of each other so matching positions line up.\n- **Why it matters:** Alignment is the foundation of comparative genomics — it reveals which parts of DNA are shared across species, pointing to regions that are evolutionarily important.\n- **Example:**\n  ```text\n  Seq A: A C G T A\n  Seq B: A C - T A\n  ```\n  Here, `-` is a gap inserted so comparable positions align.",
        "Multiple alignment": "- **What it is:** the same idea as alignment, but for 3 or more sequences at once.\n- **Why it matters:** Comparing many species at once gives much stronger evidence for conservation than comparing just two. If a motif is identical in human, mouse, dog, AND chicken, it's almost certainly functionally important.\n- **Example:**\n  ```text\n  Human:  A C G T A\n  Mouse:  A C - T A\n  Dog:    A C G T -\n  ```\n  A MAF file stores many such multi-sequence alignment blocks.",
        "MAF (Multiple Alignment Format)": "- **What it is:** a text format that stores multiple-alignment blocks and genome coordinates.\n- **Why this format:** MAF was designed specifically for multi-species alignments because it stores not just the DNA letters, but also the coordinate mapping between each species' genome — information that simpler formats like FASTA don't carry.\n- **Example concept:** one block may represent a conserved region across human, mouse, and dog.",
        "Reference sequence": "- **What it is:** the primary genome that the other sequences in an alignment block are compared against, usually the first sequence in a MAF block.\n- **Why it matters:** MAFin reports all results relative to the reference genome's coordinates. The other species are used to calculate conservation, but the *locations* in your output always refer to the reference.\n- **Example:** in an alignment of Human, Mouse, and Dog, the Human genome is typically the reference. Searching 'reference' only looks at the Human DNA.",
        "BED file": "- **What it is:** a standard text format used by bioinformaticians to define specific regions on a genome.\n- **Why it matters:** BED is the universal language for 'this region on this chromosome.' Genome browsers (UCSC, IGV) expect this format, so MAFin's BED output can be directly loaded for visualization.\n- **Example:** it always contains at least three columns (Chromosome, Start position, End position), and MAFin adds the Motif name, Score, and Strand.",
        "DNA": "- **What it is:** a sequence made of four chemical letters: A (Adenine), C (Cytosine), G (Guanine), and T (Thymine). These are the building blocks of all genetic information.\n- **Why it matters:** DNA is what MAFin searches through. Every search — K-mer, Regex, or PWM — is ultimately looking for specific patterns in these four letters.\n- **Example:** `ATGCCGTA`",
        "RNA": "- **What it is:** similar to DNA, but uses U (Uracil) instead of T (Thymine). RNA is the working copy cells make from DNA.\n- **Why it matters:** MAFin works with DNA sequences (A, C, G, T). If your data contains RNA notation (U instead of T), convert it before searching.\n- **Example:** `AUGCCGUA`",
        "K-MER (exact sequence search)": "- **What it is:** an exact DNA string of length K. The name literally means 'a string of K letters.'\n- **Why it matters:** K-mer search is the fastest method because it only requires exact letter-by-letter comparison with no rule evaluation. Use it when you know the precise sequence.\n- **Example:** `ATGC` is a 4-mer; `AAAATTTTGGGGCCCC` is a 16-mer.\n- **Use when:** you want exact matches only.",
        "Regex (pattern search)": "- **What it is:** a rule-based pattern language that allows flexibility at specific positions.\n- **Why it matters:** In biology, the same functional motif often has slight variations between species. Instead of listing every variant as a separate K-mer, a single regex captures them all.\n- **Example:** `A[GT]T` matches `AGT` or `ATT`.\n- **Use when:** some positions can vary.",
        "PWM / JASPAR (motif profile search)": "- **What it is:** a matrix-based model of base preferences at each motif position. JASPAR is a public database of experimentally-determined transcription factor binding profiles.\n- **Why it matters:** Real biological binding sites don't follow rigid rules — they have *preferences*. A PWM captures this reality: position 1 might strongly prefer G (40%%) but tolerate C (30%%). Neither K-mers nor regex can express this kind of graded preference.\n- **Example idea:** position 1 may strongly prefer `A`, position 2 may prefer `C/G`.\n- **Use when:** you want biologically similar matches, not only exact strings.",
        "P-value (confidence threshold)": "- **What it is:** a statistical measure of how likely a motif match occurred by pure random chance. A lower number means a stricter, higher-confidence search.\n- **Why it matters:** A high-scoring PWM match might still be random coincidence in a large genome. The p-value accounts for sequence composition and motif length to separate real biological signals from noise.\n- **Example:** a threshold of `1e-4` (0.0001) means there is only a 1 in 10,000 chance that the found sequence matches the motif randomly.",
        "Background frequencies": "- **What it is:** the natural, average occurrence rate of the letters A, C, G, and T in the specific genome you are searching.\n- **Why it matters:** If your genome is 60%% AT, then an AT-rich motif match is less surprising than in a GC-rich genome. The background frequencies adjust for this bias, ensuring the p-value is calibrated to your specific genome.\n- **Example:** `0.25 0.25 0.25 0.25` assumes all four bases are equally common. If a genome is AT-rich, it might look closer to `0.30 0.20 0.20 0.30`.",
        "Motif": "- **What it is:** a recurring biological pattern in DNA that may tolerate small variations.\n- **Why it matters:** Motifs mark functional regions — gene switches, protein binding sites, regulatory elements. Finding them reveals which parts of DNA actually *do* something rather than being passive.\n- **Example:** a transcription factor binding pattern shared across many locations in the genome.",
        "Reverse complement": "- **What it is:** the sequence you get by reading the opposite DNA strand. The bases swap (A↔T, C↔G) and the order reverses.\n- **Why it matters:** DNA is double-stranded, and genes can sit on either strand. A motif on the reverse strand is invisible to a forward-only search. Searching both strands ensures you don't miss biological signals.\n- **Example:** reverse complement of `ATGC` is `GCAT`.",
        "Conservation Percentage (Output)": "- **What it is:** a measure of how conserved a motif is across all aligned genomes.\n- **Why it matters:** Regions conserved across millions of years of evolution are likely functionally important — nature has preserved them because mutations there would be harmful. Low conservation suggests the region is free to mutate without consequence.\n- **Example:** `100%%` means the motif was found perfectly intact across every single species in the alignment block. `50%%` means only half the bases matched.",
        "Similarity Vector (Output)": "- **What it is:** a visual binary string showing exactly which bases matched or mismatched in a specific genome compared to the reference.\n- **Why it matters:** While the conservation percentage gives you a summary number, the similarity vector shows you *exactly where* the differences are — essential for understanding which specific positions are variable.\n- **Example:** `111101100` where `1` is a match, `0` is a mismatch, and `-` is a gap.",
        "Strand (Output)": "- **What it is:** denotes whether the motif was found on the forward (`+`) or reverse (`-`) DNA strand.\n- **Why it matters:** Knowing the strand is essential for biological interpretation. Many genes and regulatory elements are strand-specific — a motif on the wrong strand might not be biologically relevant.\n- **Example:** if you search both strands, hits on the opposite strand will be labeled `-`.",
        "Chromosome": "- **What it is:** a long, continuous thread of DNA. Genomes are broken up into multiple chromosomes.\n- **Why it matters:** Chromosomes are how genomes are organized. When MAFin reports a hit on 'chr24', it tells you exactly which physical DNA molecule the motif was found on.\n- **Example:** Humans have 23 pairs of chromosomes. In a MAF file, `chr1` refers to chromosome 1.",
        "0-indexed Coordinates": "- **What it is:** a counting system used in computer science (and BED files) where the first position is counted as `0` instead of `1`.\n- **Why this convention:** It makes length calculation trivial — length = End − Start (no +1 needed). This avoids the off-by-one errors that are common with 1-indexed systems.\n- **Example:** If a motif is at the very beginning of a chromosome, its start position is `0`.",
        "Information Content Score": "- **What it is:** a score used specifically for PWM searches that measures how well a sequence matches the mathematical probability model.\n- **Why it matters:** Unlike K-mer or Regex searches (which are yes/no), PWM matches are graded. The information content score tells you *how good* the match is, letting you prioritize the strongest candidates.\n- **Example:** A score of `8.45` means it's an excellent match to the motif profile, while a score below `6.0` might be a weak match.",
        "CSV Format (Output)": "- **What it is:** a spreadsheet-friendly format (Comma-Separated Values) best for opening in Excel or statistical software.\n- **Why it matters:** Unlike BED (which only has coordinates), CSV includes per-species similarity vectors and detailed scoring — everything you need for downstream biological analysis.\n- **Example:** Each match is a row, with columns for coordinates, motifs, and similarity vectors.",
        "JSON Format (Output)": "- **What it is:** a hierarchical, machine-readable format (JavaScript Object Notation) best for automated scripts and data pipelines.\n- **Why it matters:** If you're writing code to process results (Python, R), JSON's structured format is far easier to parse than tab-separated text files.\n- **Example:** Contains raw, structured data perfect for writing Python parsing scripts.",
        "Nucleotide": "- **What it is:** a single building block of DNA or RNA. There are four nucleotides in DNA: Adenine (A), Cytosine (C), Guanine (G), and Thymine (T).\n- **Why it matters:** Nucleotides are the fundamental unit of measurement in genomics. When we say a 'K-mer', we mean K nucleotides. When we say 'conservation', we're comparing nucleotide by nucleotide.\n- **Example:** The sequence `ATGC` is 4 nucleotides long. When we say a '16-mer', we mean a sequence of 16 nucleotides."
    }
    
    for term, definition in glossary_terms.items():
        if not search_term or search_term in term.lower() or search_term in definition.lower():
            with st.expander(term):
                st.markdown(definition)

elif nav_selection == "📖 Architecture":
    st.title("📖 How MAFin Works")
    st.markdown("Understanding the engine behind MAFin. Choose a topic below to learn more.")
    
    hw_tab1, hw_tab2, hw_tab3 = st.tabs(["🏗️ Architecture", "🧬 Conservation Scoring", "🔄 Reverse Complement"])
    
    with hw_tab1:
        st.subheader("System Architecture")
        with st.expander("🧠 Understanding the Diagram Details", expanded=True):
            st.markdown(
                "The diagram below shows the exact steps MAFin takes under the hood. Here is a breakdown of the terminology used:\n\n"
                "- **Split File into Chunks:** *Divide and conquer.* MAF files can be enormous (gigabytes of aligned genomes). Instead of processing the entire file sequentially, MAFin divides it into chunks at alignment block boundaries so multiple CPU cores can work simultaneously.\n"
                "- **Build Aho-Corasick Automaton / Compile Regex / Compute PWM Thresholds:** *Prepare the search engine.* Before scanning begins, MAFin prepares the appropriate search tool. For K-mers, it builds an Aho-Corasick automaton (provided by the `pyahocorasick` library) that can search for *thousands of sequences simultaneously* in a single pass — dramatically faster than searching one at a time. For Regex, it compiles the pattern rules. For PWMs, it calculates the statistical score threshold from your p-value by sampling random sequences.\n"
                "- **Launch Parallel Workers:** *Speed through parallelism.* Each chunk is assigned to a separate CPU core (using Python's `multiprocessing` module). All workers run the same search logic independently, writing their hits to temporary files.\n"
                "- **Search Target Sequence:** *Find the motif.* Within each alignment block, MAFin searches the target sequence (the reference genome by default, or a specific genome if you chose one). It strips out alignment gaps first, then scans the continuous DNA letters.\n"
                "- **Scan Aligned Species & Compute Conservation:** *Measure evolutionary preservation.* When a motif is found, MAFin checks the same aligned position in every other species. It builds a similarity vector (`1`=match, `0`=mismatch, `-`=gap) and calculates the average conservation percentage.\n"
                "- **Merge Results from All Workers:** *Reassemble the puzzle.* After all parallel workers finish, MAFin collects their temporary output files and merges them into unified, per-genome result sets.\n"
                "- **Generate Output Files:** *Deliver the answer.* MAFin always produces a BED file (genomic coordinates). If you requested a detailed report, it also generates CSV (spreadsheet-friendly) and JSON (script-friendly) files."
            )
            
        components.html("""
        <style>
            .edgeLabel, .edgeLabel span {
                background-color: #1e293b !important;
                color: #f8fafc !important;
            }
        </style>
        <div class="mermaid">
        graph TD
            A[Upload MAF File] --> B[Split File into Chunks]
            C[Define Search Pattern] --> D{Search Engine Type}
            
            D -->|Exact K-mers| E[Build Aho-Corasick Automaton]
            D -->|Regex| F[Compile Regex Patterns]
            D -->|JASPAR Matrix| G[Compute PWM Thresholds]
            
            B --> H[Launch Parallel Workers]
            E --> H
            F --> H
            G --> H
            
            H --> I[Parse Next Alignment Block]
            I --> J[Search Target Sequence in Block]
            J --> K{Motif Found?}
            K -->|No| N{More Blocks in Chunk?}
            K -->|Yes| L[Scan Aligned Species at Match Position]
            L --> M[Compute Conservation & Similarity Vector]
            M --> W[Write Hit to Temp File]
            W --> N
            N -->|Yes| I
            N -->|No| O[Merge Results from All Workers]
            O --> P[Generate Output Files]
            P --> Q[(BED + CSV + JSON Reports)]
            
            classDef default fill:#1e293b,stroke:#475569,stroke-width:2px,color:#f8fafc;
            classDef startNode fill:#4f46e5,stroke:#3730a3,stroke-width:2px,color:white;
            classDef engine fill:#9333ea,stroke:#7e22ce,stroke-width:2px,color:white;
            classDef parallel fill:#0891b2,stroke:#0e7490,stroke-width:2px,color:white;
            classDef endNode fill:#059669,stroke:#047857,stroke-width:2px,color:white;
            
            class A,C startNode;
            class E,F,G engine;
            class H parallel;
            class Q endNode;
        </div>
        <script type="module">
            import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
            mermaid.initialize({ 
                startOnLoad: true, 
                theme: 'dark',
                themeVariables: {
                    edgeLabelBackground: '#1e293b',
                    primaryTextColor: '#f8fafc'
                }
            });
        </script>
        """, height=900, scrolling=True)
        
    with hw_tab2:
        st.subheader("How we calculate Conservation")
        st.markdown(
            "Once MAFin finds your motif in the reference sequence, it checks all other aligned species to see if they share that same sequence. "
            "It does this by comparing the **ungapped sequences**, while strictly respecting the evolutionary gaps (`-`) defined in your MAF file. "
            "Why treat gaps specially? Because gaps represent evolutionary insertions or deletions (indels) — the DNA isn't *different*, it's simply *absent* in one species. Counting a gap as a 'mismatch' would be misleading."
        )
        
        st.markdown("**The Scoring Rules:**")
        col_r1, col_r2, col_r3, col_r4 = st.columns(4)
        col_r1.success("✅ **Match (1)**\nBases are identical.")
        col_r2.error("❌ **Mismatch (0)**\nBases are different.")
        col_r3.warning("➖ **Gap (-)**\nReference has a gap.")
        col_r4.info("⏭️ **Skip**\nBoth have gaps.")
        
        st.markdown("---")
        
        st.markdown("### Step-by-Step Example")
        st.markdown("**Motif:** `ATCGAC`")
        
        col_c1, col_c2 = st.columns(2)
        with col_c1:
            st.code("Reference: A - T C G A - C\nCompared:  A G T G - A - C", language="text")
        
        st.markdown("**Position-by-Position Analysis:**")
        df_hw = pd.DataFrame([
            ["A", "A", "Match", "1"],
            ["-", "G", "Gap", "-"],
            ["T", "T", "Match", "1"],
            ["C", "G", "Mismatch", "0"],
            ["G", "-", "Mismatch", "0"],
            ["A", "A", "Match", "1"],
            ["-", "-", "Skip", ""],
            ["C", "C", "Match", "1"],
        ], columns=["Reference Base", "Compared Base", "Status", "Vector Result"])
        st.table(df_hw)
        
        st.success("**Final Similarity Vector:** `1-10011` (length 7)\n\n**Conservation Percentage:** (4 matches / 7 positions) * 100% ≈ **57.14%**")
        
        st.markdown("### Genome Coordinates")
        st.markdown("In addition to the similarity vector, MAFin supplies the genomic coordinates of the motif. Consider the example above (using 0-indexed coordinates): after comparing the sequences (where gaps in both sequences are skipped), the ungapped reference motif is: `A, T, C, G, A, C` (a total of 6 bases). If this motif begins at position 1000 in the reference genome, its coordinates are:")
        st.code("Start: 1000\nEnd:   1005", language="text")
        st.info("🚨 **IMPORTANT:** The `End` coordinate in BED output files is **non-inclusive**. This means the exact same hit above will be written to the BED file as `Start: 1000, End: 1006`. Why non-inclusive? It's a BED format convention that makes length calculation trivial: **length = End − Start** (no need for +1). This avoids the off-by-one errors that plague inclusive coordinate systems.")

    with hw_tab3:
        st.subheader("Handling the Reverse Strand")
        st.markdown("""
        DNA is **double-stranded** — two complementary strands run in opposite directions. The bases always pair:
        - **A** ↔ **T** (Adenine pairs with Thymine)
        - **C** ↔ **G** (Cytosine pairs with Guanine)
        
        The **reverse complement** of a sequence is obtained by reversing the strand and swapping each base with its pair.
        For example, the reverse complement of `ATGC` is `GCAT` (reverse → `CGTA`, then swap: C→G, G→C, T→A, A→T → `GCAT`).
        
        Often, a biological motif might exist on the *opposite* strand. When you select **Include Reverse Complement: Yes**, MAFin handles this automatically depending on your input type:
        """)
        
        col_rc1, col_rc2, col_rc3 = st.columns(3)
        
        with col_rc1:
            with st.container(border=True):
                st.markdown("#### 🧬 Exact K-mers")
                st.markdown("We simply compute the reverse complement of your input string. This is straightforward because reversing a fixed string is unambiguous — there's only one possible reverse complement.")
                st.code("Input:  CCG\nSearch: CGG", language="text")
                
        with col_rc2:
            with st.container(border=True):
                st.markdown("#### 🧩 Regex Patterns")
                st.markdown("We reverse complement the *Reference sequence itself*, then search it with the original regex. Why? Because reversing regex *rules* would scramble their meaning — `A[GT]C` reversed as a string would become `C]TG[A`, which is nonsense. Reversing the target instead keeps the pattern logic intact.")
                st.code("Ref:    ATCGGCA\nRevRef: TGCCGAT\nRegex:  C{2}G", language="text")
                
        with col_rc3:
            with st.container(border=True):
                st.markdown("#### 📊 PWMs (JASPAR)")
                st.markdown("We mathematically reverse the probability matrix by swapping columns and rows (A↔T, C↔G). This works because a PWM is a pure mathematical object — reversing its rows and swapping complement pairs produces the exact reverse-strand equivalent.")
                st.code("Original Matrix\n      ↓\nReversed Matrix", language="text")


elif nav_selection == "🚀 Run Motif Search":
    st.title("🚀 Run Motif Search")
    st.markdown("Follow the steps below to configure and run your sequence search. The steps mirror the scientific workflow: first define your data, then your question, then your constraints, then execute.")

    if "analysis_complete" not in st.session_state:
        st.session_state.analysis_complete = False
    if "workspace_dir" not in st.session_state:
        st.session_state.workspace_dir = None

    st.subheader("📁 Step 1: Upload Alignment Data")
    col_u1, col_u2 = st.columns([2, 1])
    with col_u1:
        maf_file = st.file_uploader(
            "Upload MAF File (Required)",
            type=["maf", "txt"],
            key="maf_file_upload",
            help="Multiple Alignment Format (.maf) file containing your aligned genomes."
        )
    with col_u2:
        st.markdown("<br><br>", unsafe_allow_html=True)
        demo_maf = CURRENT_DIR / "test" / "test_alignment.maf"
        if demo_maf.exists():
            with open(demo_maf, "rb") as f:
                st.download_button(
                    label="📥 Download Sample Data",
                    data=f,
                    file_name="sample_alignment.maf",
                    mime="text/plain",
                    help="Don't have a MAF file? Download this sample to test MAFin.",
                    use_container_width=True
                )

    is_valid_maf = False
    if maf_file is not None:
        header = maf_file.getvalue()[:10].decode("utf-8", errors="ignore")
        if not header.startswith("##maf"):
            st.error("❌ Invalid format. File must start with '##maf'.")
        else:
            is_valid_maf = True

    st.markdown("---")
    st.subheader("🎯 Step 2: Define Search Pattern")
    
    search_type = st.radio(
        "Select Search Type",
        ["Exact K-mers", "Regular Expression (Regex)", "Position Weight Matrix (PWM)"],
        horizontal=True,
        key="search_type_radio",
        help="• **Exact K-mers**: Finds exact, 100% identical sequence matches only. Fastest method.\n\n• **Regular Expression (Regex)**: Allows for simple variation using pattern rules (e.g. `A[GT]C` matches `AGC` or `ATC`).\n\n• **Position Weight Matrix (PWM)**: Uses a probabilistic model (.jaspar file) to find complex biological motifs (like transcription factor binding sites) where some bases are strongly preferred but others are flexible."
    )

    search_input = None
    kmers_input = None
    regex_input = None
    jaspar_file = None
    search_type_param = None

    st.markdown("<br>", unsafe_allow_html=True)

    if search_type == "Exact K-mers":
        kmers_input = st.text_area(
            "Enter Exact K-mers (one per line)",
            placeholder="AAAATTTTGGGGCCCC\nTTTTTTTTTTTTTTTT",
            height=100,
            key="kmers_input_area",
            help="Type exact DNA sequences (e.g., 16-mers), one per line. Note: 'N' is allowed but matches literally (it is not a wildcard). For wildcard searches, use the Regex option below."
        )
        search_type_param = "kmers"
        search_input = None
        if kmers_input.strip():
            lines = [line.strip() for line in kmers_input.strip().split("\n") if line.strip()]
            if lines:
                k_length = len(lines[0])
                
                # Check for length consistency
                invalid_len_lines = [line for line in lines if len(line) != k_length]
                
                # Check for invalid DNA characters
                dna_chars = set("ACGTNacgtn")
                invalid_chars = set()
                for line in lines:
                    invalid_chars.update(set(line) - dna_chars)
                
                if invalid_len_lines:
                    st.error(f"❌ **Invalid Input:** All K-mers must be exactly the same length. You started with a {k_length}-mer, but also entered a {len(invalid_len_lines[0])}-mer (`{invalid_len_lines[0]}`).")
                elif invalid_chars:
                    st.error(f"❌ **Invalid DNA Sequence:** Sequences can only contain A, C, G, T, or N. Found invalid characters: `{', '.join(sorted(list(invalid_chars)))}`.")
                else:
                    search_input = kmers_input.strip()
        elif kmers_input:
            st.warning("⚠️ Please enter at least one valid sequence. (Only whitespace detected)")

    elif search_type == "Regular Expression (Regex)":
        regex_input = st.text_area(
            "Enter Regex Patterns (one per line)",
            placeholder="CTG[CC]+CGCA\nAGT",
            height=100,
            key="regex_input_area",
            help="Type patterns using standard regular expressions. For example, 'A[GT]C' means A followed by either G or T, followed by C."
        )
        search_input = None
        search_type_param = "regexes"
        
        if regex_input.strip():
            # Find EVERY character that is not a DNA base or whitespace
            import re as py_re
            non_dna_found = py_re.findall(r'[^ACGTNacgtn\s]', regex_input)
            
            if non_dna_found:
                unique_chars = sorted(list(set(non_dna_found)))
                meta_chars = set("[]().*+?^$|{}\\")
                
                # Identify characters that are not DNA bases and not regex meta-characters
                found_literal = [c for c in unique_chars if c not in meta_chars]
                
                if found_literal:
                    st.error(f"❌ **Invalid Genomic Characters:** `{', '.join(found_literal)}` are not standard DNA bases (A, C, G, T, N). These characters will not match anything in genomic data.")
            
            # Syntax Check
            syntax_error = False
            for line in regex_input.strip().split("\n"):
                if not line.strip(): continue
                try:
                    py_re.compile(line.strip())
                except Exception as e:
                    st.error(f"❌ **Regex Syntax Error:** `{line.strip()}` is invalid. Details: {e}")
                    syntax_error = True
                    break
            
            if not syntax_error:
                search_input = regex_input.strip()
        
        if regex_input and not regex_input.strip():
            st.warning("⚠️ Please enter at least one valid regex pattern. (Only whitespace detected)")

    else:
        jaspar_file = st.file_uploader(
            "Upload JASPAR Motif File (.jaspar)",
            type=["jaspar"],
            key="jaspar_file_upload",
            help="A file containing a probability matrix of base frequencies, usually obtained from the JASPAR database for studying transcription factors. Must start with a `>` header."
        )
        search_type_param = "jaspar_file"
        if jaspar_file:
            try:
                content = jaspar_file.getvalue().decode("utf-8")
                if not content.strip().startswith(">"):
                    st.error("❌ **Invalid JASPAR Content:** A valid JASPAR file must start with a `>` header line (e.g., `>MA0001.1 AGL3`).")
                    search_input = None
                else:
                    preview = "\n".join(content.splitlines()[:3])
                    st.code(f"{preview}\n...", language="text")
                    search_input = jaspar_file
            except Exception as e:
                st.error(f"❌ **Error Reading JASPAR File:** {e}")
                search_input = None

    st.markdown("---")
    st.subheader("⚙️ Step 3: Search Scope & Settings")
    
    col_basic1, col_basic2 = st.columns(2)
    with col_basic1:
        search_in = st.selectbox(
            "Search Target", 
            ["reference", "all"], 
            format_func=lambda x: "Reference Sequence" if x == "reference" else "All Aligned Sequences",
            help="• **Reference Sequence**: Only searches the primary (top) species in the alignment block. It then checks if the other aligned species share that matched sequence. This is faster and usually the preferred method.\n\n• **All Aligned Sequences**: Searches every single species in the block. Use this if you want to find the motif even if it's missing from the primary reference sequence."
        )
    with col_basic2:
        reverse_complement = st.radio(
            "Strands to Search", 
            ["no", "yes"], 
            format_func=lambda x: "Forward Strand Only" if x == "no" else "Both Strands (+/-)",
            horizontal=True, 
            help="DNA is double-stranded.\n\n• **Forward Strand Only**: Searches only the sequence exactly as it appears in the file.\n\n• **Both Strands**: Automatically searches the opposite, complementary strand as well (A↔T, C↔G). Since biological motifs can occur on either strand, selecting 'Both Strands' ensures you don't miss potential matches, though it does increase search time."
        )

    pvalue_threshold = 1e-4
    if search_type == "Position Weight Matrix (PWM)":
        pvalue_threshold = st.number_input(
            "P-value Threshold",
            min_value=1e-10,
            max_value=1e-1,
            value=1e-4,
            format="%.2e",
            help="A lower p-value demands stronger statistical evidence that the match isn't due to random chance. For example, 1e-4 means there is only a 1 in 10,000 probability the match is coincidental. Lower = stricter = fewer but more confident results.",
        )

    detailed_report = st.checkbox(
        "Generate Detailed Output (CSV/JSON) - recommended for downstream analysis", 
        value=True,
        help="The BED file alone only gives you coordinates. The CSV adds per-species conservation detail and similarity vectors that reveal exactly which bases match or differ — essential for any downstream biological interpretation. Check this to receive the full report."
    )

    with st.expander("⚙️ Advanced System Settings (Optional)", expanded=False):
        st.info("Most users will not need to change these system-level optimizations.")
        col_adv1, col_adv2 = st.columns(2)
        with col_adv1:
            cpu_count = multiprocessing.cpu_count()
            num_processes = st.slider(
                f"Analysis Speed (up to {cpu_count} available cores)", 
                1, 
                cpu_count, 
                max(1, cpu_count // 2), 
                help=f"MAFin detected {cpu_count} CPU cores on your machine. Allocating more cores will speed up processing for large MAF files, but may slow down other running apps."
            )
            bg_freqs_input = "0.25 0.25 0.25 0.25"
            if search_type == "Position Weight Matrix (PWM)":
                bg_freqs_input = st.text_input(
                    "Nucleotide Background Frequencies (A C G T)",
                    value="0.25 0.25 0.25 0.25",
                    help="The natural occurrence rate of bases in the target genome. Defaults assume equal distribution.",
                )
        with col_adv2:
            genome_ids_input = st.text_area(
                "Filter by Genome IDs (one per line)",
                placeholder="hg38\nmm10",
                height=100,
                key="genome_ids_input",
                help="Type specific genome IDs to limit the search, one per line. If left blank, all genomes are searched."
            )

    st.markdown("---")
    st.subheader("🚀 Step 4: Run Analysis")

    ready_to_run = is_valid_maf and (search_input is not None)

    if not ready_to_run:
        st.info("💡 Please upload a MAF file and define your search pattern above to continue.")
    else:
        assert maf_file is not None
        assert search_input is not None
        with st.container(border=True):
            st.markdown(f"**Ready to search!** You are looking for a **{search_type}** in **{'the Reference Sequence' if search_in == 'reference' else 'All Aligned Sequences'}**.")
            if reverse_complement == "yes":
                st.markdown("The search will include the **Reverse Complement**.")
            
    st.markdown("---")

    col_run, col_clear = st.columns([3, 1])

    with col_run:
        run_clicked = st.button(
            "🚀 Run MAFin Analysis",
            use_container_width=True,
            type="primary",
            disabled=not ready_to_run,
        )

    with col_clear:
        if st.button("🗑️ Clear Workspace", use_container_width=True):
            if st.session_state.workspace_dir and os.path.exists(st.session_state.workspace_dir):
                shutil.rmtree(st.session_state.workspace_dir, ignore_errors=True)
            st.session_state.workspace_dir = None
            st.session_state.analysis_complete = False
            st.toast("Workspace cleared successfully!", icon="🧹")
            st.rerun()

    if run_clicked:
        assert maf_file is not None
        if st.session_state.workspace_dir and os.path.exists(st.session_state.workspace_dir):
            shutil.rmtree(st.session_state.workspace_dir, ignore_errors=True)
        st.session_state.workspace_dir = tempfile.mkdtemp(prefix="mafin_analysis_")
        st.session_state.analysis_complete = False

        with st.status("Running MAFin analysis...", expanded=True) as status:
            try:
                status.write("Initializing secure temporary workspace...")
                if search_type_param == "kmers" and kmers_input is not None:
                    search_data = kmers_input.strip().encode("utf-8")
                elif search_type_param == "regexes" and regex_input is not None:
                    search_data = regex_input
                elif search_type_param == "jaspar_file" and jaspar_file is not None:
                    search_data = (jaspar_file.getbuffer(), jaspar_file.name)
                else:
                    raise ValueError("Unsupported search input")

                if genome_ids_input and genome_ids_input.strip():
                    genome_ids_data = genome_ids_input.strip().encode("utf-8")
                else:
                    genome_ids_data = None

                result, cmd = execute_mafin_analysis(
                    tmpdir=st.session_state.workspace_dir,
                    maf_file_obj=maf_file,
                    maf_name=maf_file.name,
                    search_type_param=search_type_param,
                    search_data=search_data,
                    search_in=search_in,
                    rev_comp=reverse_complement,
                    procs=num_processes,
                    pval=pvalue_threshold,
                    bg_freqs=bg_freqs_input,
                    detailed=detailed_report,
                    genome_ids_data=genome_ids_data,
                )

                status.write(f"Executed: {' '.join(cmd)}")

                if result.returncode == 0:
                    status.update(label="Analysis completed successfully!", state="complete", expanded=False)
                    st.session_state.analysis_complete = True
                    st.session_state.analysis_time = result.elapsed_time
                    st.session_state.analysis_memory = result.max_memory_mb
                else:
                    status.update(label="Analysis failed.", state="error", expanded=True)
                    st.error(result.stderr)
                    st.session_state.analysis_complete = False

            except Exception as exc:
                status.update(label="Critical error encountered.", state="error", expanded=True)
                st.error(f"Details: {exc}")
                st.session_state.analysis_complete = False

    if st.session_state.analysis_complete and st.session_state.workspace_dir and os.path.exists(st.session_state.workspace_dir):
        render_analysis_results(
            st.session_state.workspace_dir, 
            st.session_state.get('analysis_time'), 
            st.session_state.get('analysis_memory')
        )

