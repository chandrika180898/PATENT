# DNA Structural Motif-Based Biomarker Platform (Prototype in Python)

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import io
from Bio import SeqIO

# ------------------ Helper Functions ------------------

def detect_g_quadruplex(seq):
    """Detects potential G-quadruplex motifs."""
    pattern = r'(G{3,}\w{1,7}){3,}G{3,}'
    return [(m.start(), m.end()) for m in re.finditer(pattern, str(seq))]

def detect_z_dna(seq):
    """Detect Z-DNA potential zones (CG repeats)."""
    pattern = r'(CG){6,}'
    return [(m.start(), m.end()) for m in re.finditer(pattern, str(seq))]

def detect_cruciform(seq):
    """Detect potential inverted repeats forming cruciforms."""
    pattern = r'([ATCG]{4,})(?=.{0,10}\1[::-1])'
    return [(m.start(), m.end()) for m in re.finditer(pattern, str(seq))]

def detect_triplex(seq):
    """Detect potential homopurine stretches for triplex formation."""
    pattern = r'[AG]{10,}'
    return [(m.start(), m.end()) for m in re.finditer(pattern, str(seq))]

def calculate_similarity_score(detected_motifs, known_motifs):
    """Simple overlap-based similarity score between query and known disease motifs."""
    match_count = 0
    for motif_type in detected_motifs:
        for m1 in detected_motifs[motif_type]:
            if any(abs(m1[0] - m2[0]) < 10 for m2 in known_motifs.get(motif_type, [])):
                match_count += 1
    return round((match_count / (sum(len(v) for v in known_motifs.values()) + 1)) * 100, 2)

# ------------------ Streamlit UI ------------------

st.set_page_config(page_title="DNA Motif Biomarker Predictor", layout="wide")

st.title("🔬 DNA Structural Motif-Based Biomarker Prediction")
st.markdown("""
This tool analyzes **cis-regulatory DNA sequences** for **non-B DNA motifs** associated with diseases such as cancer, Alzheimer's, etc.
""")

# Upload FASTA file
uploaded_file = st.file_uploader("📂 Upload your promoter sequence in FASTA format:", type=["fasta", "fa"])

# Load known disease motifs (mocked for now)
known_disease_motifs = {
    "G4": [(100, 120), (300, 320)],
    "Z-DNA": [(150, 165)],
    "Cruciform": [(200, 215)],
    "Triplex": [(250, 270)]
}

if uploaded_file:
    # Decode binary file to text stream
    stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
    try:
        fasta_sequences = SeqIO.parse(stringio, "fasta")

        for record in fasta_sequences:
            st.subheader(f"🧬 Analyzing: {record.id}")
            sequence = str(record.seq)

            # Detect motifs
            motifs = {
                "G4": detect_g_quadruplex(sequence),
                "Z-DNA": detect_z_dna(sequence),
                "Cruciform": detect_cruciform(sequence),
                "Triplex": detect_triplex(sequence)
            }

            st.markdown("### 🧪 Detected Motifs")
            for motif, positions in motifs.items():
                st.write(f"**{motif}:** {len(positions)} regions found")

            # Plot motif density
            st.markdown("### 📊 Motif Density Plot")
            motif_density = [0] * len(sequence)
            for positions in motifs.values():
                for start, end in positions:
                    for i in range(start, min(end, len(sequence))):
                        motif_density[i] += 1

            fig, ax = plt.subplots()
            ax.plot(motif_density, color='darkgreen')
            ax.set_xlabel("Position")
            ax.set_ylabel("Motif Count")
            ax.set_title("Motif Density Across Sequence")
            st.pyplot(fig)

            # Compare with known disease motifs
            similarity = calculate_similarity_score(motifs, known_disease_motifs)
            st.markdown(f"### 🧠 Disease Likelihood Score: **{similarity}%** match with known disease motifs")

            if similarity > 70:
                st.success("✅ High similarity to known disease-associated promoters.")
            elif similarity > 40:
                st.warning("⚠️ Moderate similarity. Experimental validation recommended.")
            else:
                st.info("ℹ️ Low similarity to known disease motifs.")

    except Exception as e:
        st.error(f"❌ Error while processing FASTA file: {e}")
else:
    st.info("👆 Please upload a FASTA file to begin analysis.")
