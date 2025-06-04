import streamlit as st
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import re

# ------------------ Helper Functions ------------------

def detect_g_quadruplex(seq):
    pattern = r'(G{3,}\w{1,7}){3,}G{3,}'
    return [(m.start(), m.end()) for m in re.finditer(pattern, str(seq))]

def detect_z_dna(seq):
    pattern = r'(CG){6,}'
    return [(m.start(), m.end()) for m in re.finditer(pattern, str(seq))]

def detect_cruciform(seq):
    results = []
    for i in range(len(seq) - 8):
        left = seq[i:i+4]
        right = seq[i+4:i+8][::-1].translate(str.maketrans('ATCG', 'TAGC'))
        if left == right:
            results.append((i, i+8))
    return results

def detect_triplex(seq):
    pattern = r'[AG]{10,}'
    return [(m.start(), m.end()) for m in re.finditer(pattern, str(seq))]

def calculate_similarity_score(detected_motifs, known_motifs):
    match_count = 0
    for motif_type in detected_motifs:
        for m1 in detected_motifs[motif_type]:
            if any(abs(m1[0] - m2[0]) < 10 for m2 in known_motifs.get(motif_type, [])):
                match_count += 1
    return round((match_count / (sum(len(v) for v in known_motifs.values()) + 1)) * 100, 2)

# ------------------ Streamlit UI ------------------

st.title("🔬 DNA Structural Motif-Based Biomarker Prediction")
st.markdown("""
This tool analyzes **cis-regulatory DNA sequences** for **non-B DNA motifs** associated with diseases such as cancer, Alzheimer's, etc.
""")

uploaded_file = st.file_uploader("📂 Upload your promoter sequence file (any format):")

known_disease_motifs = {
    "G4": [(100, 120), (300, 320)],
    "Z-DNA": [(150, 165)],
    "Cruciform": [(200, 215)],
    "Triplex": [(250, 270)]
}

if uploaded_file:
    try:
        # Attempt to parse as FASTA, fallback if fails
        try:
            fasta_sequences = SeqIO.parse(uploaded_file, "fasta")
            records = list(fasta_sequences)
            if not records:
                raise ValueError("No sequences found in FASTA format.")
        except Exception:
            # Not fasta, try to read as plain text sequence
            uploaded_file.seek(0)
            seq = uploaded_file.read().decode("utf-8").strip().replace("\n", "").replace(" ", "").upper()
            records = [seq]

        for idx, record in enumerate(records):
            if hasattr(record, "id"):  # FASTA record
                st.subheader(f"🧬 Analyzing: {record.id}")
                sequence = str(record.seq)
            else:  # plain sequence string
                st.subheader(f"🧬 Analyzing: Sequence {idx+1}")
                sequence = record

            motifs = {
                "G4": detect_g_quadruplex(sequence),
                "Z-DNA": detect_z_dna(sequence),
                "Cruciform": detect_cruciform(sequence),
                "Triplex": detect_triplex(sequence)
            }

            st.markdown("### 🧪 Detected Motifs")
            for motif, positions in motifs.items():
                st.write(f"**{motif}:** {len(positions)} regions found")

            motif_density = [0] * len(sequence)
            for positions in motifs.values():
                for start, end in positions:
                    for i in range(start, min(end, len(sequence))):
                        motif_density[i] += 1

            fig, ax = plt.subplots()
            ax.plot(motif_density)
            ax.set_xlabel("Position")
            ax.set_ylabel("Motif Count")
            ax.set_title("Motif Density Across Sequence")
            st.pyplot(fig)

            similarity = calculate_similarity_score(motifs, known_disease_motifs)
            st.markdown(f"### 🧠 Disease Likelihood Score: **{similarity}%** match with known disease motifs")

            if similarity > 70:
                st.success("High similarity to known disease-associated promoters.")
            elif similarity > 40:
                st.warning("Moderate similarity. Experimental validation recommended.")
            else:
                st.info("Low similarity to known disease motifs")

    except Exception as e:
        st.error(f"❌ Error while processing file: {e}")
