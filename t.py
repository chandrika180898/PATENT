import streamlit as st
import pandas as pd
import re
import math
from sklearn.ensemble import RandomForestClassifier
import altair as alt
import pickle
import smtplib
from email.message import EmailMessage

# ------------------------- PAGE CONFIG --------------------------
st.set_page_config(
    page_title="Welcome to the DNA Motif & Perplexity Analyzer",
    layout="wide",
    page_icon="🧬"
)

# ------------------------- CSS Styling --------------------------
st.markdown("""
    <style>
    .main {
        background-color: #F0F2F6;
    }
    .css-1rs6os.edgvbvh3 {
        font-size: 18px;
    }
    .block-container {
        padding-top: 2rem;
    }
    .stDownloadButton {
        background-color: #4CAF50;
        color: white;
        padding: 0.5em 1em;
        border-radius: 8px;
        text-align: center;
    }
    </style>
""", unsafe_allow_html=True)

# ------------------------- REGISTRATION & LOGIN --------------------------
if "registered_users" not in st.session_state:
    st.session_state["registered_users"] = {}

if "logged_in" not in st.session_state:
    st.session_state["logged_in"] = False

with st.sidebar:
    if not st.session_state["logged_in"]:
        mode = st.radio("Select Mode", ["Register", "Login"])

        if mode == "Register":
            st.markdown("### 📝 Register")
            reg_username = st.text_input("Choose a Username")
            reg_password = st.text_input("Choose a Password", type="password")
            reg_email = st.text_input("Email")
            if st.button("Register"):
                if not reg_username or not reg_password or not reg_email:
                    st.warning("Please fill all fields")
                elif reg_username in st.session_state["registered_users"]:
                    st.warning("🚫 Username already exists!")
                else:
                    st.session_state["registered_users"][reg_username] = {
                        "password": reg_password,
                        "email": reg_email
                    }

                    # Mock sending confirmation email
                    try:
                        msg = EmailMessage()
                        msg.set_content(f"Welcome {reg_username}, you have successfully registered for the DNA Motif Analyzer.")
                        msg["Subject"] = "Registration Successful"
                        msg["From"] = "dnatool@genomics.com"
                        msg["To"] = reg_email

                        # Simulated email (no SMTP server used)
                        print(f"Sending email to {reg_email}: {msg.get_content()}")
                        st.success("✅ Registered successfully! Please switch to Login tab.")
                    except Exception as e:
                        st.warning(f"Registered, but email failed: {e}")

            st.stop()

        elif mode == "Login":
            st.markdown("### 🔐 Login")
            username = st.text_input("Username")
            password = st.text_input("Password", type="password")
            login_button = st.button("Login")
            if login_button:
                user_data = st.session_state["registered_users"].get(username)
                if user_data and user_data["password"] == password:
                    st.session_state["logged_in"] = True
                    st.session_state["current_user"] = username
                    st.success(f"Welcome back, {username}!")
                else:
                    st.warning("Invalid username or password")
                    st.stop()
            if not st.session_state.get("logged_in"):
                st.warning("Please login to continue")
                st.stop()
    else:
        st.markdown(f"### 👤 Logged in as {st.session_state['current_user']}")
        if st.button("Logout"):
            st.session_state["logged_in"] = False
            st.session_state["current_user"] = None
            st.experimental_rerun()

# ------------------------- SIDEBAR NAVIGATION --------------------------
if st.session_state["logged_in"]:
    st.sidebar.title("🧭 Navigation")
    page = st.sidebar.radio("Choose Page", ["🏠 Home", "📂 Upload & Analyze", "📊 Results", "📥 Download Report", "📞 Contact"])
else:
    # User not logged in, no further UI
    st.stop()

# ------------------------- HELPER FUNCTIONS --------------------------

def calculate_perplexity(sequence, k=3):
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    kmer_counts = pd.Series(kmers).value_counts()
    total = sum(kmer_counts)
    probs = kmer_counts / total
    entropy = -sum(p * math.log2(p) for p in probs)
    return 2 ** entropy

def detect_g_quadruplex(seq):
    return len(re.findall(r'(G{3,}\w{1,7}){3,}G{3,}', seq))

def detect_z_dna(seq):
    return len(re.findall(r'(CG){6,}', seq))

def reverse_complement(seq):
    complement = str.maketrans('ATGCatgc', 'TACGtacg')
    return seq.translate(complement)[::-1]

def detect_cruciform(seq, min_len=4, max_len=6, spacer=10):
    count = 0
    seq = seq.upper()
    for size in range(min_len, max_len + 1):
        for i in range(len(seq) - 2 * size - spacer + 1):
            left = seq[i:i+size]
            right = seq[i+size+spacer:i+2*size+spacer]
            if reverse_complement(left) == right:
                count += 1
    return count

def detect_tata_box(seq):
    return len(re.findall(r'TATA[AT]A[AT]', seq))

def detect_direct_repeats(seq):
    return len(re.findall(r'(.{3,6})\1+', seq))

def extract_features(seq):
    seq = seq.upper()
    return {
        'Perplexity': calculate_perplexity(seq),
        'G-Quadruplex': detect_g_quadruplex(seq),
        'Z-DNA': detect_z_dna(seq),
        'Cruciform': detect_cruciform(seq),
        'TATA-Box': detect_tata_box(seq),
        'Direct Repeats': detect_direct_repeats(seq),
        'Sequence Length': len(seq)
    }

# ------------------------- PAGES --------------------------

if page == "🏠 Home":
    st.markdown("## 🧬 Welcome to the **DNA Motif & Perplexity Analyzer**")
    st.markdown("""
    This web-based tool is built for analyzing DNA sequences and detecting **non-B DNA structures**, 
    calculating **perplexity**, and identifying key **regulatory motifs** such as:
    
    - 🔁 **Direct Repeats**
    - 🧷 **G-Quadruplexes**
    - 🔀 **Z-DNA**
    - 🎯 **TATA Boxes**
    - 🧬 **Cruciform Structures**

    Go to **📂 Upload & Analyze** to begin!
    """)

elif page == "📂 Upload & Analyze":
    st.markdown("## 📂 Upload & Analyze DNA Sequences")

    uploaded_file = st.file_uploader("📄 Upload a `.txt` file with sequences (one per line)", type=["txt"])

    if uploaded_file:
        st.success(f"📁 File uploaded: `{uploaded_file.name}`")

        content = uploaded_file.read().decode("utf-8").strip().splitlines()
        records = []
        for i, line in enumerate(content):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            seq_id = parts[0]
            sequence = "".join(parts[1:]) if len(parts) > 1 else parts[0]
            records.append((seq_id, sequence))

        with st.spinner("🔍 Extracting features..."):
            feature_rows = []
            for seq_id, sequence in records:
                try:
                    features = extract_features(sequence)
                    features["ID"] = seq_id
                    feature_rows.append(features)
                except Exception as e:
                    st.error(f"Error processing {seq_id}: {e}")

            df = pd.DataFrame(feature_rows)
            df = df[["ID"] + [col for col in df.columns if col != "ID"]]  # Reorder

            # Dummy ML prediction (always zeros)
            clf = RandomForestClassifier()
            clf.fit(df.drop(columns=["ID"]), [0]*len(df))
            df["Prediction"] = clf.predict(df.drop(columns=["ID"]))

            st.session_state["results_df"] = df
            st.success("✅ Analysis complete! Check the 📊 Results tab.")

elif page == "📊 Results":
    st.markdown("## 📊 Results Summary")
    if "results_df" in st.session_state:
        st.dataframe(st.session_state["results_df"], use_container_width=True)
    else:
        st.warning("📂 Please upload and analyze sequences first.")

elif page == "📥 Download Report":
    st.markdown("## 📥 Download Your Results")
    if "results_df" in st.session_state:
        csv = st.session_state["results_df"].to_csv(index=False).encode("utf-8")
        st.download_button("📥 Download as CSV", csv, file_name="motif_results.csv", mime="text/csv")
    else:
        st.warning("📂 No data available to download yet.")

elif page == "📞 Contact":
    st.markdown("## 📞 Contact Information")
    st.markdown("""
    - 👨‍🔬 **Dr. Y V Rajesh**  
      📧 yvrajesh_bt@kluniversity.in

    - 👩‍🔬 **G. Aruna Sesha Chandrika**  
      📧 chandrikagummadi1@gmail.com
    """)
