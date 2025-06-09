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

with st.sidebar:
    # If logged in, show logout button and user info
    if st.session_state.get("logged_in"):
        st.markdown(f"### 👋 Welcome, {st.session_state.get('current_user')}")
        if st.button("Logout"):
            # Clear login session state
            for key in ["logged_in", "current_user"]:
                if key in st.session_state:
                    del st.session_state[key]
            # Rerun app safely after logout
            st.experimental_rerun()
    else:
        # Not logged in, show Register/Login tabs
        mode = st.radio("Select Mode", ["Register", "Login"])

        if mode == "Register":
            st.markdown("### 📝 Register")
            reg_username = st.text_input("Choose a Username")
            reg_password = st.text_input("Choose a Password", type="password")
            reg_email = st.text_input("Email")
            if st.button("Register"):
                if reg_username in st.session_state["registered_users"]:
                    st.warning("🚫 Username already exists!")
                else:
                    st.session_state["registered_users"][reg_username] = {
                        "password": reg_password,
                        "email": reg_email
                    }

                    # Send confirmation email (mocked for security)
                    try:
                        msg = EmailMessage()
                        msg.set_content(f"Welcome {reg_username}, you have successfully registered for the DNA Motif Analyzer.")
                        msg["Subject"] = "Registration Successful"
                        msg["From"] = "dnatool@genomics.com"
                        msg["To"] = reg_email

                        # Simulated email (no actual SMTP used for safety)
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
                    st.experimental_rerun()
                else:
                    st.warning("Invalid username or password")
                    st.stop()
            if not st.session_state.get("logged_in"):
                st.warning("Please login to continue")
                st.stop()

# Continue with the rest of your app after successful login...
