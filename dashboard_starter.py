# ui/dashboard_starter.py

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from simulation.convergence_model import ConvergencePredictor

st.set_page_config(page_title="CEED Dashboard", layout="centered")
st.title("🛰️ CEED: Cascading Energetic Event Disruption")
st.caption("Monitoring nonlinear energy convergence. From the comfort of your GitHub.")

# Initialize model
predictor = ConvergencePredictor()

# Sidebar options
years = st.sidebar.slider("Prediction Duration (years)", 1, 5, 3)
run_sim = st.sidebar.button("Run Simulation")

if run_sim:
    st.subheader("📈 Energy Evolution Over Time")
    t, solution = predictor.predict_convergence(years=years)
    total_energy, phases = predictor.classify_phases(solution)

    labels = ['Solar', 'Magnetic', 'Atmospheric', 'Oceanic']
    colors = ['gold', 'blue', 'green', 'teal']

    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(4):
        ax.plot(t, solution[:, i], label=labels[i], color=colors[i])
    ax.plot(t, total_energy, label="Total Energy", color='black', linewidth=2)
    ax.set_xlabel("Time (Years from now)")
    ax.set_ylabel("Energy Level")
    ax.set_title("System Energy Levels")
    ax.legend()
    st.pyplot(fig)

    st.subheader("🔍 Final Status")
    final_E = total_energy[-1]
    st.write(f"**Total System Energy (end of sim):** {final_E:.2f}")

    phase_map = {1: "Phase 1: Elevated", 2: "Phase 2: Coupling", 
                 3: "Phase 3: Nonlinear", 4: "Phase 4: Runaway"}
    final_phase = phase_map[phases[-1]]
    st.write(f"**Predicted System Phase:** `{final_phase}`")
