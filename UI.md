# UI

`dashboard_starter.py` is an early-alpha Streamlit front end for the
convergence model.

## Prerequisites

Streamlit is included in the project requirements:

```bash
pip install -r requirements.txt
```

Or install just what the dashboard needs:

```bash
pip install streamlit matplotlib numpy scipy
```

## Running

From the repository root (the import of `simulation.convergence_model`
resolves relative to it):

```bash
streamlit run dashboard_starter.py
```

Pick a prediction duration in the sidebar and press **Run Simulation**.
