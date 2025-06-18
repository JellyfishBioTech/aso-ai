import streamlit as st
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np


st.set_page_config(layout="wide")

@st.cache_data
def load_data():
    df = pd.read_csv("/Users/julia_patsiukova/Downloads/aso_dashboard/ASO_table_no_mod.csv")  # Replace with your actual filename
    df = df.dropna(subset=["efficiency", "concentration"])
    df["concentration"] = pd.to_numeric(df["concentration"], errors="coerce")
    df["efficiency"] = pd.to_numeric(df["efficiency"], errors="coerce")
    return df

df = load_data()

# Standardize source_type for filter
df["source_type"] = df["source_type"].str.strip().str.lower()

# Clean species values into buckets
def clean_species(val):
    val = str(val).lower()
    if "human" in val and "mouse" in val:
        return "human and mouse"
    elif "human" in val:
        return "human"
    elif "mouse" in val:
        return "mouse"
    elif "monkey" in val:
        return "monkey"
    elif "dog" in val:
        return "dog"
    elif "unspecified" in val:
        return "unspecified"
    else:
        return "unspecified"

df["species_clean"] = df["species"].apply(clean_species)

# Standardize modification values
def clean_modification(val):
    val = str(val).lower()
    if "2ome" in val or "2'-ome" in val:
        return "2OMe"
    elif "pmo" in val and "modif" in val:
        return "modified PMO"
    elif "ppmo" in val:
        return "modified PMO (PPMO)"
    elif "unmod" in val:
        return "unmodified PMO"
    elif "tmo" in val or "tc-dna" in val:
        return "Other"
    elif val == "nan":
        return "unspecified"
    else:
        return val.strip()

df["modification_clean"] = df["modification"].apply(clean_modification)

st.title("ASO Dashboard Explorer")

st.sidebar.header("🔍 Filter Options")

# source_type = st.sidebar.selectbox("Source Type", options=df["source_type"].dropna().unique())
# species = st.sidebar.multiselect("Species", sorted(df["species"].dropna().unique()))
# cell_line = st.sidebar.multiselect("Cell Line", sorted(df["cell_line"].dropna().unique()))
# delivery_approach = st.sidebar.multiselect("Delivery Approach", sorted(df["delivery_approach"].dropna().unique()))
# modification = st.sidebar.multiselect("Modification", sorted(df["modification"].dropna().unique()))

# Source Type
source_type_options = ["All"] + sorted(df["source_type"].dropna().unique())
selected_source_type = st.sidebar.selectbox("Source Type", options=source_type_options)

# Target gene
target_gene = st.sidebar.multiselect("Target Gene", sorted(df["target_gene"].dropna().unique()))

# Species (cleaned + single-select)
species_options = sorted(df["species_clean"].dropna().unique())
selected_species = st.sidebar.selectbox("Species", options=["All"] + species_options)

# DNA/RNA (replacing Cell Line)
selected_molecule = st.sidebar.selectbox("DNA or RNA", options=["All"] + sorted(df["DNA.RNA"].dropna().unique()))

# Exon or Intron (replacing Delivery Approach)
selected_exon_intron = st.sidebar.multiselect("Exon or Intron", sorted(df["exon_or_intron"].dropna().unique()))

# Modification (cleaned)
modification = st.sidebar.multiselect("Modification", sorted(df["modification_clean"].dropna().unique()))

conc_min, conc_max = st.sidebar.slider(
    "Concentration Range", float(df["concentration"].min()), float(df["concentration"].max()),
    (float(df["concentration"].min()), float(df["concentration"].max()))
)
eff_min, eff_max = st.sidebar.slider(
    "Efficiency Range (%)", float(df["efficiency"].min()), float(df["efficiency"].max()),
    (float(df["efficiency"].min()), float(df["efficiency"].max()))
)

filtered_df = df.copy()

if selected_source_type != "All":
    filtered_df = filtered_df[filtered_df["source_type"] == selected_source_type]

if selected_species != "All":
    filtered_df = filtered_df[filtered_df["species_clean"] == selected_species]

if selected_molecule != "All":
    filtered_df = filtered_df[filtered_df["DNA.RNA"] == selected_molecule]

if selected_exon_intron:
    filtered_df = filtered_df[filtered_df["exon_or_intron"].isin(selected_exon_intron)]

if modification:
    filtered_df = filtered_df[filtered_df["modification_clean"].isin(modification)]

filtered_df = filtered_df[
    (filtered_df["concentration"] >= conc_min) & (filtered_df["concentration"] <= conc_max) &
    (filtered_df["efficiency"] >= eff_min) & (filtered_df["efficiency"] <= eff_max)
]

# st.subheader("1. Average Efficiency by Source Type")
# # fig1 = px.bar(
# #     df.groupby("source_type")["efficiency"].mean().reset_index(),
# #     x="source_type", y="efficiency", color="source_type"
# # )
# df_bar = df.groupby("source_type", as_index=False)["efficiency"].mean()
# fig1 = px.bar(df_bar, x="source_type", y="efficiency", color="source_type", text="efficiency")
# fig1.update_layout(yaxis_title="Mean Efficiency (%)")
# fig1.update_traces(texttemplate='%{text:.2f}', textposition='outside')
# st.plotly_chart(fig1, use_container_width=True)

st.subheader("1. Efficiency Distribution by Source Type")

# fig1 = px.histogram(
#     df, 
#     x="efficiency", 
#     color="source_type", 
#     marginal="box",         # Adds box plot on top
#     nbins=50, 
#     opacity=0.7,
#     histnorm='percent'      # Optional: normalize to %
# )
# fig1.update_layout(
#     xaxis_title="Efficiency (%)",
#     yaxis_title="Distribution (%)",
#     barmode="overlay"
# )

# Prepare the subplot
fig1 = make_subplots(rows=1, cols=2, shared_yaxes=True, subplot_titles=["Journal Articles", "Patents"])

# Journal articles
journal = df[df["source_type"] == "journal article"]
fig1.add_trace(
    go.Histogram(x=journal["efficiency"], nbinsx=50, name="Journal", marker_color="#1f77b4", histnorm='percent'),
    row=1, col=1
)

# Patents
patent = df[df["source_type"] == "patent"]
fig1.add_trace(
    go.Histogram(x=patent["efficiency"], nbinsx=50, name="Patent", marker_color="#aec7e8", histnorm='percent'),
    row=1, col=2
)

fig1.update_layout(
    title_text="Efficiency Distribution by Source Type",
    xaxis_title="Efficiency (%)",
    xaxis2_title="Efficiency (%)",
    yaxis_title="Distribution (%)",
    showlegend=False
)
st.plotly_chart(fig1, use_container_width=True)

st.subheader("2. Efficiency vs. Concentration")
# fig2 = px.scatter(
#     filtered_df, x="concentration", y="efficiency", color="modification",
#     hover_data=["delivery_approach", "cell_line", "link"]
# )
# fig2 = px.scatter(
#     filtered_df, x="concentration", y="efficiency", color="modification",
#     hover_data=["target_gene", "source_type", "delivery_approach", "cell_line", "link"]
# )
# fig2.update_layout(xaxis_type="log", xaxis_title="Concentration (log scale)", yaxis_title="Efficiency (%)")
# Bin concentration (log scale)
# filtered_df["conc_bin"] = pd.qcut(np.log10(filtered_df["concentration"]), q=5, duplicates="drop")
# Remove rows with NaN or non-positive concentration values
concentration_data = filtered_df[filtered_df["concentration"] > 0].copy()

# Apply log10 and create bins
concentration_data["log_conc"] = np.log10(concentration_data["concentration"])
concentration_data["conc_bin"] = pd.qcut(
    concentration_data["log_conc"], 
    q=5, 
    duplicates="drop"
).astype(str)

# Aggregate mean efficiency per concentration level
grouped = concentration_data.groupby("concentration")["efficiency"].mean().reset_index()

fig2 = px.scatter(
    grouped,
    x="concentration",
    y="efficiency",
    title="Mean Efficiency by Concentration",
    log_x=True
)
fig2.update_layout(
    xaxis_title="Concentration (log scale)",
    yaxis_title="Mean Efficiency (%)"
)

st.plotly_chart(fig2, use_container_width=True)

# st.subheader("3. ASO Counts per Target Gene (by Source Type)")
# fig3 = px.histogram(filtered_df, x="target_gene", color="source_type", barmode="group")
# st.plotly_chart(fig3, use_container_width=True)

st.subheader("3. ASO Counts per Target Gene (by Source Type)")

# Limit to top 10 most frequent target genes
top_genes = filtered_df["target_gene"].value_counts().nlargest(10).index
df_top_genes = filtered_df[filtered_df["target_gene"].isin(top_genes)]

fig3 = px.histogram(
    df_top_genes,
    x="target_gene",
    color="source_type",
    barmode="group"
)
fig3.update_layout(
    yaxis_title="Count of ASOs",
    xaxis_title="Target Gene"
)
st.plotly_chart(fig3, use_container_width=True)

# st.subheader("4. Species Distribution")
# fig4 = px.pie(filtered_df, names="species", title="Species Breakdown")
# st.plotly_chart(fig4, use_container_width=True)

st.subheader("4. Species Distribution")

species_counts = filtered_df["species"].value_counts().reset_index()
species_counts.columns = ["species", "count"]

fig4 = px.pie(
    species_counts,
    values="count",
    names="species",
    title="Species Breakdown",
    hole=0.3
)
fig4.update_traces(textinfo='percent+label')
st.plotly_chart(fig4, use_container_width=True)

# st.subheader("Filtered Data Table")
# st.dataframe(filtered_df[["target_gene", "species", "modification", "concentration", "efficiency", "link"]])
st.subheader("Full ASO Dataset (Unfiltered)")
st.dataframe(df[["target_gene", "number_exon_intron", "exon_or_intron", "species", "cell_line",	
                 "delivery_approach", "aso_type",	"oligo_sequence", "modification", "concentration", 
                 "efficiency", "DNA.RNA", "link", "top1_hit",	"top1_real", "top1_predicted", "no_mod_ASO"]])

# with st.expander("Show Advanced Columns"):
#     st.dataframe(filtered_df[["oligo_sequence", "top1_predicted", "no_mod_ASO"]])
