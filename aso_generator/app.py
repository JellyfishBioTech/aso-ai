# app.py
import os
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from aso_generator import ASOGenerator

st.set_page_config(page_title="ASO Generator", layout="wide")

st.title("ASO Generator")

with st.sidebar:
    st.header("Design Parameters")
    # transcript_id = st.text_input("Ensembl Transcript ID", value="ENST00000357033")
    input_mode = st.radio("Choose input type:", ["Ensembl Transcript ID", "Gene Name"])

    if input_mode == "Ensembl Transcript ID":
        transcript_id = st.text_input("Transcript ID (e.g., ENST00000357033)")
        gene_name = None
    else:
        gene_name = st.text_input("Gene Name (e.g., TP53)")
        transcript_id = None
        selected_transcript_display = None
        if gene_name:
            try:
                gen_temp = ASOGenerator(gene_name=gene_name)
                gen_temp.load_mane_transcript(gene_name)
                transcript_options = [t["display"] for t in gen_temp.transcript_choices]
                selected_transcript_display = st.selectbox("Choose Transcript", transcript_options)
                if selected_transcript_display:
                    transcript_id = next(t["id"] for t in gen_temp.transcript_choices if t["display"] == selected_transcript_display)
            except Exception as e:
                st.error(f"Could not fetch transcripts for gene {gene_name}: {e}")

    aso_type = st.selectbox("ASO Type", ["SSO", "GAPMER"])
    region = st.selectbox("Target Region", ["exon", "full"])
    target_exon = st.number_input("Target Exon", min_value=1, value=3) if region == "exon" else None

    st.markdown("---")
    window_min = st.slider("Minimum Window Size", 15, 25, 18)
    window_max = st.slider("Maximum Window Size", 18, 30, 25)
    gc_min = st.slider("Minimum GC%", 30, 70, 40)
    gc_max = st.slider("Maximum GC%", 40, 80, 65)
    homopolymer_len = st.number_input("Max Homopolymer Length", min_value=2, max_value=6, value=4)

    run_button = st.button("Generate ASOs")

if run_button:
    with st.spinner("Running ASO Generator... This may take a minute"):
        try:
            gen = ASOGenerator(
                transcript_id=transcript_id,
                gene_name=gene_name,
                aso_type=aso_type,
                region=region,
                target_exon=target_exon,
                window_min=window_min,
                window_max=window_max,
                GC_min=gc_min,
                GC_max=gc_max,
                homopolymer_len=homopolymer_len,
            )
            gen.run()

            df = pd.read_csv(os.path.join(gen.output_dir, "ASO_candidates.csv"))

            st.success(f"{len(df)} ASO candidates generated.")
            st.markdown(f"**Gene:** {gen.gene_name or 'N/A'}")
            st.dataframe(df)

            csv_data = df.to_csv(index=False).encode('utf-8')
            st.download_button("⬇Download CSV", csv_data, file_name="ASO_candidates.csv")

            # st.subheader("GC% Distribution")

            # fig, ax = plt.subplots()
            # df["GC%"].plot(kind="hist", bins=15, edgecolor="black", alpha=0.7, ax=ax)
            # ax.axvline(df["GC%"].mean(), color='red', linestyle='dashed', linewidth=1)
            # ax.set_title("Distribution of GC% in ASO Candidates")
            # ax.set_xlabel("GC%")
            # ax.set_ylabel("Frequency")
            # st.pyplot(fig)

            # st.subheader("📏 ASO Length Boxplot")

            # fig, ax = plt.subplots()
            # ax.boxplot(df["Length"], vert=False)
            # ax.set_title("Boxplot of ASO Length")
            # ax.set_xlabel("Length (nt)")
            # st.pyplot(fig)

            # st.subheader("Tm Distribution of ASOs")

            # fig, ax = plt.subplots()
            # ax.hist(df["Tm"], bins=20, color='skyblue', edgecolor='black')
            # ax.set_title("Histogram of Melting Temperature (Tm)")
            # ax.set_xlabel("Tm (°C)")
            # ax.set_ylabel("Count")
            # st.pyplot(fig)

            # st.subheader("GC% vs Melting Temperature (Tm)")
            # fig5, ax5 = plt.subplots()
            # sns.regplot(data=df, x="GC%", y="Tm", ax=ax5, scatter_kws={"s": 50})
            # ax5.set_title("GC% vs Tm with Regression Line")
            # st.pyplot(fig5)

            # st.subheader("ASO Positions on Target Exon")

            # fig6 = gen.plot_chromosome_map(df)
            # st.pyplot(fig6)

            # st.subheader("🔎 Filtering Pipeline Overview")

            # fig7 = gen.plot_filtering_barplot()
            # st.pyplot(fig7)

            # Построение графиков в две колонки
            col1, col2 = st.columns(2)

            with col1:
                st.subheader("GC% Distribution")
                fig1, ax1 = plt.subplots()
                df["GC%"].plot(kind="hist", bins=15, edgecolor="black", alpha=0.7, ax=ax1)
                ax1.axvline(df["GC%"].mean(), color='red', linestyle='dashed', linewidth=1)
                ax1.set_title("Distribution of GC% in ASO Candidates")
                ax1.set_xlabel("GC%")
                ax1.set_ylabel("Frequency")
                st.pyplot(fig1)

            with col2:
                st.subheader("ASO Length Boxplot")
                fig2, ax2 = plt.subplots()
                ax2.boxplot(df["Length"], vert=False)
                ax2.set_title("Boxplot of ASO Length")
                ax2.set_xlabel("Length (nt)")
                st.pyplot(fig2)

            col3, col4 = st.columns(2)

            with col3:
                st.subheader("Tm Distribution of ASOs")
                fig3, ax3 = plt.subplots()
                ax3.hist(df["Tm"], bins=20, color='skyblue', edgecolor='black')
                ax3.set_title("Histogram of Melting Temperature (Tm)")
                ax3.set_xlabel("Tm (°C)")
                ax3.set_ylabel("Count")
                st.pyplot(fig3)

            with col4:
                st.subheader("GC% vs Melting Temperature (Tm)")
                fig4, ax4 = plt.subplots()
                sns.regplot(data=df, x="GC%", y="Tm", ax=ax4, scatter_kws={"s": 50})
                ax4.set_title("GC% vs Tm with Regression Line")
                st.pyplot(fig4)

            col5, col6 = st.columns(2)

            with col5:
                st.subheader("ASO Positions on Target Exon")
                fig5 = gen.plot_chromosome_map(df)
                st.pyplot(fig5)

            with col6:
                st.subheader("Filtering Pipeline Overview")
                fig6 = gen.plot_filtering_barplot()
                st.pyplot(fig6)



        except Exception as e:
            st.error(f"An error occurred:\n{e}")