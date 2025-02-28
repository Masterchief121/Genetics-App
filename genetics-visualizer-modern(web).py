import streamlit as st
import random
import matplotlib.pyplot as plt
import numpy as np

# Streamlit app setup
st.set_page_config(page_title="GeneViz", layout="wide")
plt.style.use('dark_background')

# Sidebar for navigation
st.sidebar.title("GeneViz")
app_choice = st.sidebar.radio(
    "Select Tool",
    ["Home", "DNA Visualizer", "Inheritance Simulator", "Evolution Simulator", "DNA Translator"]
)

# Status message
status_placeholder = st.empty()
status_placeholder.markdown("GeneViz - A Genetics Visualization App")

# Main content
if app_choice == "Home":
    st.title("Welcome to GeneViz")
    st.markdown("""
    A Genetics Visualization App. Select a tool from the sidebar to begin:
    - **DNA Visualizer**: Visualize DNA double helix structure.
    - **Inheritance Simulator**: Generate Punnett squares for genetic crosses.
    - **Evolution Simulator**: Simulate allele frequency changes over generations.
    - **DNA Translator**: Translate DNA sequences into proteins.
    """)
    status_placeholder.markdown("GeneViz - A Genetics Visualization App")

elif app_choice == "DNA Visualizer":
    st.header("DNA Visualizer")
    status_placeholder.markdown("DNA Visualizer")

    # DNA input with key for state management
    if "dna_seq" not in st.session_state:
        st.session_state.dna_seq = "ATGC"
    dna_seq = st.text_input("DNA Sequence", value=st.session_state.dna_seq, key="dna_seq_input")

    col1, col2 = st.columns(2)
    with col1:
        if st.button("Random DNA", key="random_dna"):
            length = random.randint(5, 15)
            bases = ['A', 'T', 'G', 'C']
            st.session_state.dna_seq = ''.join(random.choice(bases) for _ in range(length))
            # No need to rerun manually; key syncs it
    with col2:
        visualize = st.button("Visualize", key="visualize_dna")

    if visualize and dna_seq:
        sequence = dna_seq.upper()
        if not all(base in "ATGC" for base in sequence):
            st.error("DNA sequence can only contain A, T, G, C")
        else:
            # Plotting logic (unchanged)
            fig, ax = plt.subplots(figsize=(4, 3))
            colors = {'A': '#FF5757', 'T': '#57FF57', 'G': '#5757FF', 'C': '#FFFF57'}
            for i, base in enumerate(sequence):
                ax.add_patch(plt.Rectangle((i, 0), 0.8, 0.8, color=colors[base], alpha=0.8))
                ax.text(i + 0.4, 0.4, base, ha='center', va='center', fontweight='bold', fontsize=12, color='black')
                complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[base]
                ax.add_patch(plt.Rectangle((i, -1), 0.8, 0.8, color=colors[complement], alpha=0.8))
                ax.text(i + 0.4, -0.6, complement, ha='center', va='center', fontweight='bold', fontsize=12, color='black')
                ax.plot([i + 0.4, i + 0.4], [0, -0.2], 'white', alpha=0.7)
            ax.set_xlim(-0.5, len(sequence) + 0.5)
            ax.set_ylim(-1.5, 1.5)
            ax.set_title("DNA Double Helix", color='white', fontsize=14)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            st.pyplot(fig)

            gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
            st.markdown(f"**Length:** {len(sequence)} bp | **GC:** {gc_content:.1f}% | **AT:** {100 - gc_content:.1f}%")

elif app_choice == "Inheritance Simulator":
    # Unchanged from previous version
    st.header("Inheritance Simulator")
    status_placeholder.markdown("Inheritance Simulator")
    trait = st.text_input("Trait", value="Eye Color")
    col1, col2 = st.columns(2)
    with col1:
        parent1 = st.text_input("Parent 1 Genotype", value="Aa")
    with col2:
        parent2 = st.text_input("Parent 2 Genotype", value="Aa")
    if st.button("Generate Punnett Square"):
        p1, p2 = parent1, parent2
        if len(p1) != 2 or len(p2) != 2 or not all(c.isalpha() for c in p1 + p2):
            st.error("Genotypes must be two letters (e.g., 'Aa')")
        else:
            p1_gametes = [p1[0], p1[1]]
            p2_gametes = [p2[0], p2[1]]
            fig, ax = plt.subplots(figsize=(4, 4))
            for i, p1g in enumerate(p1_gametes):
                for j, p2g in enumerate(p2_gametes):
                    genotype = ''.join(sorted([p1g, p2g], key=lambda x: (x.lower(), x)))
                    color = '#4e54c8' if genotype[0].isupper() or genotype[1].isupper() else '#e74c3c'
                    rect = plt.Rectangle((j, -i - 1), 1, 1, color=color, alpha=0.7)
                    ax.add_patch(rect)
                    ax.text(j + 0.5, -i - 0.5, genotype, ha='center', va='center', fontweight='bold', color='white', fontsize=14)
            for i, g in enumerate(p1_gametes):
                ax.text(-0.3, -i - 0.5, g, ha='center', va='center', fontweight='bold', color='white', fontsize=14)
            for j, g in enumerate(p2_gametes):
                ax.text(j + 0.5, 0.3, g, ha='center', va='center', fontweight='bold', color='white', fontsize=14)
            ax.set_xlim(-0.5, 2)
            ax.set_ylim(-3, 1)
            ax.set_title(f"Punnett Square: {trait}", color='white', fontsize=14)
            ax.set_xticks([])
            ax.set_yticks([])
            st.pyplot(fig)
            genotypes = {}
            phenotypes = {'Dominant': 0, 'Recessive': 0}
            for p1g in p1_gametes:
                for p2g in p2_gametes:
                    genotype = ''.join(sorted([p1g, p2g], key=lambda x: (x.lower(), x)))
                    genotypes[genotype] = genotypes.get(genotype, 0) + 1
                    if genotype[0].isupper() or genotype[1].isupper():
                        phenotypes['Dominant'] += 1
                    else:
                        phenotypes['Recessive'] += 1
            st.markdown(f"**Genotype Ratio:** {'  '.join(f'{g}:{c}/4' for g, c in genotypes.items())}")
            st.markdown(f"**Phenotype Ratio:** {phenotypes['Dominant']}:{phenotypes['Recessive']} ({phenotypes['Dominant']/4*100:.0f}%:{phenotypes['Recessive']/4*100:.0f}%)")

elif app_choice == "Evolution Simulator":
    # Unchanged from previous version
    st.header("Evolution Simulator")
    status_placeholder.markdown("Evolution Simulator")
    p0 = st.slider("Initial Allele Frequency (A)", 0.0, 1.0, 0.6, 0.01)
    s = st.slider("Selection Strength", 0.0, 1.0, 0.2, 0.01)
    generations = st.slider("Generations", 5, 20, 10, 1)
    if st.button("Simulate Evolution"):
        p_values = [p0]
        q_values = [1 - p0]
        p_current = p0
        for _ in range(1, generations + 1):
            q_current = 1 - p_current
            AA = p_current ** 2
            Aa = 2 * p_current * q_current
            aa = q_current ** 2
            w_bar = AA + Aa + aa * (1 - s)
            p_next = (AA + Aa / 2) / w_bar
            p_values.append(p_next)
            q_values.append(1 - p_next)
            p_current = p_next
        fig, ax = plt.subplots(figsize=(4, 4))
        gens = range(generations + 1)
        ax.plot(gens, p_values, color='#4e54c8', linewidth=3, label='Dominant (A)')
        ax.plot(gens, q_values, color='#e74c3c', linewidth=3, label='Recessive (a)')
        ax.set_xlabel('Generation')
        ax.set_ylabel('Allele Frequency')
        ax.set_title('Allele Frequencies Over Time')
        ax.legend()
        ax.grid(True, alpha=0.3, linestyle='--')
        st.pyplot(fig)
        st.markdown(f"**Selection Coefficient:** {s:.2f}")
        st.markdown(f"**Initial:** A={p0:.2f}, a={1-p0:.2f} | **Final:** A={p_values[-1]:.2f}, a={q_values[-1]:.2f}")

elif app_choice == "DNA Translator":
    st.header("DNA Translator")
    status_placeholder.markdown("DNA Translator")

    # DNA input with key for state management
    if "translator_seq" not in st.session_state:
        st.session_state.translator_seq = "ATGCCATAG"
    dna_seq = st.text_input("Enter DNA Sequence", value=st.session_state.translator_seq, key="translator_seq_input")

    col1, col2 = st.columns(2)
    with col1:
        if st.button("Generate Random", key="random_codon"):
            codon_count = random.randint(2, 5)
            bases = ['A', 'T', 'G', 'C']
            st.session_state.translator_seq = ''.join(random.choice(bases) for _ in range(codon_count * 3))
    with col2:
        translate = st.button("Translate DNA", key="translate_dna")

    if translate and dna_seq:
        dna_sequence = dna_seq.upper()
        if not all(base in "ATGC" for base in dna_sequence):
            st.error("DNA sequence can only contain A, T, G, C")
        else:
            rna_sequence = dna_sequence.replace('T', 'U')
            st.subheader("DNA:")
            dna_colors = {'A': '#FF5757', 'T': '#57FF57', 'G': '#5757FF', 'C': '#FFFF57'}
            dna_cols = st.columns(len(dna_sequence) + len(dna_sequence) // 3)
            for i, base in enumerate(dna_sequence):
                with dna_cols[i + i // 3]:
                    st.markdown(f"<div style='background-color:{dna_colors[base]};text-align:center;padding:10px;'>{base}</div>", unsafe_allow_html=True)
                if (i + 1) % 3 == 0 and i < len(dna_sequence) - 1:
                    with dna_cols[i + i // 3 + 1]:
                        st.markdown("<div style='width:5px'></div>", unsafe_allow_html=True)

            st.subheader("RNA:")
            rna_colors = {'A': '#FF5757', 'U': '#FFA500', 'G': '#5757FF', 'C': '#FFFF57'}
            rna_cols = st.columns(len(rna_sequence) + len(rna_sequence) // 3)
            for i, base in enumerate(rna_sequence):
                with rna_cols[i + i // 3]:
                    st.markdown(f"<div style='background-color:{rna_colors[base]};text-align:center;padding:10px;'>{base}</div>", unsafe_allow_html=True)
                if (i + 1) % 3 == 0 and i < len(rna_sequence) - 1:
                    with rna_cols[i + i // 3 + 1]:
                        st.markdown("<div style='width:5px'></div>", unsafe_allow_html=True)

            if len(rna_sequence) % 3 == 0 and len(rna_sequence) > 0:
                codon_table = {
                    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
                }
                aa_colors = {
                    'A': '#8FBC8F', 'G': '#8FBC8F', 'V': '#8FBC8F', 'L': '#8FBC8F', 'I': '#8FBC8F', 'M': '#8FBC8F',
                    'F': '#87CEEB', 'W': '#87CEEB', 'Y': '#87CEEB', 'S': '#FAFAD2', 'T': '#FAFAD2', 'N': '#FAFAD2', 'Q': '#FAFAD2',
                    'K': '#FFC0CB', 'R': '#FFC0CB', 'H': '#FFC0CB', 'D': '#F08080', 'E': '#F08080', 'P': '#FFD700', 'C': '#FFD700', '*': '#DCDCDC'
                }
                aa_names = {
                    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His',
                    'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
                    'Y': 'Tyr', 'V': 'Val', '*': 'STOP'
                }
                protein = ""
                codons = []
                for i in range(0, len(rna_sequence), 3):
                    codon = rna_sequence[i:i + 3]
                    if len(codon) == 3:
                        codons.append(codon)
                        protein += codon_table.get(codon, '?')
                st.subheader("Protein:")
                protein_cols = st.columns(len(protein) * 2 - 1)
                for i, (codon, aa) in enumerate(zip(codons, protein)):
                    with protein_cols[i * 2]:
                        st.markdown(f"<div style='text-align:center'>{codon}</div><div style='text-align:center'>↓</div><div style='background-color:{aa_colors.get(aa, '#FFFFFF')};text-align:center;padding:10px;'>{aa}</div><div style='text-align:center;color:#AAAAAA'>{aa_names.get(aa, '?')}</div>", unsafe_allow_html=True)
                    if i < len(protein) - 1:
                        with protein_cols[i * 2 + 1]:
                            st.markdown("<div style='width:10px'></div>", unsafe_allow_html=True)
            else:
                st.markdown("Translation requires a DNA sequence with length divisible by 3")
