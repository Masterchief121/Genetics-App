import tkinter as tk
from tkinter import ttk, messagebox
import random
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np


class GeneVizApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GeneViz")
        self.root.geometry("480x600")
        self.root.configure(bg="#121212")

        # Create a header
        header = tk.Frame(root, bg="#121212")
        header.pack(fill="x", pady=20)

        title = tk.Label(header, text="GeneViz", font=("Arial", 24, "bold"), fg="white", bg="#121212")
        title.pack()

        # Create a frame for app buttons (home screen)
        self.home_screen = tk.Frame(root, bg="#121212")
        self.home_screen.pack(fill="both", expand=True, padx=20, pady=20)

        # Create app view frame (for individual app screens)
        self.app_view = tk.Frame(root, bg="#121212")

        # Create app buttons
        self.create_app_button(self.home_screen, "DNA Visualizer", self.open_dna_visualizer, "#4e54c8", 0, 0)
        self.create_app_button(self.home_screen, "Inheritance", self.open_inheritance_simulator, "#16a085", 0, 1)
        self.create_app_button(self.home_screen, "Evolution", self.open_evolution_simulator, "#e74c3c", 1, 0)
        self.create_app_button(self.home_screen, "DNA Translator", self.open_dna_translator, "#f39c12", 1, 1)

        # Add a status message at the bottom
        self.status = tk.Label(root, text="GeneViz - A Genetics Visualization App", fg="white", bg="#1C1C1E")
        self.status.pack(side="bottom", fill="x", pady=10)

    def create_app_button(self, parent, text, command, bg_color, row, col):
        parent.grid_columnconfigure(0, weight=1)
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)
        parent.grid_rowconfigure(1, weight=1)

        button = tk.Button(
            parent,
            text=text,
            font=("Arial", 14),
            bg=bg_color,
            fg="white",
            padx=20,
            pady=30,
            bd=0,
            command=command
        )
        button.grid(row=row, column=col, padx=10, pady=10, sticky="nsew")
        return button

    def go_home(self):
        self.app_view.pack_forget()
        self.home_screen.pack(fill="both", expand=True, padx=20, pady=20)
        self.status.config(text="GeneViz - A Genetics Visualization App")

    def open_app_screen(self, title, status_text):
        # Hide home screen
        self.home_screen.pack_forget()

        # Clear app view
        for widget in self.app_view.winfo_children():
            widget.destroy()

        # Configure and show app view
        self.app_view.pack(fill="both", expand=True, padx=20, pady=20)
        self.status.config(text=status_text)

        # Add title and back button
        header = tk.Frame(self.app_view, bg="#121212")
        header.pack(fill="x", pady=(0, 20))

        back_btn = tk.Button(header, text="← Back", bg="#121212", fg="#4e54c8", bd=0, command=self.go_home)
        back_btn.pack(side="left")

        title_label = tk.Label(header, text=title, font=("Arial", 18, "bold"), fg="white", bg="#121212")
        title_label.pack(pady=10)

        return self.app_view

    def open_dna_visualizer(self):
        app_view = self.open_app_screen("DNA Visualizer", "DNA Visualizer")

        # Create DNA input controls
        input_frame = tk.Frame(app_view, bg="#1C1C1E", padx=10, pady=10)
        input_frame.pack(fill="x", pady=10)

        tk.Label(input_frame, text="DNA Sequence:", fg="white", bg="#1C1C1E", font=("Arial", 12)).pack(side="left",
                                                                                                       padx=(0, 10))

        self.dna_seq_var = tk.StringVar(value="ATGC")
        entry = tk.Entry(input_frame, textvariable=self.dna_seq_var, width=20, font=("Arial", 12))
        entry.pack(side="left", fill="x", expand=True)

        # Buttons
        button_frame = tk.Frame(app_view, bg="#121212")
        button_frame.pack(fill="x", pady=10)

        random_btn = tk.Button(
            button_frame, text="Random DNA", bg="#4e54c8", fg="white", bd=0, padx=10, pady=5,
            command=self.generate_random_dna
        )
        random_btn.pack(side="left", padx=10)

        visualize_btn = tk.Button(
            button_frame, text="Visualize", bg="#16a085", fg="white", bd=0, padx=10, pady=5,
            command=self.visualize_dna
        )
        visualize_btn.pack(side="left", padx=10)

        # Frame for DNA visualization
        self.viz_frame = tk.Frame(app_view, bg="#121212")
        self.viz_frame.pack(fill="both", expand=True, pady=10)

        # Initial message
        tk.Label(
            self.viz_frame,
            text="Enter a DNA sequence above\nand press 'Visualize'",
            font=("Arial", 14),
            fg="#AAAAAA",
            bg="#121212"
        ).pack(pady=50)

    def generate_random_dna(self):
        length = random.randint(5, 15)
        bases = ['A', 'T', 'G', 'C']
        sequence = ''.join(random.choice(bases) for _ in range(length))
        self.dna_seq_var.set(sequence)

    def visualize_dna(self):
        # Clear visualization frame
        for widget in self.viz_frame.winfo_children():
            widget.destroy()

        sequence = self.dna_seq_var.get().upper()

        # Validate sequence
        if not all(base in "ATGC" for base in sequence):
            messagebox.showerror("Error", "DNA sequence can only contain A, T, G, C")
            return

        # Create figure with appropriate styling
        plt.rcParams.update({
            'figure.facecolor': '#121212',
            'axes.facecolor': '#121212',
            'axes.edgecolor': '#666666',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'text.color': 'white',
            'figure.figsize': (4, 3)
        })

        fig, ax = plt.subplots()

        # Define colors for each base
        colors = {'A': '#FF5757', 'T': '#57FF57', 'G': '#5757FF', 'C': '#FFFF57'}

        # Plot the DNA sequence
        for i, base in enumerate(sequence):
            # Top strand
            ax.add_patch(plt.Rectangle((i, 0), 0.8, 0.8, color=colors[base], alpha=0.8))
            ax.text(i + 0.4, 0.4, base, ha='center', va='center', fontweight='bold', fontsize=12, color='black')

            # Bottom strand (complementary)
            complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}[base]
            ax.add_patch(plt.Rectangle((i, -1), 0.8, 0.8, color=colors[complement], alpha=0.8))
            ax.text(i + 0.4, -0.6, complement, ha='center', va='center', fontweight='bold', fontsize=12, color='black')

            # Connecting lines
            ax.plot([i + 0.4, i + 0.4], [0, -0.2], 'white', alpha=0.7)

        # Style the plot
        ax.set_xlim(-0.5, len(sequence) + 0.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_title("DNA Double Helix", color='white', fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Embed plot in the frame
        canvas = FigureCanvasTkAgg(fig, master=self.viz_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        # Add analysis stats
        stats_frame = tk.Frame(self.viz_frame, bg="#2C2C2E", padx=10, pady=10)
        stats_frame.pack(fill="x", side="bottom")

        gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100

        tk.Label(
            stats_frame,
            text=f"Length: {len(sequence)} bp | GC: {gc_content:.1f}% | AT: {100 - gc_content:.1f}%",
            font=("Arial", 12),
            fg="white",
            bg="#2C2C2E"
        ).pack()

    def open_inheritance_simulator(self):
        app_view = self.open_app_screen("Inheritance Simulator", "Inheritance Simulator")

        # Create controls
        control_frame = tk.Frame(app_view, bg="#1C1C1E", padx=10, pady=15)
        control_frame.pack(fill="x", pady=10)

        # First row - trait name
        tk.Label(
            control_frame,
            text="Trait:",
            font=("Arial", 14),
            fg="white",
            bg="#1C1C1E"
        ).grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.trait_var = tk.StringVar(value="Eye Color")
        trait_entry = tk.Entry(
            control_frame,
            textvariable=self.trait_var,
            font=("Arial", 14),
            bg="#2C2C2E",
            fg="white",
            bd=0,
            width=15
        )
        trait_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        # Second row - parent genotypes
        tk.Label(
            control_frame,
            text="Parent 1:",
            font=("Arial", 14),
            fg="white",
            bg="#1C1C1E"
        ).grid(row=1, column=0, padx=5, pady=5, sticky="w")

        self.parent1_var = tk.StringVar(value="Aa")
        p1_entry = tk.Entry(
            control_frame,
            textvariable=self.parent1_var,
            font=("Arial", 14),
            bg="#2C2C2E",
            fg="white",
            bd=0,
            width=5
        )
        p1_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        tk.Label(
            control_frame,
            text="Parent 2:",
            font=("Arial", 14),
            fg="white",
            bg="#1C1C1E"
        ).grid(row=1, column=2, padx=5, pady=5, sticky="w")

        self.parent2_var = tk.StringVar(value="Aa")
        p2_entry = tk.Entry(
            control_frame,
            textvariable=self.parent2_var,
            font=("Arial", 14),
            bg="#2C2C2E",
            fg="white",
            bd=0,
            width=5
        )
        p2_entry.grid(row=1, column=3, padx=5, pady=5, sticky="w")

        # Punnett square button
        punnett_button = tk.Button(
            app_view,
            text="Generate Punnett Square",
            font=("Arial", 14),
            fg="white",
            bg="#16a085",
            padx=20,
            pady=10,
            bd=0,
            command=self.generate_punnett
        )
        punnett_button.pack(pady=10)

        # Frame for visualization
        self.punnett_frame = tk.Frame(app_view, bg="#121212")
        self.punnett_frame.pack(fill="both", expand=True, pady=10)

        # Initial instruction
        tk.Label(
            self.punnett_frame,
            text="Enter parent genotypes above\nand press 'Generate Punnett Square'",
            font=("Arial", 14),
            fg="#AAAAAA",
            bg="#121212"
        ).pack(pady=50)

    def generate_punnett(self):
        # Clear visualization frame
        for widget in self.punnett_frame.winfo_children():
            widget.destroy()

        p1 = self.parent1_var.get()
        p2 = self.parent2_var.get()
        trait = self.trait_var.get()

        # Validate input
        if len(p1) != 2 or len(p2) != 2 or not all(c.isalpha() for c in p1 + p2):
            messagebox.showerror("Error", "Genotypes must be two letters (e.g., 'Aa')")
            return

        # Create gametes
        p1_gametes = [p1[0], p1[1]]
        p2_gametes = [p2[0], p2[1]]

        # Create figure with dark background
        plt.rcParams.update({
            'figure.facecolor': '#121212',
            'axes.facecolor': '#121212',
            'axes.edgecolor': '#666666',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'text.color': 'white',
            'figure.figsize': (4, 4)
        })

        fig, ax = plt.subplots()

        # Draw punnett square
        for i, p1g in enumerate(p1_gametes):
            for j, p2g in enumerate(p2_gametes):
                # Combine gametes
                genotype = ''.join(sorted([p1g, p2g], key=lambda x: (x.lower(), x)))

                # Determine phenotype color
                if genotype[0].isupper() or genotype[1].isupper():
                    color = '#4e54c8'  # dominant
                else:
                    color = '#e74c3c'  # recessive

                # Add rectangle
                rect = plt.Rectangle((j, -i - 1), 1, 1, color=color, alpha=0.7)
                ax.add_patch(rect)

                # Add text
                ax.text(j + 0.5, -i - 0.5, genotype, ha='center', va='center',
                        fontweight='bold', color='white', fontsize=14)

        # Add gamete labels with white text
        for i, g in enumerate(p1_gametes):
            ax.text(-0.3, -i - 0.5, g, ha='center', va='center',
                    fontweight='bold', color='white', fontsize=14)

        for j, g in enumerate(p2_gametes):
            ax.text(j + 0.5, 0.3, g, ha='center', va='center',
                    fontweight='bold', color='white', fontsize=14)

        # Set limits and title
        ax.set_xlim(-0.5, 2)
        ax.set_ylim(-3, 1)
        ax.set_title(f"Punnett Square: {trait}", color='white', fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])

        # Add parent labels
        ax.text(-0.3, 0.3, "Parent 1", ha='center', va='center',
                fontweight='bold', color='white', fontsize=12)
        ax.text(1, 0.7, "Parent 2", ha='center', va='center',
                fontweight='bold', color='white', fontsize=12)

        # Embed plot in the frame
        canvas = FigureCanvasTkAgg(fig, master=self.punnett_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        # Calculate ratios
        genotypes = {}
        phenotypes = {'Dominant': 0, 'Recessive': 0}

        for p1g in p1_gametes:
            for p2g in p2_gametes:
                genotype = ''.join(sorted([p1g, p2g], key=lambda x: (x.lower(), x)))
                genotypes[genotype] = genotypes.get(genotype, 0) + 1

                # Count phenotypes
                if genotype[0].isupper() or genotype[1].isupper():
                    phenotypes['Dominant'] += 1
                else:
                    phenotypes['Recessive'] += 1

        # Show results
        results_frame = tk.Frame(self.punnett_frame, bg="#2C2C2E", padx=10, pady=10)
        results_frame.pack(fill="x", side="bottom")

        results_text = f"Genotype Ratio: "
        for genotype, count in genotypes.items():
            results_text += f"{genotype}:{count}/4  "

        results_text += f"\nPhenotype Ratio: {phenotypes['Dominant']}:{phenotypes['Recessive']} ({phenotypes['Dominant'] / 4 * 100:.0f}%:{phenotypes['Recessive'] / 4 * 100:.0f}%)"

        tk.Label(
            results_frame,
            text=results_text,
            font=("Arial", 12),
            fg="white",
            bg="#2C2C2E"
        ).pack()

    def open_evolution_simulator(self):
        app_view = self.open_app_screen("Evolution Simulator", "Evolution Simulator")

        # Create controls
        control_frame = tk.Frame(app_view, bg="#1C1C1E", padx=10, pady=15)
        control_frame.pack(fill="x", pady=10)

        # Sliders for parameters
        tk.Label(
            control_frame,
            text="Initial Allele Frequency:",
            font=("Arial", 14),
            fg="white",
            bg="#1C1C1E"
        ).grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.p_var = tk.DoubleVar(value=0.6)
        p_scale = tk.Scale(
            control_frame,
            variable=self.p_var,
            from_=0.0, to=1.0,
            resolution=0.01,
            orient=tk.HORIZONTAL,
            bg="#2C2C2E",
            fg="white",
            highlightthickness=0,
            troughcolor="#666666",
            bd=0,
            width=15,
            length=150
        )
        p_scale.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        tk.Label(
            control_frame,
            text="Selection Strength:",
            font=("Arial", 14),
            fg="white",
            bg="#1C1C1E"
        ).grid(row=1, column=0, padx=5, pady=5, sticky="w")

        self.s_var = tk.DoubleVar(value=0.2)
        s_scale = tk.Scale(
            control_frame,
            variable=self.s_var,
            from_=0.0, to=1.0,
            resolution=0.01,
            orient=tk.HORIZONTAL,
            bg="#2C2C2E",
            fg="white",
            highlightthickness=0,
            troughcolor="#666666",
            bd=0,
            width=15,
            length=150
        )
        s_scale.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        tk.Label(
            control_frame,
            text="Generations:",
            font=("Arial", 14),
            fg="white",
            bg="#1C1C1E"
        ).grid(row=2, column=0, padx=5, pady=5, sticky="w")

        self.gen_var = tk.IntVar(value=10)
        gen_scale = tk.Scale(
            control_frame,
            variable=self.gen_var,
            from_=5, to=20,
            resolution=1,
            orient=tk.HORIZONTAL,
            bg="#2C2C2E",
            fg="white",
            highlightthickness=0,
            troughcolor="#666666",
            bd=0,
            width=15,
            length=150
        )
        gen_scale.grid(row=2, column=1, padx=5, pady=5, sticky="w")

        # Simulate button
        simulate_button = tk.Button(
            app_view,
            text="Simulate Evolution",
            font=("Arial", 14),
            fg="white",
            bg="#e74c3c",
            padx=20,
            pady=10,
            bd=0,
            command=self.simulate_evolution
        )
        simulate_button.pack(pady=10)

        # Frame for visualization
        self.evolution_frame = tk.Frame(app_view, bg="#121212")
        self.evolution_frame.pack(fill="both", expand=True, pady=10)

        # Initial instruction
        tk.Label(
            self.evolution_frame,
            text="Adjust the parameters above\nand press 'Simulate Evolution'",
            font=("Arial", 14),
            fg="#AAAAAA",
            bg="#121212"
        ).pack(pady=50)

    def simulate_evolution(self):
        # Clear visualization frame
        for widget in self.evolution_frame.winfo_children():
            widget.destroy()

        p0 = self.p_var.get()  # Initial frequency of dominant allele A
        s = self.s_var.get()  # Selection coefficient
        generations = self.gen_var.get()

        # Initialize arrays
        p_values = [p0]
        q_values = [1 - p0]

        # Simulate selection over generations
        p_current = p0
        for gen in range(1, generations + 1):
            q_current = 1 - p_current

            # Calculate genotype frequencies
            AA = p_current ** 2
            Aa = 2 * p_current * q_current
            aa = q_current ** 2

            # Calculate mean fitness
            w_bar = AA + Aa + aa * (1 - s)

            # Calculate allele frequencies after selection
            p_next = (AA + Aa / 2) / w_bar

            # Store values
            p_values.append(p_next)
            q_values.append(1 - p_next)

            # Update for next generation
            p_current = p_next

        # Create figure with dark background
        plt.rcParams.update({
            'figure.facecolor': '#121212',
            'axes.facecolor': '#121212',
            'axes.edgecolor': '#666666',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'text.color': 'white',
            'figure.figsize': (4, 4)
        })

        fig, ax = plt.subplots()

        # Plot allele frequencies over time
        gens = range(generations + 1)
        ax.plot(gens, p_values, 'b-', color='#4e54c8', linewidth=3, label='Dominant (A)')
        ax.plot(gens, q_values, 'r-', color='#e74c3c', linewidth=3, label='Recessive (a)')

        ax.set_xlabel('Generation')
        ax.set_ylabel('Allele Frequency')
        ax.set_title('Allele Frequencies Over Time')
        ax.legend()
        ax.grid(True, alpha=0.3, linestyle='--', color='#555555')

        # Embed plot in the frame
        canvas = FigureCanvasTkAgg(fig, master=self.evolution_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

        # Show summary
        summary_frame = tk.Frame(self.evolution_frame, bg="#2C2C2E", padx=10, pady=10)
        summary_frame.pack(fill="x", side="bottom")

        summary_text = f"Selection coefficient: {s:.2f}\n"
        summary_text += f"Initial: A={p0:.2f}, a={1 - p0:.2f} | Final: A={p_values[-1]:.2f}, a={q_values[-1]:.2f}"

        tk.Label(
            summary_frame,
            text=summary_text,
            font=("Arial", 12),
            fg="white",
            bg="#2C2C2E"
        ).pack()

    def open_dna_translator(self):
        app_view = self.open_app_screen("DNA Translator", "DNA Translator")

        # Create DNA input controls
        input_frame = tk.Frame(app_view, bg="#1C1C1E", padx=10, pady=10)
        input_frame.pack(fill="x", pady=10)

        tk.Label(
            input_frame,
            text="Enter DNA Sequence:",
            font=("Arial", 14),
            fg="white",
            bg="#1C1C1E"
        ).pack(anchor="w", pady=(0, 5))

        self.translator_seq_var = tk.StringVar(value="ATGCCATAG")
        entry = tk.Entry(
            input_frame,
            textvariable=self.translator_seq_var,
            font=("Arial", 14),
            bg="#2C2C2E",
            fg="white",
            bd=0
        )
        entry.pack(fill="x", pady=5)

        # Button frame
        button_frame = tk.Frame(app_view, bg="#121212")
        button_frame.pack(fill="x", pady=10)

        # Random DNA button
        random_btn = tk.Button(
            button_frame,
            text="Generate Random",
            font=("Arial", 12),
            fg="white",
            bg="#4e54c8",
            padx=10,
            pady=5,
            bd=0,
            command=self.generate_random_codon
        )
        random_btn.pack(side="left", padx=10)

        # Translate button
        translate_btn = tk.Button(
            button_frame,
            text="Translate DNA",
            font=("Arial", 12),
            fg="white",
            bg="#f39c12",
            padx=10,
            pady=5,
            bd=0,
            command=self.translate_dna
        )
        translate_btn.pack(side="left", padx=10)

        # Frame for translation results
        self.translation_frame = tk.Frame(app_view, bg="#121212")
        self.translation_frame.pack(fill="both", expand=True, pady=10)

        # Initial instruction
        tk.Label(
            self.translation_frame,
            text="Enter a DNA sequence above\nand press 'Translate DNA'",
            font=("Arial", 14),
            fg="#AAAAAA",
            bg="#121212"
        ).pack(pady=50)

    def generate_random_codon(self):
        # Generate a DNA sequence that's a multiple of 3 (for clean translation)
        codon_count = random.randint(2, 5)
        bases = ['A', 'T', 'G', 'C']
        sequence = ''.join(random.choice(bases) for _ in range(codon_count * 3))
        self.translator_seq_var.set(sequence)

    def translate_dna(self):
        # Clear translation frame
        for widget in self.translation_frame.winfo_children():
            widget.destroy()

        dna_sequence = self.translator_seq_var.get().upper()

        # Validate sequence
        if not all(base in "ATGC" for base in dna_sequence):
            messagebox.showerror("Error", "DNA sequence can only contain A, T, G, C")
            return

        # Convert DNA to RNA (replace T with U)
        rna_sequence = dna_sequence.replace('T', 'U')

        # Create a results panel
        results_panel = tk.Frame(self.translation_frame, bg="#1C1C1E", padx=15, pady=15)
        results_panel.pack(fill="both", expand=True, padx=10, pady=10)

        # DNA section
        tk.Label(
            results_panel,
            text="DNA:",
            font=("Arial", 14, "bold"),
            fg="#4e54c8",
            bg="#1C1C1E"
        ).pack(anchor="w", pady=(0, 5))

        # Create DNA display frame
        dna_display = tk.Frame(results_panel, bg="#2C2C2E", padx=10, pady=10)
        dna_display.pack(fill="x", pady=(0, 15))

        # Color-coded DNA display
        dna_colors = {'A': '#FF5757', 'T': '#57FF57', 'G': '#5757FF', 'C': '#FFFF57'}
        dna_frame = tk.Frame(dna_display, bg="#2C2C2E")
        dna_frame.pack()

        for i, base in enumerate(dna_sequence):
            base_frame = tk.Frame(dna_frame, bg=dna_colors[base], width=30, height=30)
            base_frame.grid(row=0, column=i, padx=2)
            tk.Label(
                base_frame,
                text=base,
                font=("Arial", 14, "bold"),
                fg="black",
                bg=dna_colors[base]
            ).pack(padx=5, pady=5)

            # Add codon separators
            if (i + 1) % 3 == 0 and i < len(dna_sequence) - 1:
                separator = tk.Frame(dna_frame, bg="#1C1C1E", width=1, height=30)
                separator.grid(row=0, column=i + 1)

        # RNA section
        tk.Label(
            results_panel,
            text="RNA:",
            font=("Arial", 14, "bold"),
            fg="#e74c3c",
            bg="#1C1C1E"
        ).pack(anchor="w", pady=(0, 5))

        rna_display = tk.Frame(results_panel, bg="#2C2C2E", padx=10, pady=10)
        rna_display.pack(fill="x", pady=(0, 15))

        # Color-coded RNA display
        rna_colors = {'A': '#FF5757', 'U': '#FFA500', 'G': '#5757FF', 'C': '#FFFF57'}
        rna_frame = tk.Frame(rna_display, bg="#2C2C2E")
        rna_frame.pack()

        for i, base in enumerate(rna_sequence):
            base_frame = tk.Frame(rna_frame, bg=rna_colors[base], width=30, height=30)
            base_frame.grid(row=0, column=i, padx=2)
            tk.Label(
                base_frame,
                text=base,
                font=("Arial", 14, "bold"),
                fg="black",
                bg=rna_colors[base]
            ).pack(padx=5, pady=5)

            # Add codon separators
            if (i + 1) % 3 == 0 and i < len(rna_sequence) - 1:
                separator = tk.Frame(rna_frame, bg="#1C1C1E", width=1, height=30)
                separator.grid(row=0, column=i + 0.5)

        # Translate to protein if length is multiple of 3
        if len(rna_sequence) % 3 == 0 and len(rna_sequence) > 0:
            # Codon to amino acid mapping
            codon_table = {
                'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
                'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
                'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
                'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
                'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
                'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
            }

            # Protein section
            tk.Label(
                results_panel,
                text="Protein:",
                font=("Arial", 14, "bold"),
                fg="#16a085",
                bg="#1C1C1E"
            ).pack(anchor="w", pady=(0, 5))

            protein_display = tk.Frame(results_panel, bg="#2C2C2E", padx=10, pady=10)
            protein_display.pack(fill="x")

            protein = ""
            codons = []

            # Extract codons and translate
            for i in range(0, len(rna_sequence), 3):
                codon = rna_sequence[i:i + 3]
                if len(codon) == 3:
                    codons.append(codon)
                    amino_acid = codon_table.get(codon, '?')
                    protein += amino_acid

            # Amino acid color coding by property
            aa_colors = {
                'A': '#8FBC8F', 'G': '#8FBC8F', 'V': '#8FBC8F', 'L': '#8FBC8F', 'I': '#8FBC8F', 'M': '#8FBC8F',
                # Hydrophobic
                'F': '#87CEEB', 'W': '#87CEEB', 'Y': '#87CEEB',  # Aromatic
                'S': '#FAFAD2', 'T': '#FAFAD2', 'N': '#FAFAD2', 'Q': '#FAFAD2',  # Polar
                'K': '#FFC0CB', 'R': '#FFC0CB', 'H': '#FFC0CB',  # Basic
                'D': '#F08080', 'E': '#F08080',  # Acidic
                'P': '#FFD700', 'C': '#FFD700',  # Special
                '*': '#DCDCDC'  # Stop codon
            }

            # Display each codon and its corresponding amino acid
            codon_aa_frame = tk.Frame(protein_display, bg="#2C2C2E")
            codon_aa_frame.pack()

            for i, (codon, aa) in enumerate(zip(codons, protein)):
                # Create codon display
                codon_frame = tk.Frame(codon_aa_frame, bg="#3A3A3C", padx=5, pady=5)
                codon_frame.grid(row=0, column=i, padx=5)

                tk.Label(
                    codon_frame,
                    text=codon,
                    font=("Arial", 12),
                    fg="white",
                    bg="#3A3A3C"
                ).pack()

                # Create arrow
                arrow_frame = tk.Frame(codon_aa_frame, bg="#2C2C2E")
                arrow_frame.grid(row=1, column=i)

                tk.Label(
                    arrow_frame,
                    text="↓",
                    font=("Arial", 14),
                    fg="white",
                    bg="#2C2C2E"
                ).pack(pady=2)

                # Create amino acid display
                aa_color = aa_colors.get(aa, '#FFFFFF')
                aa_frame = tk.Frame(codon_aa_frame, bg=aa_color, padx=8, pady=8)
                aa_frame.grid(row=2, column=i, padx=5)

                tk.Label(
                    aa_frame,
                    text=aa,
                    font=("Arial", 14, "bold"),
                    fg="black",
                    bg=aa_color
                ).pack()

                # Add amino acid name below
                aa_names = {
                    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
                    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
                    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
                    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
                    '*': 'STOP'
                }

                tk.Label(
                    codon_aa_frame,
                    text=aa_names.get(aa, '?'),
                    font=("Arial", 10),
                    fg="#AAAAAA",
                    bg="#2C2C2E"
                ).grid(row=3, column=i, pady=(2, 0))
        else:
            # If not a multiple of 3, show a message
            tk.Label(
                results_panel,
                text="Translation requires a DNA sequence with length divisible by 3",
                font=("Arial", 12),
                fg="#AAAAAA",
                bg="#1C1C1E"
            ).pack(pady=20)


if __name__ == "__main__":
    print("Starting GeneViz app...")
    root = tk.Tk()
    app = GeneVizApp(root)
    print("App initialized, entering mainloop...")
    root.mainloop()
    print("App closed")