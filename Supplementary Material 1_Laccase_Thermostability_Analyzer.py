from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from tkinter import Tk, filedialog, Label, Button, Text, Canvas, Frame, Scrollbar
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from Bio.SeqUtils import gc_fraction


class LaccaseThermostabilityAnalyzer:
    def __init__(self, master):
        self.master = master
        master.title("Laccase Thermostability Analyzer")

        self.create_ui()

    def create_ui(self):
        # Add a header
        header_frame = Frame(self.master, bg="#4CAF50", pady=10)
        header_frame.pack(fill="x")
        Label(header_frame, text="Laccase Thermostability Analyzer", font=("Helvetica", 18), fg="white", bg="#4CAF50").pack()

        # Create text areas for enzyme sequences
        self.enzyme1_text = Text(self.master, height=10, width=100)
        self.enzyme1_text.pack(pady=10)

        self.enzyme2_text = Text(self.master, height=10, width=100)
        self.enzyme2_text.pack(pady=10)

        # Buttons with modern design
        Button(self.master, text="Select Enzyme 1 File", command=self.select_enzyme1_file, bg="#2196F3", fg="white").pack()
        Button(self.master, text="Select Enzyme 2 File", command=self.select_enzyme2_file, bg="#2196F3", fg="white").pack()

        # Introduce a canvas for interactive visualization
        self.visualization_canvas = Canvas(self.master, bg="white", height=200)
        self.visualization_canvas.pack(pady=10)

        # Modern compare button
        compare_button = Button(self.master, text="Compare Thermostability", command=self.compare_thermostability,
                                bg="#4CAF50", fg="white")
        compare_button.pack()

        # Modern reset and clear buttons
        reset_button = Button(self.master, text="Reset", command=self.reset_input, bg="#FF9800", fg="white")
        reset_button.pack(side="left", padx=(0, 10))

        clear_button = Button(self.master, text="Clear", command=self.clear_results, bg="#FF5722", fg="white")
        clear_button.pack(side="left")

        # Result label with a touch of elegance
        self.result_frame = Frame(self.master)
        self.result_frame.pack()
        self.result_label = Text(self.result_frame, font=("Helvetica", 12), wrap="word", height=10, width=80)
        self.result_label.pack(side="left", fill="y")

        # Scrollbar for the result label
        result_scrollbar = Scrollbar(self.result_frame, command=self.result_label.yview)
        result_scrollbar.pack(side="right", fill="y")
        self.result_label.config(yscrollcommand=result_scrollbar.set)

        # Matplotlib figure and canvas references
        self.matplotlib_figure = None
        self.matplotlib_canvas = None

    def select_enzyme1_file(self):
        file_path = filedialog.askopenfilename(title="Select Enzyme 1 File", filetypes=[("FASTA files", "*.fasta;*.fas")])
        if file_path:
            with open(file_path, "r") as file:
                sequence = str(next(SeqIO.parse(file, "fasta")).seq)
                self.enzyme1_text.delete("1.0", "end")
                self.enzyme1_text.insert("end", sequence)
                self.result_label.insert("end", f"Enzyme 1 File Selected: {file_path}\n")

    def select_enzyme2_file(self):
        file_path = filedialog.askopenfilename(title="Select Enzyme 2 File", filetypes=[("FASTA files", "*.fasta;*.fas")])
        if file_path:
            with open(file_path, "r") as file:
                sequence = str(next(SeqIO.parse(file, "fasta")).seq)
                self.enzyme2_text.delete("1.0", "end")
                self.enzyme2_text.insert("end", sequence)
                self.result_label.insert("end", f"Enzyme 2 File Selected: {file_path}\n")

    def visualize_thermostability(self, features1, features2):
        # Create an interactive bar chart
        if self.matplotlib_canvas:
            self.matplotlib_canvas.get_tk_widget().destroy()  # Destroy the previous canvas

        fig, ax = plt.subplots(figsize=(8, 3))
        bars = ax.bar(["Aromatic", "Extremophilic"], [features1[0], features1[1]], color="#2196F3", alpha=0.7,
                      label='Enzyme 1')
        bars2 = ax.bar(["Aromatic", "Extremophilic"], [features2[0], features2[1]], color="#4CAF50", alpha=0.7,
                       label='Enzyme 2')

        # Add labels and legend
        ax.set_ylabel('Ratio')
        ax.legend()

        # Embed the matplotlib figure in Tkinter
        self.matplotlib_canvas = FigureCanvasTkAgg(fig, master=self.visualization_canvas)
        self.matplotlib_canvas.draw()
        self.matplotlib_canvas.get_tk_widget().pack(side='top', fill='both', expand=1)

    def visualize_sequence_alignment(self, sequence1, sequence2):
        alignments = pairwise2.align.globalxx(sequence1, sequence2, one_alignment_only=True, score_only=True)
        alignment_str = format_alignment(*pairwise2.align.globalxx(sequence1, sequence2, one_alignment_only=True)[0])
        self.result_label.insert("end", f"Sequence Alignment Score: {alignments}\n\nAlignment:\n{alignment_str}\n")

    def gc_content_calculation(self, sequence):
        # Add GC content calculation
        gc_content = gc_fraction(sequence)
        self.result_label.insert("end", f"GC Content: {gc_content:.2%}\n")

    def compare_thermostability(self):
        sequence1 = self.enzyme1_text.get("1.0", "end").strip()
        sequence2 = self.enzyme2_text.get("1.0", "end").strip()

        if sequence1 and sequence2:
            if self.is_valid_amino_acid_sequence(sequence1) and self.is_valid_amino_acid_sequence(sequence2):
                self.result_label.delete("1.0", "end")  # Clear previous results
                self.result_label.insert("end", "Calculating thermostability features... Please wait.\n")
                self.master.update()

                features1 = self.calculate_thermostability_features(sequence1)
                features2 = self.calculate_thermostability_features(sequence2)

                result = f"Enzyme 1 Features: Aromatic Ratio={features1[0]:.2f}, Extremophilic Ratio={features1[1]:.2f}\n"
                result += f"Enzyme 2 Features: Aromatic Ratio={features2[0]:.2f}, Extremophilic Ratio={features2[1]:.2f}\n"

                if features1[0] > features2[0] and features1[1] > features2[1]:
                    result += "Enzyme 1 is predicted to be more thermostable.\n"
                elif features1[0] < features2[0] and features1[1] < features2[1]:
                    result += "Enzyme 2 is predicted to be more thermostable.\n"
                else:
                    result += "The predicted thermostability is similar for both enzymes.\n"

                self.result_label.insert("end", result)

                # Visualize thermostability features
                self.visualize_thermostability(features1, features2)

                # Visualize sequence alignment
                self.visualize_sequence_alignment(sequence1, sequence2)

                # Additional analyses
                self.gc_content_calculation(sequence1)
                self.gc_content_calculation(sequence2)

            else:
                return  # Invalid sequence, error message already shown
        else:
            self.result_label.insert("end", "Error: Both enzyme sequences are required.\n")

    def is_valid_amino_acid_sequence(self, sequence):
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        invalid_characters = [aa for aa in sequence if aa.upper() not in valid_amino_acids]
        if invalid_characters:
            self.result_label.insert("end",
                                     f"Error: Invalid amino acid characters found: {', '.join(invalid_characters)}\n")
            return False
        return True

    def calculate_thermostability_features(self, sequence):
        aromatic_ratio = self.calculate_aromatic_ratio(sequence)
        extremophilic_ratio = self.calculate_extremophilic_ratio(sequence)
        return aromatic_ratio, extremophilic_ratio

    def calculate_aromatic_ratio(self, sequence):
        aromatic_aa = set("YWF")
        total_aa = len(sequence)
        aromatic_count = sum(1 for aa in sequence if aa.upper() in aromatic_aa)
        aromatic_ratio = aromatic_count / total_aa
        return aromatic_ratio

    def calculate_extremophilic_ratio(self, sequence):
        extremophilic_aa = set("DERK")
        total_aa = len(sequence)
        extremophilic_count = sum(1 for aa in sequence if aa.upper() in extremophilic_aa)
        extremophilic_ratio = extremophilic_count / total_aa
        return extremophilic_ratio

    def reset_input(self):
        # Clear input sequences and visualization
        self.enzyme1_text.delete("1.0", "end")
        self.enzyme2_text.delete("1.0", "end")
        if self.matplotlib_canvas:
            self.matplotlib_canvas.get_tk_widget().destroy()
        self.result_label.delete("1.0", "end")

    def clear_results(self):
        # Clear only the analysis results
        if self.matplotlib_canvas:
            self.matplotlib_canvas.get_tk_widget().destroy()
        self.result_label.delete("1.0", "end")


if __name__ == "__main__":
    root = Tk()
    app = LaccaseThermostabilityAnalyzer(root)
    root.mainloop()
