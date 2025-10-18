# 🧬 GMP: Genetic Algorithm–Based MDL Protein Compressor

This repository contains the implementation of **GMP (Genetic Algorithm–based MDL Protein Compressor)** — a hybrid framework designed to compress protein sequences by discovering biologically meaningful recurring amino acid motifs (**kAAMs**) and encoding them efficiently.

---

## 📁 Project Structure

| File | Description |
|------|--------------|
| `GMPCompress.java` | Main compression module — discovers frequent amino acid subsequences (kAAMs), builds the dictionary (`𝒟`), and compresses the input protein sequence. |
| `GMPDecompress.java` | Decompression module — reconstructs the original protein sequence from the compressed database and dictionary. |
| `README.md` | Project documentation (this file). |

---

## ⚙️ Requirements

- **Java 8** or higher  
- Sufficient memory for large protein datasets  
- Protein sequence files in plain text (`.txt`) or FASTA format (`.fasta`)

---

## 🚀 Usage

### 1️⃣ Compression

To compress a protein dataset, run:

java GMPCompress <dataset_name>

This will:

Load the dataset HI.fasta (or HI.txt)

Discover recurring amino acid motifs using a Genetic Algorithm

Store pattern–position mappings in a dictionary (𝒟)

Compress the updated database using the integrated compressor (AC2 backend + RLE or LZMA2)

Output files:

HI_compressed.bin → Compressed data

HI_dictionary.bin → Encoded dictionary
