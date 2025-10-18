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

```bash
java GMPCompress <dataset_name>
