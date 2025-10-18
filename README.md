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

### 1️ Compression

To compress a protein dataset, run GMPCompress.java file.
Give dataset name in 
String datasetName = "HI"; -->Dataset name here


This will:

Load the dataset "HI"

Discover recurring amino acid motifs using a Genetic Algorithm

Store pattern–position mappings in a dictionary (𝒟)

Compress the updated database using the integrated compressor (AC2 backend)

Output files:

HI.co → Compressed data

HI_dictionary.xz → Encoded dictionary


### 2 Decompression
Decompress HI.co using AC2 decompressor
Decoding dictionary (RLE) (HI_dictionary.xz)
Rebuilding sequence structure
Original file restored: HI_restored


## License

This project is distributed under the MIT License — free for academic and research use.

---

## Acknowledgements

GMP builds upon the concepts of Minimum Description Length (MDL) pattern mining and integrates components inspired by:

-   AC and AC2 compressors (context modeling)
-   Genetic optimization for subsequence discovery
-   LZMA2 for efficient dictionary encoding

---

## Contact

For questions, collaboration, or citation updates:

**Muhammad Zohaib Nawaz**  
📧 mahsoodzohaib@gmail.com  
🌐 https://www.philippe-fournier-viger.com
