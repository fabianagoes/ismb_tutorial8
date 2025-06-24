## About IP8

This tutorial focuses on learning various ways to represent sequences for machine learning applications. Practical 1 focuses on traditional representations such as numeric mappings and sequence composition features while Practical 2 focuses on embeddings techniques for large language models (LLMs).

## Getting ready for IP8
### What would you need?

- gmail account (to access google colaboratory). Google colaboratory lets you write and execute Python code in your browser, with access to GPUs and TPUs.
- github account  and [git bash](https://git-scm.com/downloads) (to save the codebase locally).

### Running on Google Colab

To use Google Colab simply click on the links below and run the cells. Each cell has related instructions. To run the tutorials on Google Colab:

1. Click on the notebook links provided below to open them directly in Colab.
2. Follow the step-by-step instructions included in each notebook cell.
3. Run each cell by clicking the ▶️ button on the left side of the code block.

> **Tip:** Google Colab lets you write and execute Python code for free in your browser, with access to GPUs and TPUs.

### Cloning the git repository locally (not compulsory)

If you wish to keep a copy of the repo locally

1. Open the terminal
2. Clone the repository:
   ```bash
   git clone https://github.com/fabianagoes/ismb_tutorial8.git
   ```

## Practical 1

Practical 1 focuses on learning representations to encode nucleic acid sequences and using these encodings/representations to train and evaluate simple machine learning models for classification tasks.

[Open Practical 1 in Colab](https://colab.research.google.com/github/fabianagoes/ismb_tutorial8/blob/main/tutorial_practical1_colab.ipynb)

## Practical 2

Practical 2 focuses on learning vector representations to encode nucleic acid sequences using word embedding techniques and large language models (LLMs). The learned embeddings are then used to train and evaluate some popular LLMs for sequence-based prediction tasks.

[Open Practical 2.1 in Colab](https://colab.research.google.com/github/fabianagoes/ismb_tutorial8/blob/main/tutorial_practical2_1_colab.ipynb)

[Open Practical 2.2 in Colab](https://colab.research.google.com/github/fabianagoes/ismb_tutorial8/blob/main/tutorial_practical2_2_colab.ipynb)

## Resources

### Feature Engineering tools
- iFeatureOmega: https://github.com/Superzchen/iFeature
- iLean: https://github.com/Superzchen/iLearn
- MathFeature: https://github.com/Bonidia/MathFeature

### LLMs
- DNABERT2: https://github.com/MAGICS-LAB/DNABERT_2
- Nucleotide Tranformer: https://github.com/instadeepai/nucleotide-transformer
- RNA-FM: https://github.com/ml4bio/RNA-FM

### Datasets

The datasets used in this tutorial are based on the Genome Understanding Evaluation (GUE) benchmark proposed by DNABERT2 authors. GUE is a broad and well-structured benchmark designed to evaluate models across multiple genome analysis tasks. It includes a diverse collection of datasets covering various tasks and species, making it a valuable resource for assessing model performance in genomic sequence understanding. More details available at: https://arxiv.org/abs/2306.15006.
