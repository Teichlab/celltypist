# Celltypist

**_Automated cell type annotation for scRNA-seq datasets_**

## Install

```console
pip install celltypist-dev
```

## Usage

## Python package

### Sample data and default model

```python
import celltypist

sample_input = celltypist.samples.get_sample_csv()
result = celltypist.annotate(sample_input)
result.predicted_labels_as_df().to_csv("labels.csv")
```

### Using your own data

```python
import celltypist

input_data = "/path/to/cell_by_gene_matrix.csv"
result = celltypist.annotate(input_data)
result.write_excel("annotation_result.xlsx")
```

### Using included models

Included models are downloaded automatically inside a inscript when the defualt model tries to be loaded or by using `celltypist --update-models` from the command line.

```python
import celltypist

input_data = "/path/to/cell_by_gene_matrix.csv"
result = celltypist.annotate(input_data, model="default")
result.write_excel("annotation_result.xlsx")
```

#### List models included in with the package

```python
import celltypist

available_models = celltypist.models.get_all_models()
available_models
```

### Use custom models

```python
import celltypist

indata = "/path/to/cell_by_gene_matrix.csv"
custom_model = "/path/to/custom_model.pkl"

result = celltypist.annotate(input_data, model=custom_model)
result.write_excel("annotation_result.xlsx")
```

### Full example

```python
import celltypist

indata = "/path/to/cell_by_gene_matrix.csv"
custom_model = "/path/to/custom_model.pkl"

result = celltypist.annotate(input_data, model=custom_model)
result.write_excel("annotation_result.xlsx")

print(result.summary_as_df())
result.predicted_labels_as_df().to_csv("labels.csv")
result.probability_matrix_as_df().to_csv("prob_matrix.csv")
```

## Command Line

### Basic usage

```bash
celltypist --indata=/path/to/dataset.csv --model default
```

### Advance usage

```bash
celltypist
    --indata /path/to/dataset.csv \ # input dataset
    --model /path/to/model.pkl    \ # path to model
    --outpref Dataset1_           \ # add a prefix to the output files
    --outdir /path/to/output      \ # set an output directory for the files
```
