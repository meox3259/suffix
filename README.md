***

## Usage

To run the program, use the following command-line syntax:

```bash
./your_program_name [options] -i
```

### Options

The following table lists all available command-line arguments:

| Short | Long Flag | Argument Type | Description |
| :--- | :--- | :--- | :--- |
| **`-i`** | — | `filename` | **Required**. Path to the input file. |
| **`-k`** | `--kmer` | `integer` | Sets the k-mer size used for processing. |
| **`-e`** | `--error` | `float` | Sets the expected error rate. |
| **`-n`** | `--nth` | `integer` | Sets the N threshold value. |
| **`-t`** | `--num_threads` | `integer` | Number of threads to use for parallel processing. |
| **`-v`** | `--verbose` | — | Enable verbose output mode (prints detailed logs). |
| **`-h`** | `--help` | — | Display the help message and exit. |

### Examples

**Basic usage:**
Run with default settings, specifying only input and output files:
```bash
./your_tool -i input.data
```

**Advanced usage:**
Run with a k-mer size of 31, using 8 threads, and a specific error rate:
```bash
./your_tool -i reads.fq result.txt --kmer 31 --error 0.05 --num_threads 8
```

**Verbose mode:**
Run with verbose logging enabled to debug or monitor progress:
```bash
./your_tool -i input.data -v
```