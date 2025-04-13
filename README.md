<p align="center">
  <img src="logo.png" alt="filter_kraken2_report logo" width="200"/>
</p>

# filter_kraken2_report

A Python tool for filtering [Kraken2](https://ccb.jhu.edu/software/kraken2/) reports by removing specified taxa and all their descendants. Useful for excluding known contaminants, host sequences, or unwanted clades from downstream analysis.

## Features

- Removes a specified taxon and all its children from a Kraken2 report
- Automatically recalculates percentages based on the remaining dataset
- Outputs a filtered report and a list of removed taxa

## Requirements

- Python 3.7+
- `pandas`
- `networkx`

You can install the required packages with:

```bash
pip install pandas networkx
```

or with Conda

## Usage

```bash
python filter_kraken2_report.py \
  -r input_report.txt \
  -f 9606,12333 \
  -o filtered_report.txt \
  -p removed_report.txt
```

### Arguments

| Argument           | Description                                                                 |
|--------------------|-----------------------------------------------------------------------------|
| `-r, --report`      | Input Kraken2 report file                                                   |
| `-f, --filter_ids`  | Comma-separated list of tax IDs to remove, e.g. `9606,12333`                |
| `-o, --out_fn`      | Output file for the filtered Kraken2 report                                 |
| `-p, --removed_out_fn` | Output file listing all removed entries from the original report         |

## Notes

- The script uses indentation in the Kraken2 report to reconstruct taxonomic hierarchy.
- Taxa not found in the hierarchy will be skipped with a warning.
## License

Licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.
