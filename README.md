# ngram_metrics
Evaluation scripts for automated melody extractions from the Weimar Jazz Database.

## Usage

Clone the repository and put your extracted note tracks in a directory of your choice (we recommend `data/test_set`. You can start the evaluation script on the commandline with 

``` console
Rscript ngram_metrics.R -t [TEST_SET_DIRECTORY] -o [OUTPUT_DIR] -m [MAX_N] -f [FILE_FORMAT = en, de]
```
Results will be written to `OUTPUT_DIR` (default: 'output'), using a ngrams up to length `MAX_N` (default: 10). Results are three CSV files, using the language convention as specificed in `FILE_FORMAT` (default is en), and a boxplot of F1 scores across n-gram lengths and thresholds (currently 30, 50 and 70 ms). Be patient, grab a coffee or have a smoke, as the script is rather slow.



