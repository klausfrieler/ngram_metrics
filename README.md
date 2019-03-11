# ngram_metrics
Evaluation scripts for automated melody extractions from the Weimar Jazz Database.

## Usage

Clone the repository and put your extracted note tracks in a directory of your choice (we recommend `data/test_set`). 
You need to have the `tidyverse` and `optparse`package installed to run the script. 

You can start the evaluation script on the commandline with 

``` console
Rscript run_eval.R -t [TEST_SET_DIRECTORY] -o [OUTPUT_DIR] -m [MAX_N] -f [FILE_FORMAT = en, de]
```
The script execute two kinds of evaluation. The first by matching notes and calculating recall, precision and F-scores, and the second by comparing pattern similarity for a random subset of solo pairs.

Results will be written to `OUTPUT_DIR` (default: 'output'), using n-grams up to length `MAX_N` (default: 10). Results are four CSV files, using the language convention as specificed in `FILE_FORMAT` (default is en), and a boxplot of F1 scores across n-gram lengths and a set fo thresholds (currently 30, 50 and 70 ms, see below). Be patient, grab a coffee or have a smoke, as the script is rather slow.

##Output 
###Retrieval Metrics
The script creates three CSV files based on a direct note-by-note matching procedure, which results are regarded as a classical retrieval problem. Note tracks are represented as `(onset, pitch)` pairs. For each note in the ground truth note track ('target')  all note events within a certain time window `T` of the  target onset are retrieved from the extracted note track ('query'). There can be zero, one or more  notes in this result set. Next, the pitches in the result set are compared  with the target pitch. All note events with the same pitch are  counted as  true positives (`TP`), note events with a different pitch are thus false positives (`FP`), and a empty result set is regarded as a false negative (`FN`). This procedure is repeated for each target event in each solo in the test set over a range of thresholds `T`, currently 30, 50, and 70 ms, and a range of n-gram length, as defined by the `MAX_N` parameter. The matching algorithm for n-grams of length `N` is similar to the one for note events (unigrams), i.e., n-grams in the target are matched with n-grams of the same length in the query, whereby the onset of an n-gram is defined here as the mean of onsets of all elements in a query n-gram.    

The first result file `ngram_stats_raw.csv` contains lists of the raw result. An example looks like this

``` console
threshold;direction;target_solo;target_n;target_pos;TP;FP;FN
0.03;target;BenWebster_ByeByeBlackbird_FINAL.csv;1;1;1;0;0
0.03;target;BenWebster_ByeByeBlackbird_FINAL.csv;1;2;1;0;0
0.03;target;BenWebster_ByeByeBlackbird_FINAL.csv;1;3;1;0;0
``` 
The first column contains the threshold value (in seconds). The second column is the direction of match (currently always `target`, so it can be safely ignored; this parameter is reserved for future use, where matching could be done with reversed roles of target and query note tracks). The third column is the name of the target solo. The fourth is the length and the fifth the position of the n-gram in the target solo. The next three columns contain the evaluation values, i.e.,  the number of of true positives, false positives and false negative in the results set for an n-gram. There can be more than one true positive, e.g., when the event got split up into several events by the note tracker. Likewise, there can be also more than one false positive, and even a mixture of true and false positives, e.g., when the note tracker resolved a tone with strong vibrato into a sequence of short oscillating pitches. If the result set is empty, then `TP` and `FP` are 0 and `FN` is 1.

The second result file, `ngram_stat_solo.csv`, is an aggregated version over solos of the first one. For example:

``` console
threshold;direction;target_n;target_solo;prec;rec;F1
0.03;target;1;BenWebster_ByeByeBlackbird_FINAL.csv;0.884816753926702;0.800947867298578;0.840796019900497
0.03;target;1;BenWebster_DidYouCallHerToday_FINAL.csv;0.86144578313253;0.742857142857143;0.797768479776848
```

The first four columns have the same meaning as in the first result file. The next three columns contain precision, recall, and F1 values, as defined by the usual formulas. However, since the result set of a target n-grams can contain more than one element, there are certain degrees of freedom in calculating precision, recall and F1 scores. We adopt here a rather strict policy, in defining a 'really' true positive when there is exactly one matching n-gram, i.e., only if `TP = 1` and `FP = FN = 0`. A 'true' false positive is then defined if there is at least one false positive, i.e., if `FP > 0`. False negatives need no modification, since `FN = 1` implies `TP = FP = 0` by construction. Precision, recall and F1 scores are thus not calculated  on the raw `TP`, `FP`, and `FN` values, but based on these criteria. However, experiments show that the differences in result are actually only marginal, as these cases are basically only encountered for unigrams.

Finally, the third result file `ngram_stat_sum.csv` is an aggregation of the second result file for thresholds and n-gram lengths. For example:

``` console
threshold;direction;target_n;prec_mean;rec_mean;F1_mean;prec_median;rec_median;F1_median;prec_sd;rec_sd;F1_sd;prec_max;rec_max;F1_max;prec_min;rec_min;F1_min
0.03;target;1;0.917838050364794;0.81273770010961;0.860197903718169;0.923076923076923;0.840579710144927;0.877076411960133;0.0449646778826509;0.0908927493458668;0.0677448224623553;0.983957219251337;0.978070175438597;0.975929978118162;0.784313725490196;0.423326133909287;0.552112676056338
```
Here, the first three columns have again the same meaning as in the two other files. The remaining columns contain statistical descriptor for the precision, recall and F1 value distribution (mean, median, standard deviation (sd), minimum (min), and maximum (max)).

The fourth result file is a boxplot of F1 scores as the main metric, with n-gram lengths on the x-axis and three boxes for the different thresholds each.

###Pattern Similarity Metrics

An alternative n-gram metric is based on the idea of measuring the capability of the automatically extracted note track to reproduce some important pattern metrics, specifically, **pattern commonalities** (or n-gram similarities). The pattern commonality of two subsets of solos (note tracks) is defined as the total statistical variation of n-gram distributions over a range of n-grams. Total statistical variation is simply the L1 norm of two probability distributions. Fixing an n-gram length `N`, and two set of note tracks `S1` and `S2` to be compared, we calculate the n-gram probability distribution over the **common** set of n-grams in `S1` and `S2` for each note track. Averaging the absolute difference in (estimated) probabilities between the two sets is a measure of pattern commonality. 

The metrics here is then based on comparing pattern commonality (n-gram similarity) values calculated once for the ground truth note tracks and once for the extracted note tracks. This could be done for all possible `M (M-1)/2` pairs of solos in the test set, but this would be too time-consuming. Currently, it is done on 10 folds of 10 randomly selected solos (550 pairs in total in each fold) and a range of n-gram lengths as defined by the `MAX_N` parameter. 

There are two result files produced. The first, `pat_sim_eval_raw.csv`, looks like this

``` console
N;id1;id2;common_wjd;common_db;d_sim;d_sim_rel;batch
1;WyntonMarsalis_CherokeeII_FINAL.csv;SidneyBechet_ReallyTheBlues_FINAL.csv;0.669249263689367;0.703283848538946;-0.0340345848495788;-0.0508548708174238;1
1;WyntonMarsalis_U.M.M.G._FINAL.csv;SidneyBechet_ReallyTheBlues_FINAL.csv;0.547968141913295;0.580031874967868;-0.0320637330545726;-0.0585138634932649;1
```

where `N` is the n-gram length, `id1` and `id2` are the ids of the solos of the comparison pair, `common_wjd` is the pattern commonality of this solos pair calculated based on the ground truth, and `common_db` the one based on the extracted note tracks. `d_sim` is the difference between these values, and `d_sim_rel` is this difference divided by  `common_wjd`. Finally, `batch` is the number of the fold.

The second result file `pat_sim_eval_sum.csv` contains aggregated values as well as correlations. It looks like this:

``` console
N;cor;mean_common_wjd;mean_d;rel_mean_d;sd_d;min_d;max_d
1;0.972420178410897;0.558787961670915;-0.0141644868698324;-0.025348589879204;0.0352382932487599;-0.130891169731448;0.0597043356490319
```

THe first column is the n-gram length `N`, second contains (Pearson) correlations between pattern commonality values of the ground truth and the extracted note tracks. The third column is `mean_common_wjd`, the mean pattern commonality for this `N` in the ground truth and the next column `mean_d` is the mean difference between ground truth and extracted note set pattern commonality. The next column is the relative difference, hence, `mean_d` divided by `mean_common_wjd`. The last three colums contains standard deviation (`sd_d`), mininum (`min_d`) and maximum (`max_d`) of the difference.

