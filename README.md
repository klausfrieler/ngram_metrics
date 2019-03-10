# ngram_metrics
Evaluation scripts for automated melody extractions from the Weimar Jazz Database.

## Usage

Clone the repository and put your extracted note tracks in a directory of your choice (we recommend `data/test_set`). 
You need to have the `tidyverse` and `optparse`package installed to run the script. 

You can start the evaluation script on the commandline with 

``` console
Rscript run_eval.R -t [TEST_SET_DIRECTORY] -o [OUTPUT_DIR] -m [MAX_N] -f [FILE_FORMAT = en, de]
```
The script execute two kinds of evaluation. The first by matching notes and calculating recall, precision and F-scores, and the second by comparing pattern similarity for a random subset of 550 solo paris.

Results will be written to `OUTPUT_DIR` (default: 'output'), using a ngrams up to length `MAX_N` (default: 10). Results are four  CSV files, using the language convention as specificed in `FILE_FORMAT` (default is en), and a boxplot of F1 scores across n-gram lengths and thresholds (currently 30, 50 and 70 ms). Be patient, grab a coffee or have a smoke, as the script is rather slow.

##Output 
###Retrieval Metrics
The script creates three CSV files from a direct note-by-note matching procedure, which results are regarded as a classical retrieval problem. Note tracks are represented as `(onset, pitch)$`  pairs. For each note in the ground truth note track ('target') the notes in the extracted note track ('query') within a time window `T` of the are target onset are retrieved. There can be zero, one multiples notes in this result set. Next, the pitches are compared between the target and the result set. A note event with matching pitches is counted as a true positives (`TP`), a note event with a different pitch is a  false positives (`FP`), and a empty result set is regarded as a false negative (`FN`). This procedure is reperated for each target event, for each solo in the test set and for a range of thresholds, currently 30, 50, and 70 ms, and also for a range of n-grams, default is 1 to 10 but is a parameter for the script. The matching algorithm for n-grams of length `N` is similar to the one for note evens (unigrams), matching n-grams  in the target with n-grams in the query of the same length, whereby the onset of an n-gram is taken to be the mean of onsets of all elements in the n-gram.    

The first result file `ngram_stats_raw.csv` contains lists of the raw result. An example looks like this

``` console
threshold;direction;target_solo;target_n;target_pos;TP;FP;FN
0.03;target;BenWebster_ByeByeBlackbird_FINAL.csv;1;1;1;0;0
0.03;target;BenWebster_ByeByeBlackbird_FINAL.csv;1;2;1;0;0
0.03;target;BenWebster_ByeByeBlackbird_FINAL.csv;1;3;1;0;0
``` 
The first column contains the threshold. The second the direction of match (currently only `target`, so it can be savely ignored; sometimes matching the query to the target could be interest, so this column is reserved for future use). The third column is the name of the target solo. The fourth is the length of n-grams and the fifth the position of the n-gram as an enumeration position in the target solo. The next three columns contain the analysis, the number of of true positives, false positives and false negative in the results set for this n-gram. There can more than one treu positive, e.g., when the event got split up in more event by the note tracker. Likewise, there can be also more than one false positive, and even a mixture of true and false positive. If the result was empty, then TP and FP are zero and FN is one.

The second result file, `ngram_stat_solo.csv` is an aggregated version overf solos of the first one. For example:

``` console
threshold;direction;target_n;target_solo;prec;rec;F1
0.03;target;1;BenWebster_ByeByeBlackbird_FINAL.csv;0.884816753926702;0.800947867298578;0.840796019900497
0.03;target;1;BenWebster_DidYouCallHerToday_FINAL.csv;0.86144578313253;0.742857142857143;0.797768479776848
```

The first four columns are the same as in the first result file. The next three contain precision, recall, and F1 values, as defined by the usual formulas. However, since the result set of a target n-gram can contain more than one element, there is certain freedom of calculating precision, recall and F1 scores, which however only make slight differences in the values as this basically only happens for unigrams. A really true positive is encountered when there is only one matching n-gram, i.e., `TP = 1` and `FP = FN = 0`, whereas a true a a true fals positive is defined if there is at least one false positive, hence, when `FP > 0`. False negatives need no modifcation, obviously. So, precision, recall and F1 score are calculated not based on the raw TP, FP, and FN values but based on this criteria and stored in this file. As mentioned earlier, the differences in result are marginal using either approach.

Finally, the third result file `ngram_stat_sum.csv` is an aggregation of the second result file for thresholds and n-gram length. For example:

``` console
threshold;direction;target_n;prec_mean;rec_mean;F1_mean;prec_median;rec_median;F1_median;prec_sd;rec_sd;F1_sd;prec_max;rec_max;F1_max;prec_min;rec_min;F1_min
0.03;target;1;0.917838050364794;0.81273770010961;0.860197903718169;0.923076923076923;0.840579710144927;0.877076411960133;0.0449646778826509;0.0908927493458668;0.0677448224623553;0.983957219251337;0.978070175438597;0.975929978118162;0.784313725490196;0.423326133909287;0.552112676056338
```
Here, the first three columns have again the same meaning as in the two other files. The remaining columns contain statistical descriptor for the precision, recall and F1 vale distribution (mean, median, standard deviation (sd), minimum (min), and maximum (max)).

The fourth result file is a boxplot of F1 scores, as the main metric, with n-gram length on the x-axis and three boxes each fot the different thresholds.

###Pattern Similarity Metrics

An alternative n-gram metric is based on the measuring the capability of the automatically extracted note track to reproduce some important pattern metrics, specifically, n-gram commonalities. Pattern commonality of two subsets of solos (target note tracks) is defined as the total statistical variation of n-gram distributions over a range of n-grams. Total statistical variation is the L1 norm of two probability distributions. Fixing an n-gram length `N`, and two set of note tracks `S1` and `S2` to be compared, we calculate the probability distribution over the **common** set of n-grams in `S1` and `S2`. Averaging the absolute difference in (estimated) probability between the two sets is a measure of n-gram commonality. The metrics here is now based on comparing these pattern commonality values (or n-gram similariy values) calculated once on the ground truth note tracks and once on the extracted note tracks. Ideally, this done for all possible `M (M-1)/2` pairs of solos in the test set, but this would be too time-consuming, so currently this done on 10 folds of 10 randomly selected solos and the range of n-gram lengths as defined by the `MAX_N` parameter. 

There are two result files. The first, `pat_sim_eval_raw.csv` looks like this

``` console
N;id1;id2;common_wjd;common_db;d_sim;d_sim_rel;batch
1;WyntonMarsalis_CherokeeII_FINAL.csv;SidneyBechet_ReallyTheBlues_FINAL.csv;0.669249263689367;0.703283848538946;-0.0340345848495788;-0.0508548708174238;1
1;WyntonMarsalis_U.M.M.G._FINAL.csv;SidneyBechet_ReallyTheBlues_FINAL.csv;0.547968141913295;0.580031874967868;-0.0320637330545726;-0.0585138634932649;1
```

where `N` is the n-gram length, `id1` and `id2` are the ids of the solos of the comparison pair, `common_wjd` is the pattern commonality of this solos pair calculated based on the ground truth, and `common_db` the one based on the extracted note tracks. `d_sim` is the difference between these values, and `d_sim_rel` is this difference divided by  `common_wjd`. Finally, `batch` is the number of the fold.





