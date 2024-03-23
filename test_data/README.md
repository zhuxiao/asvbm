# SV_STAT experiments

## experiment description

In this experiment, we conducted an individualized benchmarking of the identification results for chromosome 1 (chr1) from the HG002 CCS (Consensus Circular Sequence) data by using seven structural variants (SV) detection methods: cuteSV(v2.0.3), pbsv (v2.9.0), Sniffles2(v2.0.2), asvclr(v1.4.0), SVDSS(v1.0.5), DeBreak(v1.0.2), SVIM(v2.0.0). Here we will only focus on the benchmarking results of insertion and deletion variants. Subsequently, the identification results obtained from these seven different detection methods were comprehensively benchmarked against their benchmark sets through a multi-sample benchmarking approach.

## data description and preparation


In this experiment benchmarking, we used the chr1 of hg37 as reference, and human reference genome hs37d5 was downloaded as follows:
``` 
$ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz. 
```
We also used the chr1 of the benchmark set of HG002 and the detailed information of the benchmark set (HG002_SVs_Tier1_v0.6.vcf) was downloaded as follows:
```
$ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6
``` 

```
## Download source programs and test data
$ wget -c https://github.com/zhuxiao/sv_stat/releases/download/1.0.0/sv_stat_1.0.0.tar.xz
$ tar -xf sv_stat_1.0.0.tar.xz
$ cd sv_stat_1.0.0/
$ ./autogen.sh
```
The name of the chr1 file (reference) was renamed to test_chr1.fa. We uploaded the sequence of chr1 into three files, which can be downloaded and merged by following steps:
```
## Merge chr1 references and extract it
$ cd sv_stat/test_data/
$ cat test_chr1.fa.tar.xz_* > test_chr1.fa.tar.xz
$ tar -xJf test_chr1.fa.tar.xz
```
Extract the detection results from compressed packages:
```
## The identification results of chr1 by cuteSV were extracted
$ tar -xJf test_cutesv_chr1.vcf.tar.xz  
## The identification results of chr1 by pbsv were extracted
$ tar -xJf test_pbsv_chr1.vcf.tar.xz
## The identification results of chr1 by Sniffles were extracted
$ tar -xJf test_sniffles_chr1.vcf.tar.xz
## The identification results of chr1 by asvclr were extracted
$ tar -xJf test_asvclr_chr1.vcf.tar.xz
## The identification results of chr1 by SVDSS were extracted
$ tar -xJf test_SVDSS_chr1.vcf.tar.xz
## The identification results of chr1 by DeBreak were extracted
$ cat test_DeBreak_chr1_tar_a* > test_DeBreak_chr1.tar.xz
$ tar -xJf test_DeBreak_chr1.tar.xz
## The identification results of chr1 by SVIM were extracted
$ tar -xJf test_SVIM_chr1.tar.xz

## Benchmarking the VCF files mentioned in the test data using SV_STAT, run the following command:
$ sv_stat -m 50000 -T "cuteSV;pbsv;Sniffles2;asvclr;SVDSS;DeBreak;SVIM" -C "1" test_data/test_cutesv_chr1.vcf test_data/test_pbsv_chr1.vcf test_data/test_sniffles_chr1.vcf test_data/test_asvclr_chr1.vcf test_SVDSS_chr1.vcf test_DeBreak_chr1.vcf test_SVIM_chr1.vcf test_benchmark_chr1.vcf test_data/test_chr1.fa
```
SV_STAT benchmarked the provided identification set based on the input, generating corresponding file information. If needed, more detailed benchmarking information can be obtained by running test data independently.

## result description
The benchmarking results were as follows:
<table>
  <thead>
    <tr>
    	  <th></th>
    	  <th>SVs</th>
      <th>TP_benchmark</th>
      <th>TP_user</th>
      <th>FP</th>
      <th>FN</th>
      <th>Recall</th>
      <th>precision</th>
      <th>F1 score</th>
      <th>Seqcons</th>
    </tr>
  </thead>
    <tbody>
    <tr>
    <th>asvclr</th>
      <th>3893</th>
      <th>3374</th>
      <th>3276</th>
      <th>617</th>
      <th>2023</th>
      <th>0.6252</th>
      <th>0.8415</th>
      <th>0.7174</th>
      <th>0.9183</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>SVDSS</th>
      <th>3318</th>
      <th>2506</th>
      <th>2759</th>
      <th>559</th>
      <th>2891</th>
      <th>0.4643</th>
      <th>0.8315</th>
      <th>0.5959</th>
      <th>0.9360</th>
    </tr>
  </tbody>
    <tbody>
    <tr>
    <th>DeBreak</th>
      <th>3429</th>
      <th>3033</th>
      <th>2854</th>
      <th>467</th>
      <th>2364</th>
      <th>0.5620</th>
      <th>0.8594</th>
      <th>0.6796</th>
      <th>0.9318</th>
    </tr>
  </tbody>
     <tbody>
    <tr>
    <th>Sniffles2</th>
      <th>3972</th>
      <th>3318</th>
      <th>3189</th>
      <th>676</th>
      <th>2079</th>
      <th>0.6148</th>
      <th>0.8251</th>
      <th>0.7046</th>
      <th>0.9220</th>
    </tr>
  </tbody> 
   <tbody>
    <tr>
    <th>pbsv</th>
    	<th>3859</th>
      <th>3241</th>
      <th>3139</th>
      <th>635</th>
      <th>2156</th>
      <th>0.6005</th>
      <th>0.8317</th>
      <th>0.6975</th>
      <th>0.9662</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>cuteSV</th>
      <th>3384</th>
      <th>2914</th>
      <th>2734</th>
      <th>488</th>
      <th>2483</th>
      <th>0.5399</th>
      <th>0.8485</th>
      <th>0.6599</th>
      <th>0.9184</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>SVIM</th>
      <th>11490</th>
      <th>3504</th>
      <th>3452</th>
      <th>2384</th>
      <th>1893</th>
      <th>0.6493</th>
      <th>0.5915</th>
      <th>0.6190</th>
      <th>0.9549</th>
    </tr>
  </tbody>
</table>

In addition, the center distance and the region size ratio were calculated for the variation overlapping between the user set and the benchmark set. The closer the center distance is to 0, the smaller the deviation and the more accurate the identification result. and the closer the ratio is to 1, the smaller the deviation, which means more accurate identification result.

The statistical results of center distance were shown as follows:
<table>
  <thead>
    <tr>
    	  <th></th>
    	  <th>-200--151</th>
      <th>-150--101</th>
      <th>-100--51</th>
      <th>-50--1</th>
      <th>0-50</th>
      <th>51-100</th>
      <th>101-150</th>
      <th>151-200</th>
    </tr>
  </thead>
  <tbody>
    <tr>
    <th>asvclr</th>
      <th>12</th>
      <th>34</th>
      <th>46</th>
      <th>418</th>
      <th>2250</th>
      <th>279</th>
      <th>168</th>
      <th>100</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>SVDSS</th>
      <th>22</th>
      <th>27</th>
      <th>82</th>
      <th>871</th>
      <th>1619</th>
      <th>159</th>
      <th>104</th>
      <th>72</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>DeBreak</th>
      <th>15</th>
      <th>28</th>
      <th>36</th>
      <th>157</th>
      <th>1881</th>
      <th>379</th>
      <th>175</th>
      <th>149</th>
    </tr>
  </tbody>
     <tbody>
    <tr>
    <th>Sniffles2</th>
      <th>28</th>
      <th>43</th>
      <th>56</th>
      <th>250</th>
      <th>2439</th>
      <th>308</th>
      <th>142</th>
      <th>96</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>pbsv</th>
    	  <th>23</th>
      <th>30</th>
      <th>55</th>
      <th>334</th>
      <th>2378</th>
      <th>175</th>
      <th>106</th>
      <th>83</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>cuteSV</th>
      <th>18</th>
      <th>31</th>
      <th>50</th>
      <th>361</th>
      <th>2009</th>
      <th>238</th>
      <th>117</th>
      <th>77</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>SVIM</th>
      <th>43</th>
      <th>59</th>
      <th>87</th>
      <th>454</th>
      <th>2594</th>
      <th>278</th>
      <th>178</th>
      <th>131</th>
    </tr>
  </tbody>
</table>

The statistical results of the region size ratio were shown as follows:

<table>
  <thead>
    <tr>
    	  <th></th>
    	  <th>0.0-0.5</th>
      <th>0.5-0.7</th>
      <th>0.7-1.2</th>
      <th>1.2-2.0</th>
      <th>2.0-5.0</th>
      <th>5.0-10.0</th>
      <th>10.0-50.0</th>
      <th>50.0-100.0</th>
      <th>>100.0</th>
    </tr>
  </thead>
     <tbody>
    <tr>
    <th>asvclr</th>
      <th>280</th>
      <th>21</th>
      <th>3086</th>
      <th>41</th>
      <th>35</th>
      <th>10</th>
      <th>54</th>
      <th>37</th>
      <th>239</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>SVDSS</th>
      <th>347</th>
      <th>25</th>
      <th>2859</th>
      <th>31</th>
      <th>32</th>
      <th>3</th>
      <th>27</th>
      <th>18</th>
      <th>9</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>DeBreak</th>
      <th>125</th>
      <th>17</th>
      <th>1543</th>
      <th>39</th>
      <th>54</th>
      <th>16</th>
      <th>736</th>
      <th>182</th>
      <th>589</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>Sniffles2</th>
      <th>342</th>
      <th>23</th>
      <th>3279</th>
      <th>31</th>
      <th>37</th>
      <th>6</th>
      <th>40</th>
      <th>18</th>
      <th>118</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>pbsv</th>
    	<th>310</th>
      <th>21</th>
      <th>2696</th>
      <th>40</th>
      <th>53</th>
      <th>6</th>
      <th>315</th>
      <th>84</th>
      <th>237</th>
    </tr>
  </tbody>
    <tbody>
    <tr>
    <th>cuteSV</th>
      <th>276</th>
      <th>16</th>
      <th>2807</th>
      <th>33</th>
      <th>26</th>
      <th>7</th>
      <th>38</th>
      <th>19</th>
      <th>152</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>SVIM</th>
      <th>571</th>
      <th>25</th>
      <th>3506</th>
      <th>40</th>
      <th>54</th>
      <th>13</th>
      <th>82</th>
      <th>58</th>
      <th>362</th>
    </tr>
  </tbody>
</table>

Additionally, basic metrics for different structural variant (SV) size ranges were computed, primarily categorized into the following seven intervals. The results of asvclr were showed as follows:

<table>
  <thead>
    <tr>
    	  <th>region</th>
      <th>TP_benchmark</th>
      <th>TP_user</th>
      <th>FP</th>
      <th>FN</th>
      <th>Recall</th>
      <th>precision</th>
      <th>F1 score</th>
      <th>Seqcons</th>
    </tr>
  </thead>
  <tbody>
    <tr>
    	  <th>1-100bp</th>
      <th>2836</th>
      <th>2754</th>
      <th>532</th>
      <th>1885</th>
      <th>0.6007</th>
      <th>0.8381</th>
      <th>0.6998</th>
      <th>0.9102</th>
    </tr>
    <tr>
    	  <th>101-250bp</th>
      <th>119</th>
      <th>111</th>
      <th>73</th>
      <th>124</th>
      <th>0.4897</th>
      <th>0.6033</th>
      <th>0.5406</th>
      <th>0.9927</th>
    </tr>
    <tr>
    	  <th>251-500bp</th>
      <th>126</th>
      <th>124</th>
      <th>136</th>
      <th>90</th>
      <th>0.5833</th>
      <th>0.4769</th>
      <th>0.5248</th>
      <th>0.9922</th>
    </tr>
    <tr>
    	  <th>501-1000bp</th>
      <th>31</th>
      <th>31</th>
      <th>41</th>
      <th>45</th>
      <th>0.4079</th>
      <th>0.4306</th>
      <th>0.4189</th>
      <th>0.9988</th>
    </tr>
    <tr>
    	  <th>1001-5000bp</th>
      <th>41</th>
      <th>43</th>
      <th>22</th>
      <th>51</th>
      <th>0.4457</th>
      <th>0.6615</th>
      <th>0.5326</th>
      <th>0.9997</th>
    </tr>
    <tr>
    	  <th>5001-1000bp</th>
      <th>10</th>
      <th>10</th>
      <th>9</th>
      <th>16</th>
      <th>0.3846</th>
      <th>0.5263</th>
      <th>0.4444</th>
      <th>1.0000</th>
    </tr>
    <tr>
    	  <th>>=10001bp</th>
      <th>4</th>
      <th>4</th>
      <th>3</th>
      <th>19</th>
      <th>0.1739</th>
      <th>0.5714</th>
      <th>0.2667</th>
      <th>1.0000</th>
    </tr>
  </tbody>
</table>

The comprehensive benchmarking results of the four detection methods are shown in the figure as follows:

<div style="text-align:center;">
<img src="img/Performance comparison between different tools.png" alt= "Performance comparison between different tools" width= "400" heiht="300" style="margin-right:25px"> 
<img src="img/Benchmark results between different tools.png" alt= "Benchmark results between different tools" width= "400" heiht="300" style="margin-left:25px;" >
</div>
