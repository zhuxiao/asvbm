# SV_STAT experiments

## experiment description

In this experiment, we conducted an individualized evaluation of the identification results for chromosome 1 (chr1) from the HG002 CCS (Consensus Circular Sequence) data by using four structural variant (SV) detection methods: cuteSV(v2.0.3), pbsv (v2.9.0), Sniffles(v2.0.2), asvclr(vxxx).  Subsequently, the identification results obtained from these four different detection methods were comprehensively evaluated against their benchmark sets through a multi-sample evaluation approach.

## data description and preparation


In this experiment evaluation, we used the chr1 of hg37 as reference, and human reference genome hs37d5 was downloaded as follows:
``` 
$ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz. 
```
We also used the chr1 of the benchmark set of HG002, and the deteiled information of the benchmark set (HG002_SVs_Tier1_v0.6.vcf) was downloaded as follows:
```
$ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6
``` 

```
## Download source programs and test data
$ wget -c https://github.com/zhuxiao/sv_stat/releases/download/0.8.0/sv_stat_0.8.0.tar.xz
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

## Evaluating the VCF files mentioned in the test data using SV_STAT, run the following command:
$ sv_stat -m 50000 -T "cuteSV;pbsv;Sniffles;asvclr" test_data/test_cutesv_chr1.vcf test_data/test_pbsv_chr1.vcf test_data/test_sniffles_chr1.vcf test_data/test_asvclr_chr1.vcf test_data/test_chr1.fa
```
SV_STAT evaluated the provided identification set based on the input, generating corresponding file information. If needed, more detailed evaluation information can be obtained by running test data independently.

## result description
The evaluation results were as follows:
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
    </tr>
  </thead>
  <tbody>
    <tr>
    <th>cuteSV</th>
      <th>3384</th>
      <th>2828</th>
      <th>2646</th>
      <th>576</th>
      <th>2569</th>
      <th>0.52</th>
      <th>0.64</th>
      <th>0.82</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>pbsv</th>
    	  <th>3859</th>
      <th>2677</th>
      <th>2534</th>
      <th>1240</th>
      <th>2720</th>
      <th>0.50</th>
      <th>0.67</th>
      <th>0.57</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>Sniffles</th>
      <th>3972</th>
      <th>3249</th>
      <th>3080</th>
      <th>785</th>
      <th>2148</th>
      <th>0.60</th>
      <th>0.80</th>
      <th>0.68</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>asvclr</th>
      <th>3489</th>
      <th>3086</th>
      <th>2934</th>
      <th>555</th>
      <th>2311</th>
      <th>0.57</th>
      <th>0.84</th>
      <th>0.68</th>
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
    <th>cuteSV</th>
      <th>18</th>
      <th>34</th>
      <th>53</th>
      <th>374</th>
      <th>2039</th>
      <th>253</th>
      <th>133</th>
      <th>89</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>pbsv</th>
    	  <th>25</th>
      <th>34</th>
      <th>66</th>
      <th>341</th>
      <th>2425</th>
      <th>190</th>
      <th>126</th>
      <th>93</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>Sniffles</th>
      <th>29</th>
      <th>46</th>
      <th>61</th>
      <th>260</th>
      <th>2481</th>
      <th>330</th>
      <th>155</th>
      <th>111</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>asvclr</th>
      <th>11</th>
      <th>19</th>
      <th>39</th>
      <th>342</th>
      <th>2217</th>
      <th>267</th>
      <th>133</th>
      <th>74</th>
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
    <th>cuteSV</th>
      <th>256</th>
      <th>18</th>
      <th>2782</th>
      <th>34</th>
      <th>26</th>
      <th>9</th>
      <th>28</th>
      <th>16</th>
      <th>145</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>pbsv</th>
    	  <th>293</th>
      <th>29</th>
      <th>2668</th>
      <th>43</th>
      <th>47</th>
      <th>8</th>
      <th>311</th>
      <th>82</th>
      <th>224</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>Sniffles</th>
      <th>330</th>
      <th>30</th>
      <th>3249</th>
      <th>32</th>
      <th>33</th>
      <th>5</th>
      <th>30</th>
      <th>13</th>
      <th>111</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>asvclr</th>
      <th>254</th>
      <th>20</th>
      <th>2968</th>
      <th>26</th>
      <th>23</th>
      <th>3</th>
      <th>15</th>
      <th>12</th>
      <th>2</th>
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
    </tr>
  </thead>
  <tbody>
    <tr>
    	  <th>1-100bp</th>
      <th>2780</th>
      <th>2653</th>
      <th>542</th>
      <th>1941</th>
      <th>0.59</th>
      <th>0.83</th>
      <th>0.69</th>
    </tr>
    <tr>
    	  <th>101-250bp</th>
      <th>112</th>
      <th>104</th>
      <th>19</th>
      <th>31</th>
      <th>0.46</th>
      <th>0.85</th>
      <th>0.60</th>
    </tr>
    <tr>
    	  <th>251-500bp</th>
      <th>120</th>
      <th>118</th>
      <th>2</th>
      <th>96</th>
      <th>0.56</th>
      <th>0.98</th>
      <th>0.71</th>
    </tr>
    <tr>
    	  <th>501-1000bp</th>
      <th>30</th>
      <th>30</th>
      <th>6</th>
      <th>46</th>
      <th>0.40</th>
      <th>0.83</th>
      <th>0.54</th>
    </tr>
    <tr>
    	  <th>1001-5000bp</th>
      <th>16</th>
      <th>15</th>
      <th>0</th>
      <th>76</th>
      <th>0.17</th>
      <th>1.0</th>
      <th>0.30</th>
    </tr>
    <tr>
    	  <th>5001-1000bp</th>
      <th>9</th>
      <th>9</th>
      <th>1</th>
      <th>17</th>
      <th>0.35</th>
      <th>0.90</th>
      <th>0.50</th>
    </tr>
    <tr>
    	  <th>>=10001bp</th>
      <th>3</th>
      <th>3</th>
      <th>0</th>
      <th>20</th>
      <th>0.13</th>
      <th>1.0</th>
      <th>0.23</th>
    </tr>
  </tbody>
</table>

The comprehensive evaluation results of the four detection methods were shown in the figure as follows:

<div style="text-align:center;">
<img src="img/Performance comparison between different tools.png" alt= "Performance comparison between different tools" width= "400" heiht="300" style="margin-right:25px"> 
<img src="img/Benchmark results between different tools.png" alt= "Benchmark results between different tools" width= "400" heiht="300" style="margin-left:25px;" >
</div>
