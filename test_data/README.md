# ASVBM experiments

## experiment description

In this experiment, we conducted an individualized benchmarking of the identification results for chromosome 1 (chr1) from the HG002 CCS (Consensus Circular Sequence) data by using seven structural variants (SV) detection methods: cuteSV(v2.0.3), pbsv (v2.9.0), Sniffles2(v2.0.2), ASVCLR(v1.4.4), SVDSS(v1.0.5), DeBreak(v1.0.2), SVIM(v2.0.0). Here we will only focus on the benchmarking results of insertion and deletion variants. Subsequently, the identification results obtained from these seven different detection methods were comprehensively benchmarked against their benchmark sets through a multi-sample benchmarking approach.

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
$ wget -c https://github.com/zhuxiao/asvbm/releases/download/1.1.0/asvbm_1.1.0.tar.xz
$ tar -xf asvbm_1.3.0.tar.xz
$ cd asvbm_1.3.0/
$ ./autogen.sh
$ sudo apt install -y r-base r-base-dev
# Install the UpSetR package
$ R
$ install.packages("UpSetR")
```
The name of the chr1 file (reference) was renamed to test_chr1.fa. We uploaded the sequence of chr1 into three files, which can be downloaded and merged by following steps:
```
## Merge chr1 references and extract it
$ cd asvbm_1.3.0/test_data/
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

## Benchmarking the VCF files mentioned in the test data using ASVBM, run the following command:
$ asvbm stat -m 50000 -T "ASVCLR;SVDSS;DeBreak;Sniffles2;pbsv;cuteSV;SVIM" -C "1" -p 0.7 test_data/test_ASVCLR_chr1.vcf test_data/test_SVDSS_chr1.vcf test_data/test_DeBreak_chr1.vcf test_data/test_Sniffles2_chr1.vcf test_pbsv_chr1.vcf test_cuteSV_chr1.vcf test_SVIM_chr1.vcf test_benchmark_chr1.vcf test_data/test_chr1.fa
```
ASVBM benchmarked the provided identification set based on the input, generating corresponding file information. If needed, more detailed benchmarking information can be obtained by running test data independently.

## result description
The benchmarking results were as follows:
<table>
  <thead>
    <tr>
    	  <th>SV caller</th>
    	  <th>SVs</th>
      <th>TP</th>
      <th>FP</th>
      <th>FN</th>
      <th>LP</th>
      <th>Recall</th>
      <th>precision</th>
      <th>F1 score</th>
      <th>Identity</th>
    </tr>
  </thead>
    <tbody>
    <tr>
    <th>ASVCLR</th>
      <th>3998</th>
      <th>3404</th>
      <th>593</th>
      <th>1993</th>
      <th>44</th>
      <th>0.6307</th>
      <th>0.8516</th>
      <th>0.7247</th>
      <th>0.9808</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>SVDSS</th>
      <th>3318</th>
      <th>2525</th>
      <th>543</th>
      <th>2872</th>
      <th>58</th>
      <th>0.4678</th>
      <th>0.8230</th>
      <th>0.5966</th>
      <th>0.9740</th>
    </tr>
  </tbody>
    <tbody>
    <tr>
    <th>DeBreak</th>
      <th>3429</th>
      <th>3008</th>
      <th>483</th>
      <th>2389</th>
      <th>4</th>
      <th>0.5573</th>
      <th>0.8616</th>
      <th>0.6769</th>
      <th>0.9341</th>
    </tr>
  </tbody>
     <tbody>
    <tr>
    <th>Sniffles2</th>
      <th>3972</th>
      <th>3330</th>
      <th>659</th>
      <th>2067</th>
      <th>51</th>
      <th>0.6170</th>
      <th>0.8348</th>
      <th>0.7096</th>
      <th>0.9709</th>
    </tr>
  </tbody> 
   <tbody>
    <tr>
    <th>pbsv</th>
    	<th>3859</th>
      <th>3247</th>
      <th>627</th>
      <th>2150</th>
      <th>28</th>
      <th>0.6016</th>
      <th>0.8382</th>
      <th>0.7005</th>
      <th>0.9818</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>cuteSV</th>
      <th>3384</th>
      <th>2920</th>
      <th>483</th>
      <th>2477</th>
      <th>24</th>
      <th>0.5410</th>
      <th>0.8581</th>
      <th>0.6636</th>
      <th>0.9720</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>SVIM</th>
      <th>11490</th>
      <th>3515</th>
      <th>2367</th>
      <th>1882</th>
      <th>71</th>
      <th>0.6513</th>
      <th>0.5976</th>
      <th>0.6232</th>
      <th>0.9776</th>
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
    <th>ASVCLR</th>
      <th>19</th>
      <th>18</th>
      <th>55</th>
      <th>459</th>
      <th>2388</th>
      <th>278</th>
      <th>132</th>
      <th>81</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>SVDSS</th>
      <th>15</th>
      <th>21</th>
      <th>69</th>
      <th>383</th>
      <th>2108</th>
      <th>178</th>
      <th>102</th>
      <th>75/th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>DeBreak</th>
      <th>8</th>
      <th>20</th>
      <th>48</th>
      <th>223</th>
      <th>2376</th>
      <th>281</th>
      <th>115</th>
      <th>71</th>
    </tr>
  </tbody>
     <tbody>
    <tr>
    <th>Sniffles2</th>
      <th>16</th>
      <th>30</th>
      <th>55</th>
      <th>204</th>
      <th>2522</th>
      <th>315</th>
      <th>140</th>
      <th>85</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>pbsv</th>
    	  <th>17</th>
      <th>25</th>
      <th>45</th>
      <th>337</th>
      <th>2462</th>
      <th>175</th>
      <th>84</th>
      <th>62</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>cuteSV</th>
      <th>16</th>
      <th>27</th>
      <th>45</th>
      <th>331</th>
      <th>2105</th>
      <th>230</th>
      <th>108</th>
      <th>68</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>SVIM</th>
      <th>37</th>
      <th>58</th>
      <th>79</th>
      <th>438</th>
      <th>2738</th>
      <th>273</th>
      <th>143</th>
      <th>108</th>
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
    <th>AVCRL</th>
      <th>186</th>
      <th>40</th>
      <th>3477</th>
      <th>86</th>
      <th>83</th>
      <th>15</th>
      <th>20</th>
      <th>2</th>
      <th>1</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>SVDSS</th>
      <th>178</th>
      <th>43</th>
      <th>2943</th>
      <th>78</th>
      <th>89</th>
      <th>10</th>
      <th>10</th>
      <th>0</th>
      <th>0</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>DeBreak</th>
      <th>134</th>
      <th>27</th>
      <th>3185</th>
      <th>83</th>
      <th>84</th>
      <th>11</th>
      <th>13</th>
      <th>0</th>
      <th>2</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>Sniffles2</th>
      <th>206</th>
      <th>43</th>
      <th>3446</th>
      <th>89</th>
      <th>84</th>
      <th>13</th>
      <th>13</th>
      <th>0</th>
      <th>0</th>
    </tr>
  </tbody>
   <tbody>
    <tr>
    <th>pbsv</th>
    	<th>224</th>
      <th>51</th>
      <th>3219</th>
      <th>128</th>
      <th>103</th>
      <th>16</th>
      <th>20</th>
      <th>1</th>
      <th>0</th>
    </tr>
  </tbody>
    <tbody>
    <tr>
    <th>cuteSV</th>
      <th>219</th>
      <th>39</th>
      <th>2947</th>
      <th>79</th>
      <th>71</th>
      <th>10</th>
      <th>8</th>
      <th>0</th>
      <th>1</th>
    </tr>
  </tbody>
  <tbody>
    <tr>
    <th>SVIM</th>
      <th>527</th>
      <th>64</th>
      <th>3787</th>
      <th>122</th>
      <th>128</th>
      <th>27</th>
      <th>23</th>
      <th>0</th>
      <th>25</th>
    </tr>
  </tbody>
</table>

Additionally, basic metrics for different structural variant (SV) size ranges were computed, primarily categorized into the following seven intervals. The results of ASVCLR were shown as follows:

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
      <th>Identity</th>
    </tr>
  </thead>
  <tbody>
    <tr>
    	  <th>1-100bp</th>
      <th>2862</th>
      <th>2857</th>
      <th>529</th>
      <th>1859</th>
      <th>0.6062</th>
      <th>0.8438</th>
      <th>0.7055</th>
      <th>0.9708</th>
    </tr>
    <tr>
    	  <th>101-250bp</th>
      <th>121</th>
      <th>113</th>
      <th>68</th>
      <th>122</th>
      <th>0.4979</th>
      <th>0.6243</th>
      <th>0.5540</th>
      <th>0.9915</th>
    </tr>
    <tr>
    	  <th>251-500bp</th>
      <th>129</th>
      <th>127</th>
      <th>133</th>
      <th>87</th>
      <th>0.5972</th>
      <th>0.4769</th>
      <th>0.5248</th>
      <th>0.9917</th>
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
    	  <th>1001-2500bp</th>
      <th>28</th>
      <th>33</th>
      <th>15</th>
      <th>33</th>
      <th>0.4590</th>
      <th>0.6875</th>
      <th>0.5505</th>
      <th>0.9969</th>
    </tr>
        <tr>
    	  <th>2501-5000bp</th>
      <th>14</th>
      <th>22</th>
      <th>3</th>
      <th>17</th>
      <th>0.4516</th>
      <th>0.8800</th>
      <th>0.5969</th>
      <th>0.9999</th>
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
