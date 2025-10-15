#include "convert.h"
#include "util.h"
#include <typeinfo>
#include <htslib/faidx.h>

//Path query Dot
string PathqueryDot(string filename){
	string newfilename;
	size_t lastDotPos = filename.find_last_of('.');
	if(lastDotPos != string::npos and filename.substr(lastDotPos + 1).compare("bed") == 0)	return filename;
	else
		newfilename = filename.substr(0, lastDotPos) + ".bed";

	return newfilename;
}
//Path query backslash
string Pathquerybackslash(string filename){
	string newfilename;
	size_t lastSlashPos = filename.find_last_of('/');
	if (lastSlashPos != std::string::npos and lastSlashPos < filename.length() - 1) {
		newfilename =  filename.substr(lastSlashPos + 1);
		return newfilename;
	}

	return filename;
}

// convert data
/*void convertBed(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str, ref_seq, alt_seq, reg_str, lineInfo;
	vector<string> str_vec, sv_type_vec;
	size_t start_pos, endpos, start_pos2, endpos2;
	int32_t sv_len = 0, refseq_len_tmp;
	vector<int32_t> sv_len_vec;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;
	char *seq;
	faidx_t *fai;

	fai = fai_load(reffilename.c_str());

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()){
			if(line.at(0)!='#'){
				lineInfo = line;
				str_vec = split(line, "\t");

				chrname = str_vec.at(0);
				start_pos = stoi(str_vec.at(1));
				endpos = stoi(str_vec.at(2));

				// get sv type and sv_len
				sv_type_vec = getSVType(str_vec);
				//sv_len_vec = getSVLen(str_vec, sv_type_str);
				sv_len_vec = getSVLen(str_vec, sv_type_vec);

				sv_type_str = sv_type_vec.at(0);
				sv_len = sv_len_vec.at(0);

				ref_seq = str_vec.at(3);
				alt_seq = str_vec.at(4);

				if(sv_type_str.compare("DEL")==0 and sv_len_vec.at(0)>0) sv_len = -sv_len;
				if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
					chrname2 = str_vec.at(3);
					start_pos2 = stoi(str_vec.at(4));
					endpos2 = stoi(str_vec.at(5));
					ref_seq = alt_seq = "-";
				}else{
					chrname2 = "";
					start_pos2 = endpos2 = 0;
				}

				if(ref_seq.compare("N")==0 or ref_seq.compare(".")==0 or ref_seq.compare("0")==0){
					reg_str = chrname + ":" + to_string(start_pos) + "-" + to_string(endpos);
					seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
					if(strlen(seq)>=Max_SeqLen) ref_seq="-";
					else ref_seq = seq;
					free(seq);
				}
				if(alt_seq.compare("N")==0 or alt_seq.compare(".")==0 or alt_seq.compare("0")==0 or alt_seq.size()>=Max_SeqLen)
					alt_seq = "-";

				sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, ref_seq, alt_seq, lineInfo, label);
				sv_item_vec.push_back(sv_item);
			}
		}
	}
	infile.close();
	fai_destroy(fai);

	// remove mate sv items
	if(mate_filename.size()) removemateSVItems(mate_filename, sv_item_vec);

	// remove snv items
	if(snv_filename.size()) removeSNVItems(snv_filename, sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}*/

//Determine if it is present in chromosome collections
bool isExistChromosomeSet(string &chrname){
	auto it = find(chromosomeSet.begin(), chromosomeSet.end(), chrname);
	if (it != chromosomeSet.end()) 	return true;
	else	return false;
}

//To determine if it's decoy
bool isDecoyChr(string &chrname){
	bool flag = false;
	size_t pos;

	pos = chrname.find(DECOY_PREFIX);
	if(pos==0) flag = true;
	else{
		pos = chrname.find(DECOY_PREFIX2);
		if(pos==0) flag = true;
	}

	return flag;
}
void processFile(const string &infilename, bool &isVCF, bool &isBED){
	ifstream infile;
	string line;
	vector<string> headerCols;
	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}
	isVCF = false, isBED = false;
	while(getline(infile, line)){
		if(line.empty()){
			continue;
		}
		// 情况1：##开头的行（VCF元信息行）
		if(line.size() >= 2 && line[0] == '#' && line[1] == '#'){
			// 检查是否为VCF标志性元信息行
			if(line.size() >= 16 && line.substr(0, 16) == "##fileformat=VCF"){
				isVCF = true;
				break;
			}
			continue;
		}
		// 情况2：#开头但非##开头的行（标题行）
		if(line[0] == '#'){
			// 确保行长度≥2，避免访问headerCols[0]时越界（如行仅为"#"）
			if(line.size() < 2){
				continue;
			}
			headerCols = split(line, "\t");
			if(headerCols.empty()){
				continue; // 分割失败，跳过
			}
			// 判断标题行类型
			if(headerCols[0] == "#CHROM"){
				isVCF = true;
				break;
			}else if (headerCols[0] == "#Chr" || headerCols[0] == "#Chr1"){
				isBED = true;
				break;
			}
		}
	}

	infile.close();
}
/*void convertBed(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label){
	ifstream infile;
	string line, chrname, chrname2;  // 主染色体、TRA/BND涉及的第二条染色体
	string start_pos_str, endpos_str, start_pos2_str, endpos2_str;  // 位置字符串
	string sv_type_str, sv_len_str, sv_len_str1, sv_len_str2;  // 变异类型和长度
	string ref_seq, alt_seq, ref_seq1, alt_seq1, ref_seq2, alt_seq2;  // 参考/替代序列
	string reg_str, lineInfo, mate_reg;  // 基因组区域、行原始信息、mate区域
	size_t start_pos, endpos, start_pos2, endpos2;  // 位置数值（size_t适配大坐标）
	int32_t sv_len, sv_len1, sv_len2, refseq_len_tmp;  // 变异长度（支持正负）
	vector<string> str_vec;  // 行分割后的字段向量
	vector<SV_item*> sv_item_vec;  // 存储SV项目的容器
	SV_item *sv_item, *sv_item1, *sv_item2;  // SV项目指针
	char *seq;  // 从参考基因组获取的序列
	faidx_t *fai;  // 参考基因组索引（用于快速获取序列）
	vector<string> userannotationLines;  // 存储用户注释行

	fai = fai_load(reffilename.c_str());
	if(!fai){
		cerr << __func__ << ", line=" << __LINE__ << ": 无法加载参考基因组索引: " << reffilename << endl;
		exit(1);
	}
	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": 无法打开文件: " << infilename << endl;
		exit(1);
	}
	while(getline(infile, line)){
		if(line.empty()) continue;

		if(line.at(0) == '#'){
			if (label == "benchmark"){
				benchmarkannotationLines.push_back(line);  // 基准数据注释行
			}else if (label == "user"){
				userannotationLines.push_back(line);  // 用户数据注释行
			}
			continue;
		}

		// 处理数据行：非#开头的有效行
		lineInfo = line;  // 保存原始行信息，用于后续追溯
		if (label == "benchmark"){
			benchmarklineMap[lineInfo] = 0;  // 基准数据行映射（用于校验）
		}
		str_vec = split(line, "\t");  // 按制表符分割字段（BED格式强制制表符分隔）
		size_t num_cols = str_vec.size();  // 获取当前行列数


		// 校验列数：仅支持12列（非BND/TRA）或16列（BND/TRA）
		if(num_cols != 12 && num_cols != 16){
			cerr << __func__ << ": 无效BED格式（预期12或16列），行内容: " << line << endl;
			continue;
		}


		// --------------------------
		// 解析12列BED：非BND/TRA变异（如DEL、INS、INV、DUP等）
		// 格式：#Chr	Start	End	ID	SVType	SVLen	DupNum	Ref	Alt	Extra	FORMAT	sample
		// --------------------------
		if(num_cols == 12){
			chrname = str_vec[0];  // 0: 染色体名称
			start_pos_str = str_vec[1];  // 1: 起始位置（BED通常为0-based）
			endpos_str = str_vec[2];  // 2: 结束位置
			sv_type_str = str_vec[4];  // 4: 变异类型（如DEL、INS）
			sv_len_str = str_vec[5];  // 5: 变异长度
			ref_seq = str_vec[7];  // 7: 参考序列（REF）
			alt_seq = str_vec[8];  // 8: 替代序列（ALT）

			// 染色体过滤（排除decoy和无效染色体）
			if(!chromosomeSet.empty() && !isExistChromosomeSet(chrname)) continue;
			if(chromosomeSet.empty() && isDecoyChr(chrname)) continue;

			// 解析变异长度（校验有效性）
			if(sv_len_str.empty() || sv_len_str.find_first_not_of("-0123456789") != string::npos){
				sv_len = 0;  // 无效长度标记为0
			}else{
				sv_len = stoi(sv_len_str);
				// 过滤过长无效长度
				if(abs(sv_len) > MAX_VALID_SVLEN){
					sv_len = 0;
				}
			}

			// 处理参考序列（REF）：缺失时从参考基因组获取
			if(ref_seq == "N" || ref_seq == "." || ref_seq == "0"){
				start_pos = stoull(start_pos_str);
				endpos = stoull(endpos_str);
				if (endpos < start_pos) endpos = start_pos;  // 避免无效区间
				reg_str = chrname + ":" + start_pos_str + "-" + to_string(endpos);
				// 从参考基因组索引提取序列
				seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
				if(seq == nullptr || strlen(seq) >= Max_SeqLen){
					ref_seq = "-";  // 序列无效或过长，标记为-
				}else{
					ref_seq = seq;  // 有效序列赋值
				}
				free(seq);  // 释放fai_fetch分配的内存
			}

			// 处理替代序列（ALT）：过长或无效时标记为-
			if(alt_seq.size() >= Max_SeqLen || alt_seq == "."){
				alt_seq = "-";
			}

			// 转换位置为数值型
			start_pos = stoull(start_pos_str);
			endpos = stoull(endpos_str);
			// TRA/BND相关字段置空（非此类变异）
			chrname2 = "";
			start_pos2 = endpos2 = 0;
			sv_len1 = sv_len2 = 0;

			// 特殊处理DEL类型：长度取绝对值（同convertVcf逻辑）
			if(sv_type_str == "DEL"){
				sv_len = abs(sv_len);
			}

			// 创建SV项目并添加到容器
			sv_item = allocateSVItem(
					chrname, start_pos, endpos, chrname2, start_pos2, endpos2,
					sv_type_str, sv_len, ref_seq, alt_seq, lineInfo, label
			);
			sv_item_vec.push_back(sv_item);
		}


        // --------------------------
        // 解析16列BED：BND/TRA变异（涉及两条染色体的易位或断裂点）
        // 格式：#Chr1	Start1	End1	Chr2	Start2	End2	ID	SVType	SVLen1	SVLen2	MateReg	Ref1	Alt1	Ref2	Alt2	FORMAT	sample
        // --------------------------
        else if (num_cols == 16) {
            // 提取核心字段（列索引固定）
            chrname = str_vec[0];  // 0: 第一条染色体
            start_pos_str = str_vec[1];  // 1: 第一条染色体起始位置
            endpos_str = str_vec[2];  // 2: 第一条染色体结束位置
            chrname2 = str_vec[3];  // 3: 第二条染色体（TRA/BND涉及）
            start_pos2_str = str_vec[4];  // 4: 第二条染色体起始位置
            endpos2_str = str_vec[5];  // 5: 第二条染色体结束位置
            sv_type_str = str_vec[7];  // 7: 变异类型（BND或TRA）
            sv_len_str1 = str_vec[8];  // 8: 第一条染色体侧长度
            sv_len_str2 = str_vec[9];  // 9: 第二条染色体侧长度
            mate_reg = str_vec[10];  // 10: mate区域信息
            ref_seq1 = str_vec[11];  // 11: 第一条染色体侧REF
            alt_seq1 = str_vec[12];  // 12: 第一条染色体侧ALT
            ref_seq2 = str_vec[13];  // 13: 第二条染色体侧REF
            alt_seq2 = str_vec[14];  // 14: 第二条染色体侧ALT

            // 染色体过滤：两条染色体均需有效
            if ((!chromosomeSet.empty() && (!isExistChromosomeSet(chrname) || !isExistChromosomeSet(chrname2))) ||
                (chromosomeSet.empty() && (isDecoyChr(chrname) || isDecoyChr(chrname2)))) {
                continue;
            }

            // 解析变异长度（两条染色体侧）
            sv_len1 = (sv_len_str1.empty() || sv_len_str1.find_first_not_of("-0123456789") != string::npos) ? 0 : stoi(sv_len_str1);
            sv_len2 = (sv_len_str2.empty() || sv_len_str2.find_first_not_of("-0123456789") != string::npos) ? 0 : stoi(sv_len_str2);
            if (abs(sv_len1) > MAX_VALID_SVLEN) sv_len1 = 0;
            if (abs(sv_len2) > MAX_VALID_SVLEN) sv_len2 = 0;

            // 处理第一条染色体侧参考序列（REF1）
            if (ref_seq1 == "N" || ref_seq1 == "." || ref_seq1 == "0") {
                start_pos = stoull(start_pos_str);
                endpos = stoull(endpos_str);
                if (endpos < start_pos) endpos = start_pos;
                reg_str = chrname + ":" + start_pos_str + "-" + to_string(endpos);
                seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
                ref_seq1 = (seq == nullptr || strlen(seq) >= Max_SeqLen) ? "-" : seq;
                free(seq);
            }

            // 处理第二条染色体侧参考序列（REF2）
            if (ref_seq2 == "N" || ref_seq2 == "." || ref_seq2 == "0") {
                start_pos2 = stoull(start_pos2_str);
                endpos2 = stoull(endpos2_str);
                if (endpos2 < start_pos2) endpos2 = start_pos2;
                reg_str = chrname2 + ":" + start_pos2_str + "-" + to_string(endpos2);
                seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
                ref_seq2 = (seq == nullptr || strlen(seq) >= Max_SeqLen) ? "-" : seq;
                free(seq);
            }

            // 处理替代序列（ALT1和ALT2）
            if (alt_seq1.size() >= Max_SeqLen || alt_seq1 == ".") alt_seq1 = "-";
            if (alt_seq2.size() >= Max_SeqLen || alt_seq2 == ".") alt_seq2 = "-";

            // 转换位置为数值型
            start_pos = stoull(start_pos_str);
            endpos = stoull(endpos_str);
            start_pos2 = stoull(start_pos2_str);
            endpos2 = stoull(endpos2_str);

            // 对BND/TRA变异创建两个SV项目（分别对应两条染色体侧）
            sv_item1 = allocateSVItem(
                chrname, start_pos, endpos, chrname2, start_pos2, endpos2,
                sv_type_str, sv_len1, ref_seq1, alt_seq1, lineInfo, label
            );
            sv_item2 = allocateSVItem(
                chrname2, start_pos2, endpos2, chrname, start_pos, endpos,  // 反向关联主染色体
                sv_type_str, sv_len2, ref_seq2, alt_seq2, lineInfo, label
            );
            sv_item_vec.push_back(sv_item1);
            sv_item_vec.push_back(sv_item2);
        }
    }



	infile.close();
	fai_destroy(fai);
	if(!mate_filename.empty()) removemateSVItems(mate_filename, sv_item_vec);
	if(!snv_filename.empty()) removeSNVItems(snv_filename, sv_item_vec);

	output2File(outfilename, sv_item_vec, outConvertScreenFile);
	releaseSVItemVec(sv_item_vec);

	usersetsannotationLines.push_back(userannotationLines);
	usersets_num++;
}*/
void convertBed(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label) {
    ifstream infile;
    string line, lineInfo;

    // 存储两行表头的列名（原始列名和小写列名）

    vector<string> header1_lower, header2_lower;  // 小写列名（用于不区分大小写匹配）
    int header_count = 0;  // 已读取的表头行数（0/1/2）

    // 第一行表头（非TRA/BND）的字段索引结构体
    struct NormalIndex {
        size_t chr = (size_t)-1;    // #Chr列索引
        size_t start = (size_t)-1;  // Start列索引
        size_t end = (size_t)-1;    // End列索引
        size_t svtype = (size_t)-1; // SVType列索引
        size_t svlen = (size_t)-1;  // SVLen列索引
        size_t ref = (size_t)-1;    // Ref列索引
        size_t alt = (size_t)-1;    // Alt列索引
    } normal_idx;

    // 第二行表头（TRA/BND）的字段索引结构体
    struct TraIndex {
        size_t chr1 = (size_t)-1;   // #Chr1列索引
        size_t start1 = (size_t)-1; // Start1列索引
        size_t end1 = (size_t)-1;   // End1列索引
        size_t chr2 = (size_t)-1;   // Chr2列索引
        size_t start2 = (size_t)-1; // Start2列索引
        size_t end2 = (size_t)-1;   // End2列索引
        size_t svtype = (size_t)-1; // SVType列索引
        size_t svlen1 = (size_t)-1; // SVLen1列索引
        size_t svlen2 = (size_t)-1; // SVLen2列索引
        size_t ref1 = (size_t)-1;   // Ref1列索引
        size_t alt1 = (size_t)-1;   // Alt1列索引
        size_t ref2 = (size_t)-1;   // Ref2列索引
        size_t alt2 = (size_t)-1;   // Alt2列索引
    } tra_idx;

    // 其他变量定义
    string chrname, chrname2;
    string start_pos_str, endpos_str, start1_str, end1_str, start2_str, end2_str;
    string sv_type_str, sv_len_str, sv_len1_str, sv_len2_str;
    string ref_seq, alt_seq, ref1_seq, alt1_seq, ref2_seq, alt2_seq;
    string reg_str;
    size_t start_pos, endpos, start1, end1, start2, end2;
    int32_t sv_len = 0, sv_len1 = 0, sv_len2 = 0, refseq_len_tmp;
    vector<string> str_vec;
    vector<SV_item*> sv_item_vec;
    SV_item *sv_item, *sv_item1, *sv_item2;
    char *seq;
    faidx_t *fai;
    vector<string> userannotationLines;
    string empty_str;

    // 加载参考基因组索引
    fai = fai_load(reffilename.c_str());
    if (!fai) {
        cerr << __func__ << ", line=" << __LINE__ << ": 无法加载参考基因组索引: " << reffilename << endl;
        exit(1);
    }

    // 打开输入BED文件
    infile.open(infilename);
    if (!infile.is_open()) {
        cerr << __func__ << ", line=" << __LINE__ << ": 无法打开文件: " << infilename << endl;
        exit(1);
    }

    // 逐行处理BED文件
    while (getline(infile, line)) {
        if (line.empty()) continue;

        // 1. 处理注释行（支持两行表头）
        if (line.at(0) == '#' && header_count < 2) {
            header_count++;
            vector<string> cols = split(line, "\t");
            vector<string> cols_lower;

            // 转换列名为小写，用于不区分大小写匹配
            for (const auto &col : cols) {
            	string col_clean = col;
            	    // 直接移除所有回车符（ASCII=13）
            	    col_clean.erase(remove(col_clean.begin(), col_clean.end(), '\r'), col_clean.end());

            	    // 再移除其他控制字符（ASCII < 32）
            	    string col_filtered;
            	    for (char c : col_clean) {
            	        if (c >= 32) {
            	            col_filtered += c;
            	        }
            	    }

            	    // 转换为小写
            	    string col_lower;
            	    transform(col_filtered.begin(), col_filtered.end(), back_inserter(col_lower), ::tolower);
            	    if (header_count == 1) {
            	            header1_lower.push_back(col_lower);  // 第一行表头的清洗后小写列名
            	        } else if (header_count == 2) {
            	            header2_lower.push_back(col_lower);  // 第二行表头的清洗后小写列名
            	        }
            }
            // 解析第一行表头（非TRA/BND）
            if (header_count == 1) {


                // 提取非TRA/BND字段索引
                for (size_t i = 0; i < header1_lower.size(); ++i) {
                    if (header1_lower[i] == "#chr" || header1_lower[i] == "chr") {
                        normal_idx.chr = i;
                    } else if (header1_lower[i] == "start") {
                        normal_idx.start = i;
                    } else if (header1_lower[i] == "end") {
                        normal_idx.end = i;
                    } else if (header1_lower[i] == "svtype") {
                        normal_idx.svtype = i;
                    } else if (header1_lower[i] == "svlen") {
                        normal_idx.svlen = i;
                    } else if (header1_lower[i] == "ref") {
                        normal_idx.ref = i;
                    } else if (header1_lower[i] == "alt") {
                        normal_idx.alt = i;
                    }
                }

                // 验证第一行表头必填列
                if (normal_idx.chr == (size_t)-1 || normal_idx.start == (size_t)-1 || normal_idx.end == (size_t)-1) {
                    cerr << __func__ << ": 第一行表头缺少非TRA/BND必填列（Chr/Start/End）" << endl;
                    exit(1);
                }

            }
            // 解析第二行表头（TRA/BND）
            else if (header_count == 2) {


                // 提取TRA/BND字段索引
                for (size_t i = 0; i < header2_lower.size(); ++i) {
                    if (header2_lower[i] == "#chr1" || header2_lower[i] == "chr1") {
                        tra_idx.chr1 = i;
                    } else if (header2_lower[i] == "start1") {
                        tra_idx.start1 = i;
                    } else if (header2_lower[i] == "end1") {
                        tra_idx.end1 = i;
                    } else if (header2_lower[i] == "chr2") {
                        tra_idx.chr2 = i;
                    } else if (header2_lower[i] == "start2") {
                        tra_idx.start2 = i;
                    } else if (header2_lower[i] == "end2") {
                        tra_idx.end2 = i;
                    } else if (header2_lower[i] == "svtype") {
                        tra_idx.svtype = i;
                    } else if (header2_lower[i] == "svlen1") {
                        tra_idx.svlen1 = i;
                    } else if (header2_lower[i] == "svlen2") {
                        tra_idx.svlen2 = i;
                    } else if (header2_lower[i] == "ref1") {
                        tra_idx.ref1 = i;
                    } else if (header2_lower[i] == "alt1") {
                        tra_idx.alt1 = i;
                    } else if (header2_lower[i] == "ref2") {
                        tra_idx.ref2 = i;
                    } else if (header2_lower[i] == "alt2") {
                        tra_idx.alt2 = i;
                    }
                }

                // 验证第二行表头必填列
                if (tra_idx.chr1 == (size_t)-1 || tra_idx.start1 == (size_t)-1 || tra_idx.end1 == (size_t)-1 ||
                    tra_idx.chr2 == (size_t)-1 || tra_idx.start2 == (size_t)-1 || tra_idx.end2 == (size_t)-1) {
                    cerr << __func__ << ": 第二行表头缺少TRA/BND必填列（Chr1/Start1/End1/Chr2/Start2/End2）" << endl;
                    exit(1);
                }
            }

            // 保存注释行
            if (label == "benchmark") {
                benchmarkannotationLines.push_back(line);
            } else if (label == "user") {
                userannotationLines.push_back(line);
            }
            continue;
        }

        // 2. 处理数据行（需至少读取1行表头）
        if (header_count == 0) {
            cerr << __func__ << ": 数据行前未找到标题行，格式错误" << endl;
            continue;
        }

        lineInfo = line;
        if (label == "benchmark") {
            benchmarklineMap[lineInfo] = 0;
        }
        str_vec = split(line, "\t");
        size_t num_cols = str_vec.size();

        // 3. 动态匹配数据行对应的表头（根据列数）
        bool is_normal = (header_count >= 1 && num_cols == header1_lower.size());  // 匹配第一行表头（非TRA/BND）
        bool is_tra = (header_count >= 2 && num_cols == header2_lower.size());    // 匹配第二行表头（TRA/BND）
        // 4. 处理非TRA/BND变异（匹配第一行表头）
        if (is_normal) {
            // 检查必填列是否有效
            if (normal_idx.chr >= num_cols || normal_idx.start >= num_cols || normal_idx.end >= num_cols) {
                cerr << __func__ << ": 非TRA/BND数据行缺少必填列，跳过: " << line << endl;
                continue;
            }

            // 提取字段（使用第一行表头索引）
            chrname = str_vec[normal_idx.chr];
            start_pos_str = str_vec[normal_idx.start];
            endpos_str = str_vec[normal_idx.end];
            sv_type_str = (normal_idx.svtype != (size_t)-1 && normal_idx.svtype < num_cols) ? str_vec[normal_idx.svtype] : "-";
            sv_len_str = (normal_idx.svlen != (size_t)-1 && normal_idx.svlen < num_cols) ? str_vec[normal_idx.svlen] : "-";
            ref_seq = (normal_idx.ref != (size_t)-1 && normal_idx.ref < num_cols) ? str_vec[normal_idx.ref] : "-";
            alt_seq = (normal_idx.alt != (size_t)-1 && normal_idx.alt < num_cols) ? str_vec[normal_idx.alt] : "-";
            // 染色体过滤
            if (!chromosomeSet.empty() && !isExistChromosomeSet(chrname)) continue;
            if (chromosomeSet.empty() && isDecoyChr(chrname)) continue;

            // 解析SVLen（核心修复：正确提取非TRA/BND的SVLen值）
            sv_len = 0;
            sv_len_str.erase(remove(sv_len_str.begin(), sv_len_str.end(), '\r'), sv_len_str.end());
            if (sv_len_str != "-") {
                // 检查是否为有效数字
                if (sv_len_str.find_first_not_of("-0123456789") == string::npos) {
                    try {
                        sv_len = stoi(sv_len_str);
                        // 过滤超过阈值的长度
                        if (abs(sv_len) > MAX_VALID_SVLEN) {
                            sv_len = 0;
                        }
                    } catch (...) {
                        sv_len = 0;  // 转换失败时设为0
                    }
                }
                // 否则保持0（非数字）
            }
            // 处理参考序列
            if (ref_seq == "N" || ref_seq == "." || ref_seq == "0" || ref_seq == "-") {
                try {
                    start_pos = stoull(start_pos_str);
                    endpos = stoull(endpos_str);
                } catch (...) {
                    cerr << __func__ << ": 位置格式错误，跳过: " << line << endl;
                    continue;
                }
                if (endpos < start_pos) endpos = start_pos;
                reg_str = chrname + ":" + start_pos_str + "-" + to_string(endpos);
                seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
                ref_seq = (seq == nullptr || strlen(seq) >= Max_SeqLen) ? "-" : seq;
                free(seq);
            }

            // 处理替代序列
            if (alt_seq.size() >= Max_SeqLen || alt_seq == "." || alt_seq == "0") {
                alt_seq = "-";
            }

            // 转换位置为数值型
            try {
                start_pos = stoull(start_pos_str);
                endpos = stoull(endpos_str);
            } catch (...) {
                cerr << __func__ << ": 位置转换失败，跳过: " << line << endl;
                continue;
            }

            // DEL类型长度取绝对值
            if (sv_type_str != "-" && (sv_type_str == "DEL" || sv_type_str == "del")) {
                sv_len = abs(sv_len);
            }

            // 创建SV_item
            empty_str.clear();
            sv_item = allocateSVItem(
                chrname, start_pos, endpos, empty_str, 0, 0,
                sv_type_str, sv_len, ref_seq, alt_seq, lineInfo, label
            );
            sv_item_vec.push_back(sv_item);



        }
        // 5. 处理TRA/BND变异（匹配第二行表头）
        else if (is_tra) {
            // 检查必填列是否有效
            if (tra_idx.chr1 >= num_cols || tra_idx.start1 >= num_cols || tra_idx.end1 >= num_cols ||
                tra_idx.chr2 >= num_cols || tra_idx.start2 >= num_cols || tra_idx.end2 >= num_cols) {
                cerr << __func__ << ": TRA/BND数据行缺少必填列，跳过: " << line << endl;
                continue;
            }

            // 提取字段（使用第二行表头索引）
            chrname = str_vec[tra_idx.chr1];
            start1_str = str_vec[tra_idx.start1];
            end1_str = str_vec[tra_idx.end1];
            chrname2 = str_vec[tra_idx.chr2];
            start2_str = str_vec[tra_idx.start2];
            end2_str = str_vec[tra_idx.end2];
            sv_type_str = (tra_idx.svtype != (size_t)-1 && tra_idx.svtype < num_cols) ? str_vec[tra_idx.svtype] : "-";
            sv_len1_str = (tra_idx.svlen1 != (size_t)-1 && tra_idx.svlen1 < num_cols) ? str_vec[tra_idx.svlen1] : "-";
            sv_len2_str = (tra_idx.svlen2 != (size_t)-1 && tra_idx.svlen2 < num_cols) ? str_vec[tra_idx.svlen2] : "-";
            ref1_seq = (tra_idx.ref1 != (size_t)-1 && tra_idx.ref1 < num_cols) ? str_vec[tra_idx.ref1] : "-";
            alt1_seq = (tra_idx.alt1 != (size_t)-1 && tra_idx.alt1 < num_cols) ? str_vec[tra_idx.alt1] : "-";
            ref2_seq = (tra_idx.ref2 != (size_t)-1 && tra_idx.ref2 < num_cols) ? str_vec[tra_idx.ref2] : "-";
            alt2_seq = (tra_idx.alt2 != (size_t)-1 && tra_idx.alt2 < num_cols) ? str_vec[tra_idx.alt2] : "-";

            // 染色体过滤
            if ((!chromosomeSet.empty() && (!isExistChromosomeSet(chrname) || !isExistChromosomeSet(chrname2))) ||
                (chromosomeSet.empty() && (isDecoyChr(chrname) || isDecoyChr(chrname2)))) {
                continue;
            }

            // 解析SVLen1和SVLen2
            sv_len1 = 0;
            if (sv_len1_str != "-" && sv_len1_str.find_first_not_of("-0123456789") == string::npos) {
                try {
                    sv_len1 = stoi(sv_len1_str);
                    if (abs(sv_len1) > MAX_VALID_SVLEN) sv_len1 = 0;
                } catch (...) {
                    sv_len1 = 0;
                }
            }

            sv_len2 = 0;
            if (sv_len2_str != "-" && sv_len2_str.find_first_not_of("-0123456789") == string::npos) {
                try {
                    sv_len2 = stoi(sv_len2_str);
                    if (abs(sv_len2) > MAX_VALID_SVLEN) sv_len2 = 0;
                } catch (...) {
                    sv_len2 = 0;
                }
            }

            // 处理参考序列1
            if (ref1_seq == "N" || ref1_seq == "." || ref1_seq == "0" || ref1_seq == "-") {
                try {
                    start1 = stoull(start1_str);
                    end1 = stoull(end1_str);
                } catch (...) {
                    cerr << __func__ << ": Chr1位置格式错误，跳过: " << line << endl;
                    continue;
                }
                if (end1 < start1) end1 = start1;
                reg_str = chrname + ":" + start1_str + "-" + to_string(end1);
                seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
                ref1_seq = (seq == nullptr || strlen(seq) >= Max_SeqLen) ? "-" : seq;
                free(seq);
            }

            // 处理参考序列2
            if (ref2_seq == "N" || ref2_seq == "." || ref2_seq == "0" || ref2_seq == "-") {
                try {
                    start2 = stoull(start2_str);
                    end2 = stoull(end2_str);
                } catch (...) {
                    cerr << __func__ << ": Chr2位置格式错误，跳过: " << line << endl;
                    continue;
                }
                if (end2 < start2) end2 = start2;
                reg_str = chrname2 + ":" + start2_str + "-" + to_string(end2);
                seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
                ref2_seq = (seq == nullptr || strlen(seq) >= Max_SeqLen) ? "-" : seq;
                free(seq);
            }

            // 处理替代序列
            if (alt1_seq.size() >= Max_SeqLen || alt1_seq == "." || alt1_seq == "0") alt1_seq = "-";
            if (alt2_seq.size() >= Max_SeqLen || alt2_seq == "." || alt2_seq == "0") alt2_seq = "-";

            // 转换位置为数值型
            try {
                start1 = stoull(start1_str);
                end1 = stoull(end1_str);
                start2 = stoull(start2_str);
                end2 = stoull(end2_str);
            } catch (...) {
                cerr << __func__ << ": 位置转换失败，跳过: " << line << endl;
                continue;
            }

            // 创建两个SV_item（对应两个断点）
            sv_item1 = allocateSVItem(
                chrname, start1, end1, chrname2, start2, end2,
                sv_type_str, sv_len1, ref1_seq, alt1_seq, lineInfo, label
            );
            sv_item2 = allocateSVItem(
                chrname2, start2, end2, chrname, start1, end1,
                sv_type_str, sv_len2, ref2_seq, alt2_seq, lineInfo, label
            );
            sv_item_vec.push_back(sv_item1);
            sv_item_vec.push_back(sv_item2);

        }
        // 6. 列数不匹配任何表头
        else {
            cerr << __func__ << ": 数据行列数（" << num_cols << "）不匹配任何表头（第一行：" << header1_lower.size()
                 << "列，第二行：" << header2_lower.size() << "列），跳过: " << line << endl;
            continue;
        }
    }

    // 后续处理
    infile.close();
    fai_destroy(fai);

    if (!mate_filename.empty()) removemateSVItems(mate_filename, sv_item_vec);
    if (!snv_filename.empty()) removeSNVItems(snv_filename, sv_item_vec);

    output2File(outfilename, sv_item_vec, outConvertScreenFile);
    releaseSVItemVec(sv_item_vec);

    usersetsannotationLines.push_back(userannotationLines);
    usersets_num++;
}


// convert vcf result
void convertVcf(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label){
	ifstream infile;
	string line, chrname, start_pos_str, endpos_str, endpos_str_1, endpos_str_2, chrname2, start_pos_str2, endpos_str2, sv_type_str, sv_type_str1, sv_type_str2, sv_len_str, sv_len_str1, sv_len_str2, bnd_str, bnd_pos_str, str_tmp;
	string ref_seq, alt_seq, reg_str, lineInfo;
	size_t start_pos, endpos, endpos_1, endpos_2, start_pos2, endpos2;
	size_t i, end_pos_ref, end_pos_ref1, end_pos_ref2;
	int32_t sv_len, sv_len1, sv_len2, ref_seq_len, refseq_len_tmp;
	bool is_svlen_flag;
	vector<int32_t> sv_len_vec;
	vector<string> str_vec, sv_type_vec, info_vec, sub_info_vec, bnd_pos_vec;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item, *sv_item1, *sv_item2;
	char *seq;
	faidx_t *fai;
	vector<string> userannotationLines;

	fai = fai_load(reffilename.c_str());
	is_svlen_flag = false;

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()){
			if(line.at(0)!='#'){
				lineInfo = line;

				if(label.compare("benchmark")==0)	benchmarklineMap[lineInfo] = 0;
				str_vec = split(line, "\t");

				// get locations
				chrname = str_vec.at(0);
				if(chromosomeSet.size()>0){
					if(isExistChromosomeSet(chrname) == false) continue;
				}else{
					if(isDecoyChr(chrname) == true) continue;
				}
				start_pos_str = str_vec.at(1);
//				if(start_pos_str!="39059")
//					continue;
//				cout<<"start_pos_str:"<<start_pos_str<<endl;
				ref_seq_len = str_vec.at(3).size();
				ref_seq = str_vec.at(3);
				alt_seq = str_vec.at(4);

				if(alt_seq.at(0) == '<') alt_seq = "-";

				chrname2 = endpos_str = sv_type_str = sv_len_str = "";
				info_vec = split(str_vec.at(7), ";");

				for(i=0; i<info_vec.size(); i++){
					sub_info_vec = split(info_vec.at(i), "=");
					if(sub_info_vec.at(0).compare("CHR2")==0)
						chrname2 = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("END")==0)
						endpos_str = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("SVTYPE")==0)
						sv_type_str = sub_info_vec.at(1);
					if(sub_info_vec.at(0).compare("SVLEN")==0){
						sv_len_str = sub_info_vec.at(1);

						if(sv_len_str.find_first_not_of("-0123456789")!=string::npos) sv_len_str = "";
						else{ // ignore too large invalid variant length, e.g. SVLEN=999999999 in Sniffles output
							sv_len = stoi(sv_len_str);
							is_svlen_flag = true;
							if(sv_len>MAX_VALID_SVLEN) { sv_len_str = ""; sv_len = 0;}
						}
					}
				}

				//Skip mutations with incomplete information
//				if(sv_type_str.size()==0 and sv_len_str.size()==0 and endpos_str.size()==0)	continue;

				// get sv type and length
				if(sv_type_str.size()==0){
					sv_type_vec = getSVType(str_vec);
					if(sv_type_vec.size()==1){
						sv_type_str = sv_type_vec.at(0);
					}else{
						sv_type_str1 = sv_type_vec.at(0);
						sv_type_str2 = sv_type_vec.at(1);
					}
				}
				if(sv_len_str.size()==0){
					sv_len_vec = getSVLen(str_vec, sv_type_str);
					if(sv_len_vec.size()==1){
						sv_len_str = to_string(sv_len_vec.at(0));
					}else{
						if(sv_len_vec.size()!=0){
							sv_len_str1 = to_string(sv_len_vec.at(0));
							sv_len_str2 = to_string(sv_len_vec.at(1));
						}else continue;
					}
				}
				if(endpos_str.size()==0){
					if(sv_len_vec.size()==1){
							end_pos_ref = stoi(start_pos_str) + ref_seq_len - 1;
							endpos_str = to_string(end_pos_ref);
					}else{
						if(sv_len_str.size()>0)	endpos_str = to_string(stoi(start_pos_str) + stoi(sv_len_str) - 1);
						else endpos_str = to_string(stoi(start_pos_str) + abs(sv_len_vec.at(0)) - 1);
						end_pos_ref1 = stoi(start_pos_str) + ref_seq_len - 1;
						end_pos_ref2 = stoi(start_pos_str) + ref_seq_len - 1;
						endpos_str_1 = to_string(end_pos_ref1);
						endpos_str_2 = to_string(end_pos_ref2);
					}
				}

				if(sv_type_str.size()>0){
					if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){ // TRA or BND
						if(sv_type_str.compare("BND")==0){ // BND
							bnd_str = str_vec.at(4);
							if(bnd_str.at(0)==']' or bnd_str.at(0)=='[')
								bnd_pos_str = bnd_str.substr(1, bnd_str.size()-3);
							else if(isBase(bnd_str.at(0)))
								bnd_pos_str = bnd_str.substr(2, bnd_str.size()-3);
							bnd_pos_vec = split(bnd_pos_str, ":");
							chrname2 = bnd_pos_vec.at(0);
							endpos_str = bnd_pos_vec.at(1);
						}

						start_pos_str2 = endpos_str;
						endpos_str2 = start_pos_str2;
						endpos_str = start_pos_str;
					}else { // INS, DEL, INV, DUP, UNC
						start_pos_str2 = "0";
						endpos_str2 = "0";
					}
					if(sv_len_vec.size()==1 or is_svlen_flag == true){
						start_pos = stoi(start_pos_str);
						endpos = stoi(endpos_str);
						start_pos2 = stoi(start_pos_str2);
						endpos2 = stoi(endpos_str2);
						if(sv_len_str.size()>0) sv_len = stoi(sv_len_str);
						else sv_len = 0;

						if(sv_type_str.compare("DEL")==0) sv_len = abs(sv_len);

						if(ref_seq.compare("N")==0 or ref_seq.compare(".")==0 or ref_seq.compare("0")==0){
							if(endpos < (size_t)stoi(start_pos_str))	endpos = (size_t)stoi(start_pos_str);
							reg_str = chrname + ":" + start_pos_str + "-" + to_string(endpos);
							seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
							if(seq == nullptr or strlen(seq) >= Max_SeqLen ) ref_seq = "-";
							else ref_seq = seq;
							free(seq);
						}
						if(alt_seq.size()>=Max_SeqLen or sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0 or alt_seq.compare(".")==0) alt_seq="-";
						sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, ref_seq, alt_seq, lineInfo, label);
						sv_item_vec.push_back(sv_item);





					}else{
						start_pos = stoi(start_pos_str);
						endpos_1 = stoi(endpos_str_1);
						endpos_2 = stoi(endpos_str_2);
						start_pos2 = stoi(start_pos_str2);
						endpos2 = stoi(endpos_str2);
						sv_len1 = stoi(sv_len_str1);
						sv_len2 = stoi(sv_len_str2);

						if(sv_type_str1.compare("DEL")==0) sv_len1 = abs(sv_len1);
						if(sv_type_str2.compare("DEL")==0) sv_len2 = abs(sv_len2);

						ref_seq = alt_seq = "-";

						sv_item1 = allocateSVItem(chrname, start_pos, endpos_1, chrname2, start_pos2, endpos2, sv_type_str1, sv_len1, ref_seq, alt_seq, lineInfo, label);
						sv_item2 = allocateSVItem(chrname, start_pos, endpos_2, chrname2, start_pos2, endpos2, sv_type_str2, sv_len2, ref_seq, alt_seq, lineInfo, label);
						sv_item_vec.push_back(sv_item1);
						sv_item_vec.push_back(sv_item2);


					}
				}else{
					cout << line << endl;
					cout << "missing SVTYPE information" << endl;
				}
			}else{
				if(label.compare("benchmark")==0)	benchmarkannotationLines.push_back(line);
				if(label.compare("user")==0)	{ userannotationLines.push_back(line);}
			}
		}
	}
	infile.close();
	fai_destroy(fai);

	// remove mate sv items
	if(mate_filename.size()) removemateSVItems(mate_filename, sv_item_vec);

	// remove snv items
	if(snv_filename.size()) removeSNVItems(snv_filename, sv_item_vec);


	//if(mate_filename.size()==0)	fp_vec.push_back(sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
	usersetsannotationLines.push_back(userannotationLines);
	usersets_num += 1;


}

// convert data
void convertCsv(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str, ref_seq, alt_seq, reg_str, lineInfo;
	vector<string> str_vec, sv_type_vec;
	size_t start_pos, endpos, start_pos2, endpos2;
	int32_t sv_len = 0, refseq_len_tmp;
	vector<int32_t> sv_len_vec;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;
	char *seq;
	faidx_t *fai;

	fai = fai_load(reffilename.c_str());

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()>0){
			lineInfo = line;
			if(line.at(0)!='#'){
				str_vec = split(line, ",");
				if(str_vec.at(0).size()>0){

					chrname = str_vec.at(0);
					if(isDigitString(str_vec.at(1))==false or isDigitString(str_vec.at(2))==false)
						continue;
					start_pos = stoi(str_vec.at(1));
					endpos = stoi(str_vec.at(2));

					// get sv type and length
					sv_type_vec = getSVType(str_vec);
					sv_len_vec = getSVLen(str_vec, sv_type_vec);

					ref_seq = str_vec.at(3);
					alt_seq = str_vec.at(4);

					sv_type_str = sv_type_vec.at(0);
					sv_len = sv_len_vec.at(0);
					//if(sv_type_str.compare("DEL")==0 and sv_len_vec.at(0)>0) sv_len = -sv_len;
					//if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
					if(sv_type_str.compare("DEL")==0 and sv_len_vec.at(0)>0) sv_len = -sv_len;
					if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
						chrname2 = str_vec.at(3);
						start_pos2 = stoi(str_vec.at(4));
						endpos2 = stoi(str_vec.at(5));
					}else{
						chrname2 = "";
						start_pos2 = endpos2 = 0;
					}

					if(ref_seq.compare("N")==0 or ref_seq.compare(".")==0 or ref_seq.compare("0")==0){
						reg_str = chrname + ":" + to_string(start_pos) + "-" + to_string(endpos);
						seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
						if(strlen(seq)>=Max_SeqLen) ref_seq="-";
						else ref_seq = seq;
						free(seq);
					}
					if(alt_seq.compare("N")==0 or alt_seq.compare(".")==0 or alt_seq.compare("0")==0 or alt_seq.size()>Max_SeqLen)
						alt_seq = "-";

					//sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len);
					sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, ref_seq, alt_seq, lineInfo, label);
					sv_item_vec.push_back(sv_item);
				}
			}
		}
	}
	infile.close();

	// remove mate sv items
	if(mate_filename.size()) removemateSVItems(mate_filename, sv_item_vec);

	// remove snv items
	if(snv_filename.size()) removeSNVItems(snv_filename, sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// convert data
void convertNm(const string &infilename, const string &outfilename, const string &reffilename, string &mate_filename, string &snv_filename, string &label){
	ifstream infile;
	string line, chrname, chrname2, sv_type_str, ref_seq, alt_seq, reg_str, lineInfo;
	vector<string> str_vec, sv_type_vec;
	int32_t sv_len, ref_seq_len, start_pos, endpos, start_pos2, endpos2, refseq_len_tmp;
	vector<int32_t> sv_len_vec;
	vector<SV_item*> sv_item_vec;
	SV_item *sv_item;
	char *seq;
	faidx_t *fai;

	fai = fai_load(reffilename.c_str());

	sv_len = 0;

	infile.open(infilename);
	if(!infile.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << infilename << endl;
		exit(1);
	}

	// convert
	while(getline(infile, line)){
		if(line.size()>0){
			lineInfo = line;
			if(line.at(0)!='#'){
				str_vec = split(line, "\t");

				chrname = str_vec.at(0);
				start_pos = stoi(str_vec.at(1));

				ref_seq_len = str_vec.at(3).size();
				endpos = stoi(str_vec.at(1)) + ref_seq_len - 1; // endPos

				ref_seq = str_vec.at(3);
				alt_seq = str_vec.at(4);

				if(alt_seq.at(0) == '<') alt_seq = "-";

				// get sv type
				sv_type_vec = getSVType(str_vec);

				// get sv length
				//sv_len_vec = getSVLen(str_vec, sv_type_str);
				sv_len_vec = getSVLen(str_vec, sv_type_vec);

				sv_type_str = sv_type_vec.at(0);
				sv_len = sv_len_vec.at(0);
				if(sv_type_str.compare("DEL")==0 and sv_len_vec.at(0)>0) sv_len = -sv_len;
				if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("BND")==0){
					chrname2 = str_vec.at(3);
					start_pos2 = stoi(str_vec.at(4));
					endpos2 = stoi(str_vec.at(5));
				}else{
					chrname2 = "";
					start_pos2 = endpos2 = 0;
				}
//				chrname2 = "";
//				start_pos2 = endpos2 = 0;

				if(sv_type_str.compare("DEL")==0) sv_len = -sv_len;

				if(ref_seq.compare("N")==0 or ref_seq.compare(".")==0 or ref_seq.compare("0")==0){
					reg_str = chrname + ":" + to_string(start_pos) + "-" + to_string(endpos);
					seq = fai_fetch(fai, reg_str.c_str(), &refseq_len_tmp);
					if(strlen(seq)>=Max_SeqLen) ref_seq="-";
					else ref_seq = seq;
					free(seq);
				}
				if(alt_seq.size()>=Max_SeqLen) alt_seq="-";

				sv_item = allocateSVItem(chrname, start_pos, endpos, chrname2, start_pos2, endpos2, sv_type_str, sv_len, ref_seq, alt_seq, lineInfo, label);
				sv_item_vec.push_back(sv_item);
			}
		}
	}
	infile.close();
	fai_destroy(fai);

	// remove mate sv items
	if(mate_filename.size()) removemateSVItems(mate_filename, sv_item_vec);

	// remove snv items
	if(snv_filename.size()) removeSNVItems(snv_filename, sv_item_vec);

	// save to file
	output2File(outfilename, sv_item_vec, outConvertScreenFile);

	// release memory
	releaseSVItemVec(sv_item_vec);
}

// get sv type
vector<string> getSVType(vector<string> &str_vec){
	string sv_type_str, sv_type_str1, sv_type_str2, str_tmp;
	size_t i, j;
	bool is_seq_flag_ref, is_seq_flag_alt, is_comma_flag;
	vector<string> comma_str_vec;
	vector<string> info_vec, sv_type_vec;
	vector<int32_t> sv_len_vec;
	int32_t sv_type_pos, SVTYPE_pos, vec_idx, str_pos, start_pos, end_pos;

	sv_type_str = "";
	for(i=3; i<str_vec.size(); i++){
		str_tmp = str_vec.at(i);
		if(str_tmp.compare("INS")==0 or str_tmp.compare("DEL")==0 or str_tmp.compare("DUP")==0 or str_tmp.compare("INV")==0 or str_tmp.compare("TRA")==0 or str_tmp.compare("BND")==0 or str_tmp.compare("INVTRA")==0 or str_tmp.compare("MIX")==0 or str_tmp.compare("MNP")==0 or str_tmp.compare("SNV")==0){
			sv_type_str = str_tmp;
		}else if(str_tmp.compare("<INS>")==0 or str_tmp.compare("<DEL>")==0 or str_tmp.compare("<DUP>")==0 or str_tmp.compare("<INV>")==0 or str_tmp.compare("<TRA>")==0 or str_tmp.compare("<BND>")==0 or str_tmp.compare("<INVTRA>")==0 or str_tmp.compare("<MIX>")==0 or str_tmp.compare("<MNP>")==0 or str_tmp.compare("<SNV>")==0){
			sv_type_str = str_tmp.substr(1, str_tmp.size()-2);
		}else if(str_tmp.compare("insertion")==0 or str_tmp.compare("deletion")==0 or str_tmp.compare("duplication")==0 or str_tmp.compare("inversion")==0 or str_tmp.compare("translocation")==0 or str_tmp.compare("snv")==0){
			sv_type_str = str_tmp;
		}
		if(sv_type_str.size()>0) {
			sv_type_vec.push_back(sv_type_str);
			break;
		}
	}

	if(sv_type_str.size()==0){
		vec_idx = -1; str_pos = -1;
		for(i=3; i<str_vec.size(); i++){
			str_tmp = str_vec.at(i);
			sv_type_pos = str_tmp.find("sv_type=");
			if(sv_type_pos!=-1){
				vec_idx = i;
				str_pos = sv_type_pos + 8;
				break;
			}else{
				SVTYPE_pos = str_tmp.find("SVTYPE=");
				if(SVTYPE_pos!=-1){
					vec_idx = i;
					str_pos = SVTYPE_pos + 7;
					break;
				}
			}
		}

		if(vec_idx!=-1 and str_pos!=-1){ // found
			str_tmp = str_vec.at(vec_idx);
			start_pos = str_pos; end_pos = -1;
			for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
				if(str_tmp.at(j)==';'){
					end_pos = j - 1;
					break;
				}
			}
			if(end_pos==-1) end_pos = str_tmp.size() - 1;
			sv_type_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
			if(sv_type_str.size()>0)
				sv_type_vec.push_back(sv_type_str);
		}
	}

	if(sv_type_str.size()==0){// not found
		// check sequence
		is_seq_flag_ref = isSeq(str_vec.at(3));
		if(is_seq_flag_ref){
			is_seq_flag_alt = isSeq(str_vec.at(4));
			if(is_seq_flag_alt){
				if(str_vec.at(3).size()==1 and str_vec.at(4).size()==1)sv_type_str = "SNV";
				if(str_vec.at(3).size()<str_vec.at(4).size())sv_type_str = "INS";
				if(str_vec.at(3).size()>str_vec.at(4).size())sv_type_str = "DEL";
				sv_type_vec.push_back(sv_type_str);
			}else{
				str_tmp = str_vec.at(4);
				is_comma_flag = isComma(str_vec.at(4));
				if(is_comma_flag){
					comma_str_vec = split(str_vec.at(4), ",");
					if(str_vec.at(3).size()<comma_str_vec.at(0).size())sv_type_str1 = "INS";
					if(str_vec.at(3).size()>comma_str_vec.at(0).size())sv_type_str1 = "DEL";
					if(str_vec.at(3).size()<comma_str_vec.at(1).size())sv_type_str2 = "INS";
					if(str_vec.at(3).size()>comma_str_vec.at(1).size())sv_type_str2 = "DEL";
					sv_type_vec.push_back(sv_type_str1);
					sv_type_vec.push_back(sv_type_str2);
				}
			}
		}
	}

	return sv_type_vec;
}

// get SV length
vector<int32_t> getSVLen(vector<string> &str_vec, string &sv_type){
	int32_t sv_len, sv_len1, sv_len2, vec_idx, str_pos, sv_length_pos, SVLEN_pos, start_pos, end_pos, sv_END_pos;
	size_t i, j;
	string str_tmp, sv_len_str;
	bool flag, is_seq_flag_ref, is_seq_flag_alt, is_comma_flag;
	vector<string> comma_str_vec;
	vector<int32_t> sv_len_vec;

	flag = false;
	sv_len = 0;
	if(sv_type.compare("TRA")!=0 and sv_type.compare("BND")!=0){ // not TRA
		// sv_length, or SVLEN
		vec_idx = str_pos = -1;
		for(i=3; i<str_vec.size(); i++){
			str_tmp = str_vec.at(i);
			sv_length_pos = str_tmp.find("sv_length=");
			if(sv_length_pos!=-1){
				vec_idx = i;
				str_pos = sv_length_pos + 10;
				break;
			}else{
				SVLEN_pos = str_tmp.find("SVLEN=");
				if(SVLEN_pos!=-1){
					vec_idx = i;
					str_pos = SVLEN_pos + 6;
					break;
				}
			}
		}

		if(vec_idx!=-1 and str_pos!=-1){ // found
			str_tmp = str_vec.at(vec_idx);
			start_pos = str_pos; end_pos = -1;
			for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
				if(str_tmp.at(j)==';'){
					end_pos = j - 1;
					break;
				}
			}
			if(end_pos==-1) end_pos = str_tmp.size() - 1;
			sv_len_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
			if(sv_len_str.find_first_not_of("-0123456789")==string::npos) {
				sv_len = stoi(sv_len_str);
				sv_len_vec.push_back(sv_len);
				flag = true;
				if(sv_len>=MAX_VALID_SVLEN) { sv_len_str = ""; sv_len = 0; flag = false; }
			}else { sv_len_str = ""; flag = false; }
		}

		// 'SVLEN' was not detected, then further detect 'END'
		if(flag==false){
			vec_idx = str_pos = -1;
			for(i=3; i<str_vec.size(); i++){
				str_tmp = str_vec.at(i);
				sv_END_pos = str_tmp.find("END=");
				if(sv_END_pos!=-1){
					vec_idx = i;
					str_pos = sv_END_pos + 4;
					break;
				}
			}

			if(vec_idx!=-1 and str_pos!=-1){ // found
				str_tmp = str_vec.at(vec_idx);
				start_pos = str_pos; end_pos = -1;
				for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
					if(str_tmp.at(j)==';'){
						end_pos = j - 1;
						break;
					}
				}
				if(end_pos==-1) end_pos = str_tmp.size() - 1;
				sv_len_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
				if(sv_len_str.find_first_not_of("-0123456789")==string::npos) { sv_len = stoi(sv_len_str) - stoi(str_vec.at(1)) + 1; sv_len_vec.push_back(sv_len); flag = true; }
			}else{ // not found
				// check sequence
				is_seq_flag_ref = isSeq(str_vec.at(3));
				if(is_seq_flag_ref){
					is_seq_flag_alt = isSeq(str_vec.at(4));
					if(is_seq_flag_alt){
						sv_len = (int32_t)str_vec.at(4).size() - (int32_t)str_vec.at(3).size();
						sv_len_vec.push_back(sv_len);
						flag = true;
					}else{
						str_tmp = str_vec.at(4);
						is_comma_flag = isComma(str_vec.at(4));
						if(is_comma_flag){
							comma_str_vec = split(str_vec.at(4), ",");
							sv_len1 = (int32_t)comma_str_vec.at(0).size() - (int32_t)str_vec.at(3).size();
							sv_len2 = (int32_t)comma_str_vec.at(1).size() - (int32_t)str_vec.at(3).size();
							sv_len_vec.push_back(sv_len1);
							sv_len_vec.push_back(sv_len2);
							flag = true;
						}else{
							if(str_tmp.at(0)=='<' and str_tmp.at(str_tmp.size()-1)=='>'){
								sv_len = 0;
								flag = true;
							}
						}
					}
				}

				// no sequence
				if(flag==false){
					start_pos = stoi(str_vec.at(1));
					end_pos = stoi(str_vec.at(2));
					sv_len = end_pos - start_pos + 1;
					sv_len_vec.push_back(sv_len);
					flag = true;
				}
			}
		}
		if(flag==false) cout << "missing SVLEN information" << endl;
	}else sv_len_vec.push_back(sv_len);

	return sv_len_vec;
}

// get SV length
vector<int32_t> getSVLen(vector<string> &str_vec, vector<string> &sv_type_vec){
	int32_t sv_len, sv_len1, sv_len2, vec_idx, str_pos, sv_length_pos, SVLEN_pos, start_pos, end_pos, sv_END_pos;
	size_t i, j, n;
	string str_tmp, sv_len_str, sv_type_str;
	bool flag, is_seq_flag_ref, is_seq_flag_alt, is_comma_flag;
	vector<string> comma_str_vec;
	vector<int32_t> sv_len_vec;

	for(n=0; n<sv_type_vec.size(); n++){
		flag = false;
		sv_len = 0;
		sv_type_str = sv_type_vec.at(n);
		if(sv_type_str.compare("TRA")!=0 and sv_type_str.compare("BND")!=0){ // not TRA
			// sv_length, or SVLEN
			vec_idx = str_pos = -1;
			for(i=3; i<str_vec.size(); i++){
				str_tmp = str_vec.at(i);
				sv_length_pos = str_tmp.find("sv_length=");
				if(sv_length_pos!=-1){
					vec_idx = i;
					str_pos = sv_length_pos + 10;
					break;
				}else{
					SVLEN_pos = str_tmp.find("SVLEN=");
					if(SVLEN_pos!=-1){
						vec_idx = i;
						str_pos = SVLEN_pos + 6;
						break;
					}
				}
			}

			if(vec_idx!=-1 and str_pos!=-1){ // found
				str_tmp = str_vec.at(vec_idx);
				start_pos = str_pos; end_pos = -1;
				for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
					if(str_tmp.at(j)==';'){
						end_pos = j - 1;
						break;
					}
				}
				if(end_pos==-1) end_pos = str_tmp.size() - 1;
				sv_len_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
				if(sv_len_str.find_first_not_of("-0123456789")==string::npos) {
					sv_len = stoi(sv_len_str);
					sv_len_vec.push_back(sv_len);
					flag = true;
					if(sv_len>=MAX_VALID_SVLEN) { sv_len_str = ""; sv_len = 0; flag = false; }
				}else { sv_len_str = ""; flag = false; }
			}

			// 'SVLEN' was not detected, then further detect 'END'
			if(flag==false){
				vec_idx = str_pos = -1;
				for(i=3; i<str_vec.size(); i++){
					str_tmp = str_vec.at(i);
					sv_END_pos = str_tmp.find("END=");
					if(sv_END_pos!=-1){
						vec_idx = i;
						str_pos = sv_END_pos + 4;
						break;
					}
				}

				if(vec_idx!=-1 and str_pos!=-1){ // found
					str_tmp = str_vec.at(vec_idx);
					start_pos = str_pos; end_pos = -1;
					for(j=start_pos; j<str_tmp.size(); j++){ // search ';'
						if(str_tmp.at(j)==';'){
							end_pos = j - 1;
							break;
						}
					}
					if(end_pos==-1) end_pos = str_tmp.size() - 1;
					sv_len_str = str_tmp.substr(start_pos, end_pos-start_pos+1);
					if(sv_len_str.find_first_not_of("-0123456789")==string::npos) { sv_len = stoi(sv_len_str) - stoi(str_vec.at(1)) + 1; sv_len_vec.push_back(sv_len); flag = true; }
				}else{ // not found
					// check sequence
					is_seq_flag_ref = isSeq(str_vec.at(3));
					if(is_seq_flag_ref){
						is_seq_flag_alt = isSeq(str_vec.at(4));
						if(is_seq_flag_alt){
							sv_len = (int32_t)str_vec.at(4).size() - (int32_t)str_vec.at(3).size();
							sv_len_vec.push_back(sv_len);
							flag = true;
						}else{
							str_tmp = str_vec.at(4);
							is_comma_flag = isComma(str_vec.at(4));
							if(is_comma_flag){
								comma_str_vec = split(str_vec.at(4), ",");
								sv_len1 = (int32_t)comma_str_vec.at(0).size() - (int32_t)str_vec.at(3).size();
								sv_len2 = (int32_t)comma_str_vec.at(1).size() - (int32_t)str_vec.at(3).size();
								sv_len_vec.push_back(sv_len1);
								sv_len_vec.push_back(sv_len2);
								flag = true;
							}else{
								if(str_tmp.at(0)=='<' and str_tmp.at(str_tmp.size()-1)=='>'){
									sv_len = 0;
									sv_len_vec.push_back(sv_len);
									flag = true;
								}
							}
						}
					}

					// no sequence
					if(flag==false){
						start_pos = stoi(str_vec.at(1));
						end_pos = stoi(str_vec.at(2));
						sv_len = end_pos - start_pos + 1;
						sv_len_vec.push_back(sv_len);
						flag = true;
					}
				}
			}
			if(flag==false) cout << "missing SVLEN information" << endl;
		}else sv_len_vec.push_back(sv_len);
	}

	return sv_len_vec;
}

bool isSeq(string &seq){
	bool flag = true;
	bool comma_flag = false;
	char ch;
	for(size_t i=0; i<seq.size(); i++){
		ch = seq.at(i);
		switch(ch){
			case 'A':
			case 'a':
			case 'C':
			case 'c':
			case 'G':
			case 'g':
			case 'T':
			case 't':
			case 'N':
			case 'n':
			case 'M':
			case 'm':
			case 'R':
			case 'r':
			case 'S':
			case 's':
			case 'V':
			case 'v':
			case 'W':
			case 'w':
			case 'Y':
			case 'y':
			case 'H':
			case 'h':
			case 'K':
			case 'k':
			case 'D':
			case 'd':
			case 'B':
			case 'b':
			case '*':
				//flag = true;
				break;
			case ',':
				flag = false;
				comma_flag = true;
				break;
			default:
				flag = false;
				break;
			//default: cerr << __func__ << ": unknown base: " << ch << endl; exit(1);
		}
		if(comma_flag) flag = false;
		if(!flag) break;
	}
	//if(comma_flag) flag = false;
	return flag;
}

bool isComma(string &seq){
	bool flag = false;
	char ch;
	for(size_t i=0; i<seq.size(); i++){
		ch = seq.at(i);
		if(ch==','){
			flag = true;
			break;
		}
	}
	return flag;
}

// allocate SV item
SV_item *allocateSVItem(string &chrname, size_t startPos, size_t endPos, string &chrname2, size_t startPos2, size_t endPos2, string &sv_type_str, int32_t sv_len, string &ref_seq, string &alt_seq, string &line, string &tab){
	size_t sv_type;

	if(sv_type_str.compare("UNC")==0){
		sv_type = VAR_UNC;
	}else if(sv_type_str.compare("INS")==0 or sv_type_str.compare("insertion")==0){
		sv_type = VAR_INS;
	}else if(sv_type_str.compare("DEL")==0 or sv_type_str.compare("deletion")==0){
		sv_type = VAR_DEL;
	}else if(sv_type_str.compare("DUP")==0 or sv_type_str.compare("duplication")==0 or sv_type_str.compare("DUP:TANDEM")==0){
		sv_type = VAR_DUP;
	}else if(sv_type_str.compare("INV")==0 or sv_type_str.compare("inversion")==0){
		sv_type = VAR_INV;
	}else if(sv_type_str.compare("TRA")==0 or sv_type_str.compare("translocation")==0){
		sv_type = VAR_TRA;
	}else if(sv_type_str.compare("BND")==0){
		sv_type = VAR_BND;
	}else if(sv_type_str.compare("INVTRA")==0){
		sv_type = VAR_INV_TRA;
	}else if(sv_type_str.compare("MIX")==0){
		sv_type = VAR_MIX;
	}else if(sv_type_str.compare("MNP")==0){
		sv_type = VAR_MNP;
	}else if(sv_type_str.compare("SNV")==0){
		sv_type = VAR_SNV;
	}else if(sv_type_str.compare("CNV")==0 or sv_type_str.compare("cnv")==0){
		sv_type = VAR_CNV;
	}else{
		sv_type = VAR_UNC;
	}

	SV_item *item = new SV_item();
	item->chrname = chrname;
	item->startPos = startPos;
	item->endPos = endPos;
	item->chrname2 = chrname2;
	item->startPos2 = startPos2;
	item->endPos2 = endPos2;
	item->sv_type = sv_type;
	item->sv_len = sv_len;
	item->ref_seq = ref_seq;
	item->alt_seq = alt_seq;
	item->overlapped = false;
	item->validFlag = true;
	item->seqcons = "-";
	if(tab.compare("benchmark")==0)	item->lineInfo = line;
	else if(tab.compare("user")==0)	item->lineInfo = line;
	else item->lineInfo = "-";
	return item;
}

// remove mate SV items
void removemateSVItems(string &mate_filename, vector<SV_item*> &sv_item_vec){
	size_t i, j, mate_num;
	SV_item *item, *item2;
	ofstream redudant_file;
	string line, sv_type_str;

	cout << ">>>>>>>>> Remove duplicated mate variant items <<<<<<<<<" << endl;
	cout << "Before removing duplicated mate variant items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << ">>>>>>>>> Remove duplicated mate variant items <<<<<<<<<" << endl;
	outConvertScreenFile << "Before removing duplicated mate variant items, data size: " << sv_item_vec.size() << endl;

	for(i=0; i<sv_item_vec.size(); i++){
		item = sv_item_vec.at(i);
		if(item->validFlag and item->chrname.size()>0 and item->chrname2.size()>0){
			for(j=i+1; j<sv_item_vec.size(); j++){
				item2 = sv_item_vec.at(j);
				if(item2->validFlag and item2->chrname.size()>0 and item2->chrname2.size()>0)
					if(item->chrname.compare(item2->chrname2)==0 and item->startPos==item2->startPos2 and item->endPos==item2->endPos2
						and item->chrname2.compare(item2->chrname)==0 and item->startPos2==item2->startPos and item->endPos2==item2->endPos)
						item2->validFlag = false;
			}
		}
	}
	//outputPathname
	mate_filename = suboutputDirnamePath + "/" + mate_filename;
	redudant_file.open(mate_filename);
	if(!redudant_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << mate_filename << endl;
		exit(1);
	}

	line = "#chr\tstartPos\tendPos\tSVType\tSVLen\tRef\tAlt";
	redudant_file << line << endl;

	// remove items
	mate_num = 0;
	for(i=0; i<sv_item_vec.size(); ){
		item = sv_item_vec.at(i);
		if(item->validFlag==false){
			switch(item->sv_type){
				case VAR_UNC: sv_type_str = "UNC"; break;
				case VAR_INS: sv_type_str = "INS"; break;
				case VAR_DEL: sv_type_str = "DEL"; break;
				case VAR_DUP: sv_type_str = "DUP"; break;
				case VAR_INV: sv_type_str = "INV"; break;
				case VAR_TRA: sv_type_str = "TRA"; break;
				case VAR_BND: sv_type_str = "BND"; break;
				case VAR_INV_TRA: sv_type_str = "INVTRA"; break;
				case VAR_MIX: sv_type_str = "MIX"; break;
				case VAR_SNV: sv_type_str = "SNV"; break;
				default:
					cerr << "line=" << __LINE__ << ", invalid sv type: " << item->sv_type << endl;
					exit(1);
			}
			line = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + sv_type_str + "\t" + to_string(item->sv_len) + "\t" + item->ref_seq + "\t" + item->alt_seq;
			redudant_file << line << endl;
			mate_num ++;

			// delete mate items
			sv_item_vec.erase(sv_item_vec.begin()+i);
			delete item;
		}else i++;
	}
	redudant_file.close();

	cout << "After removing duplicated mate variant items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << "After removing duplicated mate variant items, data size: " << sv_item_vec.size() << endl;
	cout << mate_num << " duplicated mate items were saved to " << mate_filename << endl;
	outConvertScreenFile << mate_num << " duplicated mate items were saved to " << mate_filename << endl;
}

//remove SNV
void removeSNVItems(string &snv_filename, vector<SV_item*> &sv_item_vec){
	size_t i, snv_num;
	SV_item *item;
	ofstream snv_file;
	string line, sv_type_str;

	cout << ">>>>>>>>> Remove snv items <<<<<<<<<" << endl;
	cout << "Before removing snv items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << ">>>>>>>>> Remove snv items <<<<<<<<<" << endl;
	outConvertScreenFile << "Before removing snv items, data size: " << sv_item_vec.size() << endl;

//	snv_filename = outputPathname + snv_filename;
	snv_filename = suboutputDirnamePath + "/" + snv_filename;
	snv_file.open(snv_filename);
	if(!snv_file.is_open()){
		cerr << __func__ << ", line=" << __LINE__ << ": cannot open file:" << snv_filename << endl;
		exit(1);
	}

	line = "#chr\tstartPos\tendPos\tSVType\tSVLen\tRef\tAlt";
	snv_file << line << endl;

	// remove items
	snv_num = 0;
	for(i=0; i<sv_item_vec.size(); ){
		item = sv_item_vec.at(i);
		if(item->sv_len == 0 and item->sv_type == VAR_SNV){
			sv_type_str = "SNV";
			line = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + sv_type_str + "\t" + to_string(item->sv_len) + "\t" + item->ref_seq + "\t" + item->alt_seq;

			snv_file << line << endl;
			snv_num ++;

			// delete snv items
			sv_item_vec.erase(sv_item_vec.begin()+i);
			delete item;
		}else i++;
	}
	snv_file.close();

	cout << "After removing snv items, data size: " << sv_item_vec.size() << endl;
	outConvertScreenFile << "After removing snv items, data size: " << sv_item_vec.size() << endl;
	cout << snv_num << " snv items were saved to " << snv_filename << endl;
	outConvertScreenFile << snv_num << " snv items were saved to " << snv_filename << endl;
}

// release sv item vector
void releaseSVItemVec(vector<SV_item*> &sv_item_vec){
	SV_item *item;
	for(size_t i=0; i<sv_item_vec.size(); i++){
		item = sv_item_vec.at(i);
		delete item;
	}
	vector<SV_item*>().swap(sv_item_vec);
}

