#include "util.h"


// string split function
vector<string> split(const string& s, const string& delim)
{
	vector<string> elems;
	string tmp;
	size_t pos = 0;
	size_t len = s.length();
	size_t delim_len = delim.length();
	int find_pos;

	if (delim_len == 0) return elems;
	while (pos < len)
	{
		find_pos = s.find(delim, pos);
		if(find_pos < 0){
			elems.push_back(s.substr(pos, len - pos));
			break;
		}
		tmp = s.substr(pos, find_pos - pos);
		if(!tmp.empty()) elems.push_back(tmp);
		pos = find_pos + delim_len;
	}
	return elems;
}
//Split a custom interval
vector<size_t> parseSemicolonSeparatedValues(const string& str){
	vector<size_t> values;
	const char* cstr = str.c_str();
	char* endptr;

	while(*cstr){
		size_t value = static_cast<size_t>(strtol(cstr, &endptr, 10));
		if(cstr == endptr){
			cerr << "Invalid value: " << cstr << endl;
			exit(1);
		}
		values.push_back(value);
		cstr = endptr;
		if(*cstr == ';'){
			cstr++;
		}
	}

	return values;
}

// transform double to string with specified fixed precision
string double2Str(double num, int32_t fixedPrecision){
	string num_str, tmp_str;
	int32_t i, dot_idx, float_bit_num, len;
	char ch;

	tmp_str = to_string(num);
	dot_idx = tmp_str.find(".");
	if(dot_idx>0){
		float_bit_num = (int32_t)tmp_str.size() - 1 - dot_idx;
		if(float_bit_num>fixedPrecision){
			ch = tmp_str.at(dot_idx+fixedPrecision+1);
			len = dot_idx + fixedPrecision + 1;
			if(ch<'5') num_str = tmp_str.substr(0, len);
			else tmp_str.at(dot_idx+fixedPrecision) ++;
		}else{
			num_str = tmp_str;
			len = fixedPrecision - float_bit_num;
			for(i=0; i<len; i++) num_str += "0";
		}
	}else{
		num_str = tmp_str;
		if(fixedPrecision>0){
			num_str += ".";
			for(i=0; i<fixedPrecision; i++) num_str += "0";
		}
	}

	return num_str;
}

// determine whether the character is a base
bool isBase(const char ch){
	bool flag = false;
	switch(ch){
		case 'a':
		case 'c':
		case 'g':
		case 't':
		case 'A':
		case 'C':
		case 'G':
		case 'T':
		case 'N':
		case 'n':

		// mixed bases
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
		case 'b': flag = true; break;
		default: cerr << __func__ << ": unknown base: " << ch << endl; exit(1);
	}

	return flag;
}

void printSV(vector<SV_item*> &data){
	string data_out, sv_type_str;
	SV_item *item;

	for(size_t i=0; i<data.size(); i++){
		item = data.at(i);

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
			default:
				cerr << "line=" << __LINE__ << ", invalid sv type: " << item->sv_type << endl;
				exit(1);
		}

		if(item->sv_type!=VAR_TRA and item->sv_type!=VAR_BND)
			data_out = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + sv_type_str + "\t" + to_string(item->sv_len);
		else
			data_out = item->chrname + "\t" + to_string(item->startPos) + "\t" + to_string(item->endPos) + "\t" + item->chrname2 + "\t" + to_string(item->startPos2) + "\t" + to_string(item->endPos2) + "\t" + sv_type_str + "\t" + to_string(item->sv_len);
		cout << data_out << endl;
	}
}

string getProgramCmdStr(int argc, char *argv[]){
	string program_cmd_str = argv[0];
	for(int i=1; i<argc; i++) program_cmd_str = program_cmd_str + " " + argv[i];	//program_cmd_str +=" " + string(argv[i]);
	return program_cmd_str;
}

std::string getAllProgramCmdStr(int argc, char *argv[]) {
    std::string program_cmd = argv[0];

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if ((arg == "-T" || arg == "-C") && i + 1 < argc) {
            program_cmd += " " + arg + " \"" + std::string(argv[i + 1]) + "\"";
            i++;
        } else {
            if (arg.find(' ') != std::string::npos || arg.find(';') != std::string::npos) {
                arg = "\"" + arg + "\"";
            }
            program_cmd += " " + arg;
        }
    }
    return program_cmd;
}

bool isDigitString(string &str){
	bool flag = true;
	for (unsigned i = 0; i < str.size(); i++)
	   {
//	    cout<<"panduan:"<<isdigit(str[i])<<endl;
	    if (isdigit(str[i])==false)
	      {
	           flag = false;
	           break;
	       }
	    }
	return flag;
}

// reverse the sequence
void reverseSeq(string &seq){
	char tmp;
	int32_t len = seq.size();
	for(int32_t i=0; i<len/2; i++){
		tmp = seq[i];
		seq[i] = seq[len-1-i];
		seq[len-1-i] = tmp;
	}
}

// reverse complement the sequence
void reverseComplement(string &seq){
	reverseSeq(seq);

	int32_t len = seq.size();
	for(int32_t i=0; i<len; i++){
		switch(seq[i]){
			case 'a': seq[i] = 't'; break;
			case 'c': seq[i] = 'g'; break;
			case 'g': seq[i] = 'c'; break;
			case 't': seq[i] = 'a'; break;
			case 'A': seq[i] = 'T'; break;
			case 'C': seq[i] = 'G'; break;
			case 'G': seq[i] = 'C'; break;
			case 'T': seq[i] = 'A'; break;
			case 'N': seq[i] = 'N'; break;
			case 'n': seq[i] = 'n'; break;

			// mixed bases
			case 'M': seq[i] = 'K'; break;
			case 'm': seq[i] = 'k'; break;
			case 'R': seq[i] = 'Y'; break;
			case 'r': seq[i] = 'y'; break;
			case 'S': seq[i] = 'S'; break;
			case 's': seq[i] = 's'; break;
			case 'V': seq[i] = 'B'; break;
			case 'v': seq[i] = 'b'; break;
			case 'W': seq[i] = 'W'; break;
			case 'w': seq[i] = 'w'; break;
			case 'Y': seq[i] = 'R'; break;
			case 'y': seq[i] = 'r'; break;
			case 'H': seq[i] = 'D'; break;
			case 'h': seq[i] = 'd'; break;
			case 'K': seq[i] = 'M'; break;
			case 'k': seq[i] = 'm'; break;
			case 'D': seq[i] = 'H'; break;
			case 'd': seq[i] = 'h'; break;
			case 'B': seq[i] = 'V'; break;
			case 'b': seq[i] = 'v'; break;
			default: cerr << __func__ << ": unknown base: " << seq[i] << endl; exit(1);
		}
	}
}
