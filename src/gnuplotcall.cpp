#include "gnuplotcall.h"

void replaceUnderscoreWithDot(string& str) {
	for (char& c : str) {
		if (c == '_') {
			c = '.';
		}
	}
}

string getContentAfterSlash(const string& inputString) {
	size_t slashPosition = inputString.find('/');
	if (slashPosition != string::npos) {
		return inputString.substr(slashPosition + 1);
	}
	else {
		return inputString;
	}
}

void ResultPresentation(vector<string> &sv_files1, string &outputPathname, vector<string> &tool_names, string &outputBasicMetricschart, vector< vector<float> > MeticsValues, vector< vector<int> > MeticsValues1, string &ref_file){

	cout << endl <<"############# Stage 5: Generate comparison graph information: #############" << endl << endl;
	if(sv_files1.size()>1){
//		string newInfo = "tool Recall Precision F1-score Seqcons";
//		string newInfo1 = "tool TP_user TP_benchmark FP FN";
		outputBasicMetricschartPath = outputPathname + outputBasicMetricschart;
		mkdir(outputBasicMetricschartPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		multipledataset(MeticsValues, sv_files1, tool_names, outputBasicMetricschartPath);
		multipledataset(MeticsValues1, sv_files1, tool_names, outputBasicMetricschartPath);

		outputDiffRangeBasicMetricschart = outputBasicMetricschartPath + '/' + outputDiffRangeBasicMetricschart;
		mkdir(outputDiffRangeBasicMetricschart.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		//Draw multiple single graph patterns
//		ComparisonofMetricsindifferentrangesOp(outputDiffRangeBasicMetricschart, tool_names, sv_files1);
		//Draw a multi-graph mode
		vector<string> filenames;
		GenerateFileOp(outputDiffRangeBasicMetricschart, tool_names, filenames, sv_files1);
		GenerateMultiBarCharts(outputDiffRangeBasicMetricschart, filenames, tool_names, sv_files1);

		//sv information comparison of different sizes
		SVsizeratiofile = outputBasicMetricschartPath + '/' + SVsizeratiofile;
		mkdir(SVsizeratiofile.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		GenerateSVsizeRatioFileGraph(sv_files1, tool_names, SVsizeratiofile);
		//Rscript.R
		//outputBasicMetricschartPath (figures 路径)
		RscriptGeneration(outputBasicMetricschartPath);
		//UpSetR
		outputUpSetRchart = outputBasicMetricschartPath + '/' + outputUpSetRchart;
		mkdir(outputUpSetRchart.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		TPfilemanipulation(outputUpSetRchart);
		//run Rscript.R
		RunRscriptoperation(outputBasicMetricschartPath, outputUpSetRchart);

		//Common FN
		outputCommonFN = outputPathname + outputCommonFN;
		mkdir(outputCommonFN.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		GenerateSharedFNfile(outputCommonFN);
		//shared FP
		if(refine_bench_flag){
			string mate_filename = "", snv_filename = "";
			string label1 = "user", label2 = "benchmark";
			vector<SV_item*> sv_item_vec;
			string benchfilename = outputCommonFN + "/" + "benchmarkfile.bed";
			convertVcf(outputBenchfilename, benchfilename, ref_file, mate_filename, snv_filename, label2);

			for(size_t i=0; i<sharedFPFilenames.size(); i++){
				SVcallernames[i] = SVcallernames[i] + ".bed";
				string outputSVcaller = outputCommonFN + "/" + SVcallernames[i];
				convertVcf(sharedFPFilenames[i], outputSVcaller, ref_file, mate_filename, snv_filename, label1);
				sv_item_vec = loadData(outputSVcaller);
				fp_vec.push_back(sv_item_vec);
				if(SVcallernames[i].size() >= 4 && SVcallernames[i].substr(SVcallernames[i].size() - 4) == ".bed"){
					SVcallernames[i] = SVcallernames[i].substr(0, SVcallernames[i].size() - 4);
				}
			}
			string outputSharedFP = outputCommonFN + "/";
			findSharedFP(outputSharedFP, ref_file, benchfilename, sv_files1);
		}
		if(create_bench_flag){
		string outputSharedFP = outputCommonFN + "/";
		integrationBenchmark(outputSharedFP, ref_file, sv_files1);
		}
//		findCommonFN(outputCommonFN, FNfilesPath);
		//SV distribution plot
		SVdistributionDirname = outputBasicMetricschartPath + '/' + SVdistributionDirname;
		mkdir(SVdistributionDirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		//plot
		SvNumberDistributionGraph(regSizeFiles, SVdistributionDirname);
		//html
		Generatehtml(outputBasicMetricschart);
	}else{
		outputBasicMetricschartPath = outputPathname + outputBasicMetricschart;
		mkdir(outputBasicMetricschartPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		Histogram_drawing(MeticsValues,  outputPathname, outputBasicMetricschartPath);
		Histogram_drawing(MeticsValues1,  outputPathname, outputBasicMetricschartPath);
		Generatehtml(outputBasicMetricschart);
	}
	cout << endl << "Details of the benchmarking results are saved in: " << outputPathname + htmlFilename << endl;

	cout << endl << "## More information ## " << endl;
	cout << "For more detailed benchmarking results, please refer to the generated result information in the respective folders." << endl;
	cout << "For more detailed experiment information, please refer to the github repositories: https://github.com/zhuxiao/asvbm and https://github.com/zhuxiao/asvbm-experiments." << endl;
//	cout << "For more detailed evaluation results, please refer to the generated result information in the respective folders." << endl;
	cout << "If you have any problems, comments, or suggestions, please contact xzhu@ytu.edu.cn without hesitation. Thank you very much!" << endl;
}

//Compare multiple data sets
void multipledataset(vector< vector<float> > MeticsValues, vector<string> &sv_files1, vector<string> &tool_names, string &outputBasicMetricschart){
	// Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
	vector<string> scenarios;
	vector<string> metrics = { "SeqSim", "Recall", "Precision", "F1-score"};
	if(tool_names.size()>1) scenarios = tool_names;
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			replaceUnderscoreWithDot(str);
			scenarios.push_back(str);
		}
	}
	// Create a data file to hold the data
	string filename = "benchmarking_result";
	string filenamePath = outputBasicMetricschart + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file: " << filenamePath << endl;
		return ;
	}

    // Write data to a file
    /*for (size_t i = 0; i < scenarios.size(); i++) {
        dataFile << scenarios[i];
        for (size_t j = 0; j < 4; j++) {
            dataFile << " " << MeticsValues[i][j];
        }
        dataFile << endl;
    }*/
	for (size_t i = 0; i < 4; i++) {
		dataFile << metrics[i];
		for (size_t j = 0; j < scenarios.size(); j++) {
			dataFile << " " << MeticsValues[j][i];
		}
		dataFile << endl;
	}
	dataFile.close();

	string Info;
	if(tool_names.size()>1)	Info = "Tool";
	else	Info = "callsets";
	string plotCommand ="plot '" + filenamePath + "'";
	for(size_t i = 0; i<scenarios.size(); i++){
		Info += " " + scenarios[i];
		if(i==0)
			plotCommand += " using " + to_string(i+2) + ":xtic(1) title '" + scenarios[i] + "'";
		else
			plotCommand += " '' using " + to_string(i+2) + " title '" + scenarios[i] + "'";
		if(i<scenarios.size()-1)
			plotCommand += ",";
	}

	// Set the file name to save
	string outputFileName = "benchmarking_result.png";
	string outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;

	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}

	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set title 'Performance comparison between different tools' font 'Times-Roman,14'\n");
	fprintf(gnuplotPipe, "set xlabel 'Metrics' font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set ylabel 'Percentage' font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,10'\n");
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");
//    fprintf(gnuplotPipe, "set key outside right top vertical font 'Times-Roman,10'\n");
//	fprintf(gnuplotPipe, "set rmargin 12\n");

	// Set key (legend) font size
	fprintf(gnuplotPipe, "set key font 'Times-Roman,10'\n");
	// Draw a bar chart, specifying the number of columns for each data column purple  olive magenta
//    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb '#568AC6' title 'Recall', '' using 3 lc rgb '#F08700' title 'Precision', '' using 4 lc rgb '#F06767' title 'F1-score', '' using 5 lc rgb '#BBBB6D' title 'SeqSim'\n", filenamePath.c_str());
	fprintf(gnuplotPipe, "%s\n", plotCommand.c_str());
	// Close the GNUplot pipeline
	pclose(gnuplotPipe);

	AddfileInformation(filenamePath, Info);
	cout << "The performance metrics result between the user callsets is saved in: " << outputFileNamePath << endl;

}

//stat 5
void multipledataset(vector< vector<int> > MeticsValues, vector<string> &sv_files1, vector<string> &tool_names, string &outputBasicMetricschart){
	// Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
	vector<string> scenarios;
	vector<string> metrics = {"TP", "FP", "FN", "LP"};
	if(tool_names.size()>1) scenarios = tool_names;
	//Note to be modified to remove '/'
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			replaceUnderscoreWithDot(str);
			scenarios.push_back(str);
		}
	}
	// Create a data file to hold the data
	string filename = "result_classification";
	string filenamePath = outputBasicMetricschart + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file: " << filenamePath << endl;
		return ;
	}
	// Write data to a file
	for (size_t i = 0; i < 4; i++) {
		dataFile << metrics[i];
		for (size_t j = 0; j < scenarios.size(); j++) {
			dataFile << " " << MeticsValues[j][i];
		}
		dataFile << endl;
	}
	dataFile.close();
	string Info = "Tool";
	string plotCommand ="plot '" + filenamePath + "'";
	for(size_t i = 0; i<scenarios.size(); i++){
		Info += " " + scenarios[i];
		if(i==0)
			plotCommand += " using " + to_string(i+2) + ":xtic(1) title '" + scenarios[i] + "'";
		else
			plotCommand += " '' using " + to_string(i+2) + " title '" + scenarios[i] + "'";
		if(i<scenarios.size()-1)
			plotCommand += ",";
	}

	// Set the file name to save
	string outputFileName = "result_classification.png";
	string outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;

	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}

	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set title 'Benchmark results between different tools' font 'Times-Roman,14'\n");
	fprintf(gnuplotPipe, "set xlabel 'Metrics' font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set ylabel 'Count' font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,10'\n");
	fprintf(gnuplotPipe, "set yrange [0:*]\n");
	fprintf(gnuplotPipe, "set style data histograms\n");

	// Draw a bar chart, specifying the number of columns for each data column purple blue green red
	//fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb '#6FBE47' title 'TP.user', '' using 3 lc rgb '#FF6847' title 'TP.bench', '' using 4 lc rgb '#FDD700' title 'FP', '' using 5 lc rgb '#AB6520' title 'FN'\n", filenamePath.c_str());
	fprintf(gnuplotPipe, "%s\n", plotCommand.c_str());
	//fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'blue' title 'Recall', '' using 3 lc rgb 'green' title 'Precision', '' using 4 lc rgb 'red' title 'F1'\n", outputFileNamePath.c_str());

	// Close the GNUplot pipeline
	pclose(gnuplotPipe);
	AddfileInformation(filenamePath, Info);
	cout << endl << "The classification of result between different data sets is saved in: " << outputFileNamePath << endl;
}
//Generates the segmented interval
vector<string> generateScenarios() {
	vector<string> scenarios;

	if(size_div_vec.empty()){
		scenarios.push_back(">0");
		return scenarios;
	}

	scenarios.push_back("1-" + to_string(size_div_vec[0]) + "bp");

	for(size_t i = 1; i < size_div_vec.size(); ++i){
		scenarios.push_back(to_string(size_div_vec[i - 1] + 1) + "-" + to_string(size_div_vec[i]) + "bp");
	}

	scenarios.push_back(">=" + to_string(size_div_vec.back() + 1) + "bp");

	return scenarios;
}
//Modifying File Information
void AddfileInformation(string &FileNamePath, string &Info){
	ifstream inputFile(FileNamePath);
	if (!inputFile.is_open()) {
		cerr << "Unable to open file" << endl;
		return ;
	}
	// Read the original contents of the file
	vector<string> fileContents;
	string line;
	while (getline(inputFile, line)) {
		fileContents.push_back(line);
	}
	inputFile.close();
	string tempFilename = FileNamePath + ".temp";
	ofstream tempFile(tempFilename);
	if (!tempFile.is_open()) {
		std::cerr << "Unable to create temporary file" << std::endl;
		return ;
	}
	// Writes a new header line
	tempFile << Info << endl;
	// Write the contents of the original file back to the temporary file, keeping the original contents
	for (const std::string& content : fileContents) {
		tempFile << content << std::endl;
	}
	tempFile.close();
	// Delete original file
	if (remove(FileNamePath.c_str()) != 0) {
		std::cerr << "Unable to delete the original file" << std::endl;
		return ;
	}
	// Rename the temporary file to the original file
	if (rename(tempFilename.c_str(), FileNamePath.c_str()) != 0) {
		std::cerr << "Cannot rename temporary files to original files" << std::endl;
		return ;
	}
}

//stat 2
void Histogram_drawing(vector< vector<float> > MeticsValues, string &outputPathname, string &outputBasicMetricschart){
	vector<pair<string, float>> dataV;
	dataV.push_back(make_pair("SeqSim", MeticsValues[0][0]));
	dataV.push_back(make_pair("Recall", MeticsValues[0][1]));
	dataV.push_back(make_pair("Precision", MeticsValues[0][2]));
	dataV.push_back(make_pair("F1-score", MeticsValues[0][3]));

	string outputFileName = "benchmarking_result.png";
	string outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;

	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}
	fprintf(gnuplotPipe, "set term pngcairo size 800,600\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title 'Performance Metric Results Statistics' font 'Times-Roman,16'\n");
	fprintf(gnuplotPipe, "set xlabel 'performance metric' font 'Times-Roman,14'\n");
	fprintf(gnuplotPipe, "set ylabel 'Percentage'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.5\n");
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");
//	fprintf(gnuplotPipe, "plot '-' using 2:xtic(1) title 'Percentage'\n");
	fprintf(gnuplotPipe, "plot '-' using 2:xtic(1) title '%s'\n", method_name.c_str());

	for (const auto& pair : dataV) {
		fprintf(gnuplotPipe, "%s %lf\n", pair.first.c_str(), pair.second);
	}
	fprintf(gnuplotPipe, "e\n");
	pclose(gnuplotPipe);
	cout << "Performance metrics results between different data sets are saved in: " << outputFileNamePath << endl;
}

// stat 2
void Histogram_drawing(vector< vector<int> > MeticsValues1, string &outputPathname, string &outputBasicMetricschart){
		string outputFileName, outputFileNamePath;
		// Prepare the bar chart data
		vector<pair<string, float>> dataV;
		dataV.push_back(make_pair("TP", MeticsValues1[0][0]));
		dataV.push_back(make_pair("FP", MeticsValues1[0][1]));
		dataV.push_back(make_pair("FN", MeticsValues1[0][2]));
		dataV.push_back(make_pair("LP", MeticsValues1[0][3]));

		outputFileName = "result_classification.png";
		outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;
		// Create a pipe connection to GNUplot
		FILE* gnuplotPipe = popen("gnuplot -persist", "w");
		if (!gnuplotPipe) {
			cerr << "Unable to open the GNUplot pipeline" << endl;
			return ;
		}
		// Use the GNUplot command to plot the bar plot and save the output as a PNG file
		fprintf(gnuplotPipe, "set term pngcairo size 800,600\n");
		fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str());
		fprintf(gnuplotPipe, "set title 'Statistics of Benchmark Results' font 'Times-Roman,16'\n");
		fprintf(gnuplotPipe, "set xlabel 'benchmark evaluation metric' font 'Times-Roman,14'\n");
		fprintf(gnuplotPipe, "set ylabel 'Count'\n");
		fprintf(gnuplotPipe, "set style fill solid\n");
		fprintf(gnuplotPipe, "set boxwidth 1.5\n");
		fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
		fprintf(gnuplotPipe, "set yrange [0:*]\n");
		fprintf(gnuplotPipe, "set style data histograms\n");
//	fprintf(gnuplotPipe, "plot '-' using 2:xtic(1) title 'count'\n");
		fprintf(gnuplotPipe, "plot '-' using 2:xtic(1) title '%s'\n", method_name.c_str());
		// Pass the data to GNUplot
		for (const auto& pair : dataV) {
			fprintf(gnuplotPipe, "%s %lf\n", pair.first.c_str(), pair.second);
		}
		fprintf(gnuplotPipe, "e\n");
		pclose(gnuplotPipe);
		cout << endl << "The classification of results between different data sets is saved in " << outputFileNamePath << endl;
}

// stat 1
void SvNumberDistributionGraph(int max_valid_reg_thres, string &refRegSizeFinename, string &refRegSizeFinename_tmp){
	string fileName, fileNamePath, filePathstr_tmp;
	if(refRegSizeFinename.compare("ref_reg_size_benchmark")==0){
		fileName = refRegSizeFinename;
		fileNamePath = refRegSizeFinename_tmp;
		filePathstr_tmp = getContentAfterSlash(fileNamePath) + ".png";
		if(folderPng1.size()==0) folderPng1.push_back(filePathstr_tmp);
	}else if(refRegSizeFinename.compare("ref_reg_size_user")== 0 and max_valid_reg_thres == 0){
		fileName = refRegSizeFinename;
		fileNamePath = refRegSizeFinename_tmp;
		filePathstr_tmp = getContentAfterSlash(fileNamePath) + ".png";
		folderPng1.push_back(filePathstr_tmp);
	}else if(refRegSizeFinename.compare("ref_reg_size_user")== 0 and max_valid_reg_thres > 0){
		fileName = refRegSizeFinename + "_long_filtered";
		fileNamePath = refRegSizeFinename_tmp;
		regSizeFiles.push_back(fileNamePath);
		filePathstr_tmp = getContentAfterSlash(fileNamePath) + ".png";
		folderPng1.push_back(filePathstr_tmp);
	}

	ifstream inputFile(fileNamePath.c_str(), ios::in);
	if (!inputFile) {
		cerr << "Unable to open data file: " << fileNamePath << endl;
		return ;
	}
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");

	if (!gnuplotPipe) {
		cerr << "Unable to start GNUplot." << endl;
		return ;
	}
	// Set the GNUplot option
	fprintf(gnuplotPipe, "set title 'Quantitative distribution'\n");
	fprintf(gnuplotPipe, "set xlabel 'SV size'\n");
	fprintf(gnuplotPipe, "set ylabel 'Count'\n");
	fprintf(gnuplotPipe, "set xrange [-100:100]\n");
	fprintf(gnuplotPipe, "set terminal png  size 800,600\n");
//	fprintf(gnuplotPipe, "set terminal pngcairo size 800,600 dpi 300\n");
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set output '%s.png'\n", fileNamePath.c_str());

	// Draw the graph with the filename variable as the data filename
	fprintf(gnuplotPipe, "plot '%s' every ::2 using 1:2 with linespoints pt 2 linewidth 2 title 'count'\n", fileNamePath.c_str());
	// Close GNUplot
	fflush(gnuplotPipe);
	pclose(gnuplotPipe);

	inputFile.close();
}

//stat 3
void CenterdistanceAndAreasizeratio(string &sizeDifStatDirname){
	string filename1 = "sv_size_ratio_stat_long_filtered", filename2 = "sv_size_dif_stat_long_filtered", filename_pngPath_tmp;
	string title1 = "Variant size ratio", title2 = "Breakpoint distance deviation";
	string fileNamePath1 = sizeDifStatDirname + filename1, fileNamePath2 = sizeDifStatDirname + filename2;
	string filename_png = "result_details.png";
	string filename_pngPath = sizeDifStatDirname + filename_png;
	SVSizeRatioFile.push_back(fileNamePath1);
	filename_pngPath_tmp = getContentAfterSlash(filename_pngPath);
	folderPng3.push_back(filename_pngPath_tmp);
	// Create a pipe connection to the GNUplot process
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline." << endl;
		return ;
	}
	// Read the first data file
	ifstream dataFile1(fileNamePath1);

	if (!dataFile1.is_open()) {
		cerr << "Unable to open the first data file: " << fileNamePath1 << endl;
		return ;
	}

	// Prepare a container for storing data
	vector<string> xLabels1;
	vector<double> yValues1;
	// Read data
	string line1;
	while (getline(dataFile1, line1)) {
		istringstream iss(line1);
		string xLabel;
		double yValue;

		if (iss >> xLabel >> yValue) {
			xLabels1.push_back(xLabel);
			yValues1.push_back(yValue);
		}
	}
	dataFile1.close();

	ifstream dataFile2(fileNamePath2);
	if (!dataFile2.is_open()) {
		cerr << "Unable to open second data file: " << fileNamePath2 << endl;
		return ;
	}

	vector<string> xLabels2;
	vector<double> yValues2;

	string line2;
	while (getline(dataFile2, line2)) {
		istringstream iss(line2);
		string xLabel;
		double yValue;

		if (iss >> xLabel >> yValue) {
			xLabels2.push_back(xLabel);
			yValues2.push_back(yValue);
		}
	}
	dataFile2.close();
	stringstream gnuplotCommands;
	gnuplotCommands << "set terminal pngcairo size 1500,1200\n";
	gnuplotCommands << "set title '" << title1.c_str() << "' font 'Times-Roman,20'\n";
	gnuplotCommands << "set output '" << filename_pngPath.c_str()<< "'\n";
	gnuplotCommands << "set style data histograms\n";
	gnuplotCommands << "set style fill solid border -1\n";
	gnuplotCommands << "set border linecolor rgb 'black' linewidth 2\n";
	gnuplotCommands << "set boxwidth 0.5 relative\n";
	gnuplotCommands << "set xtics nomirror\n";
	gnuplotCommands << "set yrange [0:*]\n";
	gnuplotCommands << "set multiplot layout 2, 1\n";
//	gnuplotCommands << "plot '-' using 2:xtic(1) with boxes lc rgb '#2278B4' title 'count'\n";
	gnuplotCommands << "plot '-' using 2:xtic(1) with boxes lc rgb '#2278B4' title '"<< method_name << "'\n";

	fprintf(gnuplotPipe, "%s", gnuplotCommands.str().c_str());
	fflush(gnuplotPipe);
	// Send the data point of the first data file to GNUplot
	for (size_t i = 0; i < xLabels1.size(); ++i) {
		fprintf(gnuplotPipe, "'%s' %f\n", xLabels1[i].c_str(), yValues1[i]);
	}
	fprintf(gnuplotPipe, "e\n"); // End the first data point input

	// Plot the second histogram in multiplot
	gnuplotCommands.str(""); // Clear previous commands
	gnuplotCommands << "set title '" << title2.c_str() << "' font 'Times-Roman,20'\n";
	gnuplotCommands << "set xtics nomirror rotate by -45\n";
	gnuplotCommands << "set boxwidth 1.0 relative\n"; // Set the column width for the second diagram
	gnuplotCommands << "set border linecolor rgb 'black' linewidth 2\n";
//    gnuplotCommands << "set label 'Center distance difference' at graph 0.3, 0.9 left font 'Time New Roman,18\n"; // Add diagram notes to the second diagram
//	gnuplotCommands << "plot '-' using 2:xtic(1) with boxes lc rgb '#E20613'  title 'count'\n";
	gnuplotCommands << "plot '-' using 2:xtic(1) with boxes lc rgb '#E20613'  title '"<< method_name << "'\n";
	fprintf(gnuplotPipe, "%s", gnuplotCommands.str().c_str());
	fflush(gnuplotPipe);

	// Send the data points of the second data file to GNUplot
	for (size_t i = 0; i < xLabels2.size(); ++i) {
		fprintf(gnuplotPipe, "'%s' %f\n", xLabels2[i].c_str(), yValues2[i]);
	}
	fprintf(gnuplotPipe, "e\n"); // End the second data point input

	pclose(gnuplotPipe);

	cout << "\n" <<"Statistical results of breakpoint distance and variation region size ratio are saved in:" << filename_pngPath << endl;
}

//stat 4
void SVsizeAndNumstatistics(string &sizeNumStatDirname, vector< vector<float> > MeticsValues_4){
	string filename, filenamePath, outputFileName, outputFileNamePath, outputFileNamePath_tmp;
	// Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
//	vector<string> scenarios ={"1-100bp","101-250bp","251-500bp","501-1000bp","1001-2500bp", "2501-5000bp", "5001-10000bp",">10000bp"};
	vector<string> scenarios = generateScenarios();
	// Create a data file to hold the data
	filename = "calculation_result";
	filenamePath = sizeNumStatDirname + filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file: " << filenamePath << endl;
		return ;
	}

	// Write data to a file
	for (size_t i = 0; i < scenarios.size(); i++) {
		dataFile << scenarios[i];
		for (size_t j = 0; j < 4; j++) {
			dataFile << " " << MeticsValues_4[i][j];
		}
		dataFile << endl;
	}
	dataFile.close();

	// Set the file name to save
	outputFileName = "calculation_result.png";
	outputFileNamePath = sizeNumStatDirname + outputFileName;
	outputFileNamePath_tmp = getContentAfterSlash(outputFileNamePath);
	folderPng4.push_back(outputFileNamePath_tmp);
	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}

	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title 'Number statistics of classified metrics in different SV ranges' font 'Times-Roman,14'\n");
	fprintf(gnuplotPipe, "set xlabel 'SV ranges'\n");
	fprintf(gnuplotPipe, "set ylabel 'Percentage'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,10'\n");
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");

	// Draw a bar chart, specifying the number of columns for each data column
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb '#3E4998' title 'SeqSim', '' using 3 lc rgb '#D95F02' title 'Recall', '' using 4 lc rgb '#E7298A' title 'Precision', '' using 5 lc rgb '#5D94A4' title 'F1-score'\n", filenamePath.c_str());
//    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'blue' title 'Recall', '' using 3 lc rgb 'green' title 'Precision', '' using 4 lc rgb 'red' title 'F1'\n", outputFileNamePath.c_str());
	// Close the GNUplot pipeline
	pclose(gnuplotPipe);
	string newInfo = "Range SeqSim Recall Precision F1-score";
	AddfileInformation(filenamePath, newInfo);
	cout << endl << "Statistical results of basic metrics with different SV sizes are saved in: " << outputFileNamePath << endl;
}

//stat 4
void SVsizeAndNumstatistics(string &sizeNumStatDirname, vector< vector<int> > MeticsValues1_4){
	string filename, filenamePath, outputFileName, outputFileNamePath, outputFileNamePath_tmp;
// Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
//	vector<string> scenarios ={"1-100bp","101-250bp","251-500bp","501-1000bp","1001-2500bp", "2501-5000bp", "5001-10000bp",">10000bp"};
	vector<string> scenarios = generateScenarios();
	// Create a data file to hold the data
	filename = "quantitative_statistics";
	filenamePath = sizeNumStatDirname + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file: " << filenamePath << endl;
		return ;
	}

	// Write data to a file
	for (size_t i = 0; i < scenarios.size(); i++) {
		dataFile << scenarios[i];
		for (size_t j = 0; j < 4; j++) {
			dataFile << " " << MeticsValues1_4[i][j];
		}
		dataFile << endl;
	}
	dataFile.close();

	// Set the file name to save
	outputFileName = "quantitative_Statistics.png";
	outputFileNamePath = sizeNumStatDirname + outputFileName;
	outputFileNamePath_tmp=getContentAfterSlash(outputFileNamePath);
	folderPng4.push_back(outputFileNamePath_tmp);
	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}

	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title 'Quantitative Statistics'\n");
	fprintf(gnuplotPipe, "set xlabel 'SV Ranges'\n");
	fprintf(gnuplotPipe, "set ylabel 'Count'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,10'\n");
	fprintf(gnuplotPipe, "set yrange [0:*]\n");
	fprintf(gnuplotPipe, "set style data histograms\n");

	// Draw a bar chart, specifying the number of columns for each data column
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb '#A09DCB' title 'TP.user', '' using 3 lc rgb '#73B8FA' title 'TP.bench', '' using 4 lc rgb '#69BB1E' title 'FP', '' using 5 lc rgb '#EB6365' title 'FN'\n", filenamePath.c_str());
//	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'olive' title 'Precision', '' using 4 lc rgb 'magenta' title 'F1'\n", filenamePath.c_str());

	// Close the GNUplot pipeline
	pclose(gnuplotPipe);
	string newInfo = "Range TP_bench TP_user FP FN";
	AddfileInformation(filenamePath, newInfo);
	cout << "Quantitative statistics of basic metrics with different SV sizes are saved in: " << outputFileNamePath << endl << endl;
}

void ComparisonofMetricsindifferentranges(string &DiffRangeStatDirname, vector< vector<float> > meticsvalues, vector<string> &tool_names, string &FileName, string &outputfilename, vector<string> &sv_files1){
// Prepare the bar chart data, including Recall, Precision, and F1 values for n cases
	vector<string> scenarios;
	if(tool_names.size()>1) scenarios = tool_names;
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			scenarios.push_back(str);
		}
	}
	string title = "Range:" + FileName;

	// Create a data file to hold the data
	string filename = FileName;
	string filenamePath = DiffRangeStatDirname + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file" << endl;
		return ;
	}

	// Write data to a file
	for (size_t i = 0; i < scenarios.size(); i++) {
		dataFile << scenarios[i];
		for (size_t j = 0; j < 3; j++) {
			dataFile << " " << meticsvalues[i][j];
		}
		dataFile << endl;
	}
	dataFile.close();

	// Set the file name to save
	string outputFileName = outputfilename;
	string outputFileNamePath = DiffRangeStatDirname + '/' + outputFileName;

	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}

	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 800,600\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", title.c_str());
	fprintf(gnuplotPipe, "set xlabel 'Tool Name'\n");
	fprintf(gnuplotPipe, "set ylabel 'Value'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
//    fprintf(gnuplotPipe, "set xrange [-0.5:2.5]\n"); // Change the X-axis range
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");

	// Draw a bar chart, specifying the number of columns for each data column
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1'\n", filenamePath.c_str());
//    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'blue' title 'Recall', '' using 3 lc rgb 'green' title 'Precision', '' using 4 lc rgb 'red' title 'F1'\n", outputFileNamePath.c_str());
	// Close the GNUplot pipeline
	pclose(gnuplotPipe);
	string newInfo = "tool Recall Precision F1-score";
	AddfileInformation(filenamePath, newInfo);
	cout << endl << "The comparison results of different SV range indicators are saved in: " << outputFileNamePath << endl;
}

void ComparisonofMetricsindifferentrangesOp(string &outputDiffRangeBasicMetricschart, vector<string> &tool_names, vector<string> &sv_files1){
	string range = "1-100",range_0 = "Range1-100.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_0, tool_names, range, range_0, sv_files1);
	string range1 = "101-250",range_1 = "Range101-250.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_1, tool_names, range1, range_1, sv_files1);
	string range2 = "251-500",range_2 = "Range251-500.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_2, tool_names, range2, range_2, sv_files1);
	string range3 = "501-1000",range_3 = "Range501-1000.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_3, tool_names, range3, range_3, sv_files1);
	string range4 = "1001-5000",range_4 = "Range1001-5000.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_4, tool_names, range4, range_4, sv_files1);
	string range5 = "5001-10000",range_5 = "Range5001-10000.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_5, tool_names, range5, range_5, sv_files1);
	string range6 = "GreaterThan10001",range_6 = "GreaterThan10001.png";
	ComparisonofMetricsindifferentranges(outputDiffRangeBasicMetricschart, MeticsValues4_6, tool_names, range5, range_6, sv_files1);
}

//stat 5
void GenerateFileOp(string &DiffRangeStatDirname, vector<string> &tool_names, vector<string> &filenames, vector<string> &sv_files1){
//	string FileNameStr0 = "1-100";
//	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr0, MeticsValues4_0, filenames, sv_files1);
//	string FileNameStr1 = "101-250";
//	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr1, MeticsValues4_1, filenames, sv_files1);
//	string FileNameStr2 = "251-500";
//	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr2, MeticsValues4_2, filenames, sv_files1);
//	string FileNameStr3 = "501-1000";
//	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr3, MeticsValues4_3, filenames, sv_files1);
//	string FileNameStr4 = "1001-2500";
//	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr4, MeticsValues4_4, filenames, sv_files1);
//	string FileNameStr5 = "2501-5000";
//	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr5, MeticsValues4_5, filenames, sv_files1);
//	string FileNameStr6 = "5001-10000";
//	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr6, MeticsValues4_6, filenames, sv_files1);
//	string FileNameStr7 = "greaterthan10001";
//	GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr7, MeticsValues4_7, filenames, sv_files1);
	for(size_t i = 0; i < size_div_vec.size(); ++i){
		string FileNameStr;
		if (i == 0){
			FileNameStr = "1-" + to_string(size_div_vec[i]);
		} else{
			FileNameStr = to_string(size_div_vec[i - 1] + 1) + "-" + to_string(size_div_vec[i]);
		}

		GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr, MeticsValues4_num[i], filenames, sv_files1);
	}

	if(!size_div_vec.empty()){
		string FileNameStr = "greaterthan" + to_string(size_div_vec[size_div_vec.size()-1] + 1);
		GenerateFile(DiffRangeStatDirname, tool_names, FileNameStr, MeticsValues4_num[size_div_vec.size()], filenames, sv_files1);
	}
}
//stat 5
void GenerateFile(string &DiffRangeStatDirname, vector<string> &tool_names, string &FileName, vector< vector<float> > meticsvalues, vector<string> &filenames, vector<string> &sv_files1){
	vector<string> scenarios;
	vector<string> metrics = {"SeqSim", "Recall", "Precision", "F1-score"};
	if(tool_names.size()>1) scenarios = tool_names;
	//Note to be modified to remove '/'’
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			replaceUnderscoreWithDot(str);
			scenarios.push_back(str);
		}
	}
	string title = "Range:" + FileName;
	// Create a data file to hold the data
	string filename = FileName;
	string filenamePath = DiffRangeStatDirname + '/'+ filename;
	ofstream dataFile(filenamePath);
	if (!dataFile.is_open()) {
		cerr << "Unable to create data file: " << filenamePath << endl;
		return ;
	}
	// Write data to a file
	for (size_t i = 0; i < 4; i++) {
		dataFile << metrics[i];
		for (size_t j = 0; j < scenarios.size(); j++) {
			dataFile << " " << meticsvalues[j][i];
		}
		dataFile << endl;
	}
	dataFile.close();
	filenames.push_back(filenamePath);
}

//NEW stat 5
void GenerateMultiBarCharts(string &outputBasicMetricschart, vector<string>& fileNames, vector<string> &tool_names, vector<string> &sv_files1) {
	string outputFileName, outputFileNamePath, outputFileNamePath_tmp;
	vector<string> scenarios;
	size_t num_rows_plot;
	outputFileName = "different_range.png";
	outputFileNamePath = outputBasicMetricschart + '/' + outputFileName;
	outputFileNamePath_tmp=getContentAfterSlash(outputFileNamePath);
	folderPng5.push_back(outputFileNamePath_tmp);

	if(tool_names.size()>1) scenarios = tool_names;
	else {
		for(size_t i=0; i<sv_files1.size();i++){
			size_t lastSlashPos = sv_files1.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				// If '/' is not found, the entire string is printed
				 str = sv_files1.at(i-1);
			} else {
				// When '/' is reached, the part after '/' is printed
				str = sv_files1.at(i).substr(lastSlashPos + 1);
			}
			replaceUnderscoreWithDot(str);
			scenarios.push_back(str);
		}
	}

//	vector<string> titles =  {
//			        "1-100 bp", "101-250 bp", "251-500 bp",
//			        "501-1000 bp", "1001-2500 bp",  "2501-5000 bp", "5001-10000 bp", ">=10001 bp"
//			    };
	vector<string> titles = generateScenarios();
	num_rows_plot = (size_div_vec.size() + 2) / 2;
	//// Open the Gnuplot pipeline
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to start GNUplot." << endl;
		return;
	}
	// Use the GNUplot command to plot the bar plot and save the output as a PNG file
	//Fig1
	if(num_rows_plot>2){
		fprintf(gnuplotPipe, "set term pngcairo size 1500,1800\n"); // Set the output terminal to PNG
	}else if(num_rows_plot<2){
		fprintf(gnuplotPipe, "set term pngcairo size 1800,1500\n");
	}else{
		fprintf(gnuplotPipe, "set term pngcairo size 1500,1500\n");
	}
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set multiplot layout %zu, 2\n",num_rows_plot);

	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[0].c_str());
	fprintf(gnuplotPipe, "set xlabel 'Metrics'\n");
	fprintf(gnuplotPipe, "set ylabel 'Percentage'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set yrange [0:1]\n");
	fprintf(gnuplotPipe, "set ytics 0.1\n");
	fprintf(gnuplotPipe, "set style data histograms\n");
//	fprintf(gnuplotPipe, "set term pngcairo size 1500,1800\n"); // Set the output terminal to PNG
//	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
//	fprintf(gnuplotPipe, "set multiplot layout 4, 2\n");
//
//	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[0].c_str());
//	fprintf(gnuplotPipe, "set xlabel 'Metrics'\n");
//	fprintf(gnuplotPipe, "set ylabel 'Percentage'\n");
//	fprintf(gnuplotPipe, "set style fill solid\n");
//	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
//	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
//	fprintf(gnuplotPipe, "set yrange [0:1]\n");
//	fprintf(gnuplotPipe, "set ytics 0.1\n");
//	fprintf(gnuplotPipe, "set style data histograms\n");

	// Set key (legend) font size
	fprintf(gnuplotPipe, "set key font 'Times-Roman,10'\n");

	string Info;
	if(tool_names.size()>1)	Info = "Tool";
	else	Info = "callsets";
	for(size_t i = 0; i < size_div_vec.size()+1; i++){
		string plotCommand ="plot '" + fileNames[i] + "'";
		for(size_t j = 0; j < scenarios.size(); j++){
			if(i == 0)
				Info += " " + scenarios[j];
			if(j == 0)
				plotCommand += " using " + to_string(j+2) + ":xtic(1) title '" + scenarios[j] + "'";
			else
				plotCommand += " '' using " + to_string(j+2) + " title '" + scenarios[j] + "'";
			if(j < scenarios.size()-1)
				plotCommand += ",";
		}
		if(i==0)
			fprintf(gnuplotPipe, "%s\n", plotCommand.c_str());
		else{
			fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[i].c_str());
			fprintf(gnuplotPipe, "%s\n", plotCommand.c_str());
		}
	}

/*	// Draw a bar chart, specifying the number of columns for each data column
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1', '' using 5 lc rgb '#BBBB6D' title 'SeqSim'\n", fileNames[0].c_str());
//    fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'blue' title 'Recall', '' using 3 lc rgb 'green' title 'Precision', '' using 4 lc rgb 'red' title 'F1'\n", outputFileNamePath.c_str());

	//Fig2
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[1].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1', '' using 5 lc rgb '#BBBB6D' title 'SeqSim'\n", fileNames[1].c_str());

	//Fig3
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[2].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1', '' using 5 lc rgb '#BBBB6D' title 'SeqSim'\n", fileNames[2].c_str());

	//Fig4
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[3].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1', '' using 5 lc rgb '#BBBB6D' title 'SeqSim'\n", fileNames[3].c_str());

	//Fig5
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[4].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1', '' using 5 lc rgb '#BBBB6D' title 'SeqSim'\n", fileNames[4].c_str());

	//Fig6
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[5].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1', '' using 5 lc rgb '#BBBB6D' title 'SeqSim'\n", fileNames[5].c_str());

	//Fig7
	fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[6].c_str());
	fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title 'Recall', '' using 3 lc rgb 'orange' title 'Precision', '' using 4 lc rgb 'blue' title 'F1', '' using 5 lc rgb '#BBBB6D' title 'SeqSim'\n", fileNames[6].c_str());*/

	// Disable the multi-image mode
	fprintf(gnuplotPipe, "unset multiplot\n");

	// Close Gnuplot
	fflush(gnuplotPipe);
	pclose(gnuplotPipe);
	cout << endl << "The comparison of metrics in different ranges is saved in:" <<outputFileNamePath<< endl;
	for(size_t i=0; i<fileNames.size(); i++){
		AddfileInformation(fileNames[i], Info);
	}
}

//stat 5
//void GenerateSVsizeRatioFileGraph(vector<string> svfilesV, vector<string> toolnameV, string &FileSavingPath){
//	 vector<string> titles;
//	 string filename = "sv_size_ratio.png", filenamePath;
//	 filenamePath = FileSavingPath + '/' + filename;
//	 if(toolnameV.size()>0)	titles = toolnameV;
//	 else {
//	 		for(size_t i=0; i<svfilesV.size();i++){
//	 			size_t lastSlashPos = svfilesV.at(i).find_last_of('/');
//	 			string str;
//	 			if (lastSlashPos == string::npos) {
//	 				// If '/' is not found, the entire string is printed
//	 				 str = svfilesV.at(i-1);
//	 			} else {
//	 				// When '/' is reached, the part after '/' is printed
//	 				str = svfilesV.at(i).substr(lastSlashPos + 1);
//	 			}
//	 			titles.push_back(str);
//	 		}
//	 	}
////	 else titles = svfilesV;
//	 // Create a GNU Plot pipeline
//	 FILE* gnuplotPipe = popen("gnuplot -persist", "w");
//	 if (!gnuplotPipe) {
//		 std::cerr << "Unable to start GNU Plot." << std::endl;
//		 return ;
//	 }
//
//	 fprintf(gnuplotPipe, "set terminal pngcairo size 1200,1000\n");
//	 fprintf(gnuplotPipe, "set output '%s'\n", filenamePath.c_str());
//
//	 // Enter multiplot mode and set the layout
//	 fprintf(gnuplotPipe, "set multiplot layout %d,1\n", (int)SVSizeRatioFile.size());
//
//	 // Loop the subgraph
//	 for (size_t i = 0; i < SVSizeRatioFile.size(); ++i) {
//		 fprintf(gnuplotPipe, "set title '%s' font 'Arial,15'\n", titles[i].c_str());
//		 fprintf(gnuplotPipe, "set style data histograms\n");
//		 fprintf(gnuplotPipe, "set style fill solid border -1\n");
//		 fprintf(gnuplotPipe, "set boxwidth 0.5 relative\n");
//		 fprintf(gnuplotPipe, "set xtics nomirror\n");
//		 fprintf(gnuplotPipe, "set yrange [0:*]\n");
//		 fprintf(gnuplotPipe, "set xlabel 'SV size ratio range'\n");
//		 fprintf(gnuplotPipe, "set ylabel 'Count'\n");
//
//		 // subgraph
//		 fprintf(gnuplotPipe, "plot '-' using 2:xtic(1) with boxes title 'count'\n");
//
//		 ifstream dataFile(SVSizeRatioFile[i]);
//		 if (dataFile.is_open()) {
//			 string line;
//			 while (getline(dataFile, line)) {
//				 fprintf(gnuplotPipe, "%s\n", line.c_str());
//			 }
//			 dataFile.close();
//		 }
//		 fprintf(gnuplotPipe, "e\n"); // The subgraph is drawn
//	 }
//
//	 //Exit the multiplot mode
//	 fprintf(gnuplotPipe, "unset multiplot\n");
//	 // Close the GNU Plot pipeline
//	 pclose(gnuplotPipe);
//
//	 //Print file path
//	 cout << endl << "The results of the sv size ratio identified by different tools are saved in: " << filenamePath << endl;
//}

//Replace above function: stat 5
void GenerateSVsizeRatioFileGraph(vector<string> svfilesV, vector<string> toolnameV, string &FileSavingPath) {
	vector<string> annotations;
	string filename = "sv_size_ratio", filename1 = "sv_quantity_ratio", filenamePath, filename1Path, outputFileNamePath, outputFileName = "sv_size_ratio.png";
	filenamePath = FileSavingPath + '/' + filename;
	filename1Path =  FileSavingPath + '/' + filename1;
	outputFileNamePath = FileSavingPath + '/' +outputFileName;

	if (toolnameV.size() > 0) {
		annotations = toolnameV;
	} else {
		for (size_t i = 0; i < svfilesV.size(); i++) {
			size_t lastSlashPos = svfilesV.at(i).find_last_of('/');
			string str;
			if (lastSlashPos == string::npos) {
				str = svfilesV.at(i - 1);
			} else {
				str = svfilesV.at(i).substr(lastSlashPos + 1);
			}
			replaceUnderscoreWithDot(str);
			annotations.push_back(str);
		}
	}
	vector<string> header_line ={"0.0-0.5","0.5-0.7","0.7-1.2","1.2-2.0","2.0-5.0","5.0-10.0","10.0-50.0","50.0-100.0",">100.0"};
	ofstream outputFile(filenamePath);

	if (!outputFile.is_open()) {
		cerr << "Unable to open file: " << filenamePath << endl;
		return ;
	}
	// Iterate over containers and write their values and header information to a file
	for (size_t i = 0; i < SizeRatioV.size(); ++i) {
		outputFile << header_line[i] << " "; // Add title information

		for (int value : SizeRatioV[i]) {
			outputFile << value << " ";
		}
		outputFile << endl; // Insert newlines between containers
	}
	outputFile.close();
	cout << endl << "The sv size ratio of the text information of different tools is saved to: " << filenamePath << endl;

	ofstream outputFile1(filename1Path);

	if (!outputFile1.is_open()) {
		cerr << "Unable to open file: " << filename1Path << endl;
		return ;
	}
	// Iterate over containers and write their values and header information to a file
	for (size_t i = 0; i < SizeRatio_V.size(); ++i) {
		outputFile1 << header_line[i] << " "; // Add title information

		for (double value : SizeRatio_V[i]) {
			outputFile1 << value << " ";
		}
		outputFile1 << endl; // Insert newlines between containers
	}
	outputFile1.close();
	cout << endl << "The ratio of overlapping SVs of different sizes identified by different tools is saved to: " << filename1Path << endl;

	// Create a pipe connection to GNUplot
	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to open the GNUplot pipeline" << endl;
		return ;
	}
	fprintf(gnuplotPipe, "set term pngcairo\n"); // Set the output terminal to PNG
	fprintf(gnuplotPipe, "set terminal pngcairo size 1200,800\n");
	fprintf(gnuplotPipe, "set output '%s'\n", outputFileNamePath.c_str()); // Set output file name
	fprintf(gnuplotPipe, "set title 'Quantity Statistics of Overlapping SV Size Ratios from Different Tools' font 'Times-Roman,14'\n");
	fprintf(gnuplotPipe, "set xlabel 'Range of size ratio'\n");
	fprintf(gnuplotPipe, "set ylabel 'Count (histogram)'\n");
	fprintf(gnuplotPipe, "set style fill solid\n");
	fprintf(gnuplotPipe, "set boxwidth 1.0\n"); // Set the bar width
	fprintf(gnuplotPipe, "set xtics nomirror font 'Times-Roman,12'\n");
	fprintf(gnuplotPipe, "set ytics nomirror\n");
	fprintf(gnuplotPipe, "set yrange [0:*]\n");

	fprintf(gnuplotPipe, "set y2tics\n");
	fprintf(gnuplotPipe, "set y2label 'Percentage (line)'\n");
	fprintf(gnuplotPipe, "set y2range [0:1]\n");
	fprintf(gnuplotPipe, "set y2tics 0.1\n");

	// Style the bar chart based on the number of input data sets
	if(svfilesV.size()==2){
//		fprintf(gnuplotPipe, "plot '%s' using 2:xtic(1) lc rgb 'purple' title '%s', '' using 3 lc rgb 'red' title '%s' \n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str());
		fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str());
		fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 \n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str());
		cout << endl << "The identification results of the SVs region size ratio of different tools is saved to: " << outputFileNamePath << endl;
	}else if(svfilesV.size()==3){
		fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str());
		fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str());
		cout << endl << "The identification results of the SVs region size ratio of different tools is saved to: " << outputFileNamePath << endl;
	}else if(svfilesV.size()==4){
		fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str());
		fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str());
		cout << endl << "The identification results of the SVs region size ratio of different tools is saved to: " << outputFileNamePath << endl;
	}else if(svfilesV.size()==5){
		fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram, '' using ($6+0.2) lc rgb '#FAB833' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str());
		fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2, '' using 6 lc rgb '#FAB833' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str());
		cout << endl << "The identification results of the SVs region size ratio of different tools is saved to: " << outputFileNamePath << endl;
	}else if(svfilesV.size()==6){
		fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram, '' using ($6+0.2) lc rgb '#FAB833' title '%s' with histogram, '' using ($7+0.2) lc rgb '#227CDE' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str());
		fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2, '' using 6 lc rgb '#FAB833' lw 2 title '%s' with linespoints axes x1y2, '' using 7 lc rgb '#227CDE' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str());
		cout << endl << "The identification results of the SVs region size ratio of different tools is saved to: " << outputFileNamePath << endl;
	}else if(svfilesV.size()==7){
		fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram, '' using ($6+0.2) lc rgb '#FAB833' title '%s' with histogram, '' using ($7+0.2) lc rgb '#227CDE' title '%s' with histogram, '' using ($8+0.2) lc rgb '#55436E' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str(), annotations[6].c_str());
		fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2, '' using 6 lc rgb '#FAB833' lw 2 title '%s' with linespoints axes x1y2, '' using 7 lc rgb '#227CDE' lw 2 title '%s' with linespoints axes x1y2, '' using 8 lc rgb '#55436E' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str(), annotations[6].c_str());
		cout << endl << "The identification results of the SVs region size ratio of different tools is saved to: " << outputFileNamePath << endl;
	}else if(svfilesV.size()==8){
		fprintf(gnuplotPipe, "plot '%s' using ($2-0.2):xtic(1) lc rgb '#E31919' title '%s' with histogram, '' using ($3+0.2):xtic(1) lc rgb '#9BC21B' title '%s' with histogram, '' using ($4+0.2) lc rgb '#B93DB9' title '%s' with histogram, '' using ($5+0.2) lc rgb '#D18D0F' title '%s' with histogram, '' using ($6+0.2) lc rgb '#FAB833' title '%s' with histogram, '' using ($7+0.2) lc rgb '#227CDE' title '%s' with histogram, '' using ($8+0.2) lc rgb '#55436E' title '%s' with histogram, '' using ($9+0.2) lc rgb '#A7C4E1' title '%s' with histogram,\\\n", filenamePath.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str(), annotations[6].c_str(), annotations[7].c_str());
		fprintf(gnuplotPipe, "     '%s' using 2:xtic(1) lc rgb '#E31919' lw 2 title '%s' with linespoints axes x1y2, '' using 3 lc rgb '#9BC21B' lw 2 title '%s' with linespoints axes x1y2 , '' using 4 lc rgb '#B93DB9' lw 2 title '%s' with linespoints axes x1y2, '' using 5 lc rgb '#D18D0F' lw 2 title '%s' with linespoints axes x1y2, '' using 6 lc rgb '#FAB833' lw 2 title '%s' with linespoints axes x1y2, '' using 7 lc rgb '#227CDE' lw 2 title '%s' with linespoints axes x1y2, '' using 8 lc rgb '#55436E' lw 2 title '%s' with linespoints axes x1y2, '' using 9 lc rgb '#A7C4E1' lw 2 title '%s' with linespoints axes x1y2\n", filename1Path.c_str(), annotations[0].c_str(), annotations[1].c_str(), annotations[2].c_str(), annotations[3].c_str(), annotations[4].c_str(), annotations[5].c_str(), annotations[6].c_str(), annotations[7].c_str());
		cout << endl << "The identification results of the SVs region size ratio of different tools is saved to: " << outputFileNamePath << endl;
	}else
		cout << endl << "Up to eight evaluation tools for region size ratio images, exceeding six won't create images." << endl;

	string Add_Info = "#sizeRatio";
	for(size_t i=0; i<annotations.size(); i++){
		Add_Info = Add_Info + " " + annotations[i];
	}
	AddfileInformation(filenamePath, Add_Info);
	AddfileInformation(filename1Path, Add_Info);
}

void RscriptGeneration(string &filePath){
	string RscriptfilePath, Rscriptcode;
//	string filename = "UpSetR.R";
	RscriptfilePath = filePath + "/" + Rscriptfilename;
	Rscriptcode = "library(UpSetR)\n"
				"\n"
				"args <- commandArgs(trailingOnly = TRUE)\n"
				"\n"
				"show_help <- function() {\n"
				"  cat(\"Usage: RScript NewInput.R -i <input_folder> [-n <names>] [-o <output_file>]\\n\")\n"
				"  cat(\"Arguments:\\n\")\n"
				"  cat(\"  -i <input_folder>  : Path to the input folder containing the data files.\\n\")\n"
				"  cat(\"  -n <names>         : Semi-colon separated list of names for the sets.\\n\")\n"
				"  cat(\"  -o <output_file>   : Path to the output PNG file. If not specified, default is 'upset_plot.png' in the input folder.\\n\")\n"
				"  cat(\"  -h                 : Show this help message and exit.\\n\")\n"
				"}\n"
				"if (\"-h\" %in% args) {\n"
				"  show_help()\n"
				"  quit(save = \"no\", status = 0)\n"
				"}\n"
				"\n"
				"if (length(args) == 0 || !(\"-i\" %in% args)) {\n"
				"  stop(\"Use the -i parameter to specify the folder path. Use -h to get help information.\")\n"
				"}\n"
				"\n"
				"input_folder <- args[which(args == \"-i\") + 1]\n"
				"\n"
				"if (!dir.exists(input_folder)) {\n"
				"  stop(\"The specified folder does not exist.\")\n"
				"}\n"
				"\n"
				"file_paths <- list.files(path = input_folder, full.names = TRUE)\n"
				"\n"
				"file_paths <- file_paths[!grepl(\"\\\\.png$\", file_paths)]\n"
				"\n"
				"set_list <- list()\n"
				"\n"
				"for (i in 1:length(file_paths)) {\n"
				"  data <- tryCatch({\n"
				"    read.table(file_paths[i], header = FALSE, comment.char = \"#\")\n"
				"  }, error = function(e) {\n"
				"    stop(paste(\"Error reading file:\", file_paths[i], \"\\nError message:\", e))\n"
				"  })\n"
				"  \n"
				"  if (ncol(data) < 2) {\n"
				"    stop(paste(\"The file does not have at least two columns of data:\", file_paths[i]))\n"
				"  }\n"
				"  \n"
				"  data <- data[, 2]\n"
				"  \n"
				"  set_list[[i]] <- data\n"
				"}\n"
				"\n"
				"if (\"-n\" %in% args) {\n"
				"  names_arg <- args[which(args == \"-n\") + 1]\n"
				"  set_names <- unlist(strsplit(names_arg, \";\"))\n"
				"  \n"
				"  if (length(file_paths) != length(set_names)) {\n"
				"    stop(\"The number of files does not match the number of collection names\")\n"
				"  }\n"
				"  \n"
				"  names(set_list) <- set_names\n"
				"} else {\n"
				"  file_names <- basename(file_paths)\n"
				"  names(set_list) <- file_names\n"
				"}\n"
				"\n"
				"if (\"-o\" %in% args) {\n"
				"  \n"
				"  output_file <- args[which(args == \"-o\") + 1]\n"
				"  \n"
				"  if (dir.exists(output_file)) {\n"
				"    stop(\"The output path should be a file path, not a folder path.\")\n"
				"  }\n"
				"} else {\n"
				"  output_file <- file.path(input_folder, \"upset_plot.png\")\n"
				"}\n"
				"\n"
				"graphics.off()\n"
				"dpi <- 300\n"
				"png(output_file, width = 3200, height = 2000, res = dpi)\n"
				"if(length(set_list)>=6){\n"
				"  size= 2-length(set_list)*(1/5.8)\n"
				"  pointsize=2\n"
				"} else{\n"
				"  size=1.8-(length(set_list)-2)*0.2\n"
				"  \n"
				"  pointsize= 8-length(set_list)*0.5\n"
				"}\n"
				"s=upset(fromList(set_list),\n"
				"      nsets = length(set_list),\n"
				"      order.by = \"freq\",\n"
				"      main.bar.color = \"black\",\n"
				"      sets.bar.color = \"black\",\n"
				"      matrix.color = \"black\",\n"
				"      keep.order = TRUE,\n"
				"      mb.ratio = c(0.7, 0.3), \n"
				"      text.scale = c(size+0.2, 1.5, size+0.1, size, size+0.2, size),\n"
				"      point.size = pointsize\n"
				")\n"
				"s\n";
	ofstream outfile(RscriptfilePath);
	if(outfile.is_open()){
		outfile << Rscriptcode;
		outfile.close();
	}else{
		cerr << "Failed to open the Rscriptfile file, RscriptfilePath: "<< RscriptfilePath << " ." << endl;
	}
}

void TPfilemanipulation(string &filePath){
	int cmdretuern;
	string path, path_tmp, path1, path1_tmp, RscriptCmd;
	string command = "ln -s ";
	for(size_t i=0; i<TPfilesPath.size(); i++){
		path = TPfilesPath[i];
		size_t pos = path.find('/');
		if (pos != std::string::npos && pos + 1 < path.length()) {
			path_tmp = path.substr(pos + 1);
			path_tmp = "../../" + path_tmp;
		}
		command = command + path_tmp + " " + filePath + "/" + SVcallernames[i];
		cmdretuern = system(command.c_str());
		if(cmdretuern!=0){
			cout << "Failed to create symbolic link!" << endl;
		}else{
			command= "ln -s ";
		}
	}
	path1 = filePath;
	size_t pos1 = path1.find('/');
	if (pos1 != std::string::npos && pos1 + 1 < path1.length()) {
		path1_tmp = path1.substr(pos1 + 1);
	}
	RscriptfilenamePath = path1_tmp + "/" + RscriptfilePNGname;
}

void RunRscriptoperation(string &RscriptPath, string &TPfilesPath){
	int runcmdreturn;
	string runcommand, Rscriptpath, TPfilespath, Rscriptpath_tmp;
	Rscriptpath = RscriptPath, TPfilespath = TPfilesPath;
	RscriptfilePNGnamePath = TPfilesPath + "/" + RscriptfilePNGname;
	/*size_t pos = Rscriptpath.find('/');
	if (pos != std::string::npos && pos + 1 < Rscriptpath.length()) {
		Rscriptpath_tmp = Rscriptpath.substr(pos + 1);
	}*/
	Rscriptpath_tmp = Rscriptpath + "/" + Rscriptfilename;
	/*size_t pos1 = TPfilespath.find('/');
	if (pos1 != std::string::npos && pos1 + 1 < TPfilespath.length()) {
		TPfilespath_tmp = TPfilespath.substr(pos1 + 1);
	}*/
	runcommand = "Rscript " + Rscriptpath_tmp + " -i " + TPfilespath + " -o " + RscriptfilePNGnamePath;
	runcmdreturn = system(runcommand.c_str());
	if(runcmdreturn!=0){
		cout << "Rscript running exception causes UpSet diagram drawing failure!" << endl;
	}else{
//		cout << "Rscript of Upset plot runs successfully!" << endl;
		//remove Rscript file
		if(std::remove(Rscriptpath_tmp.c_str()) == 0);
//			cout << "File deleted successfully!" << endl;
		else
			cout << "Error deleting file" << endl;
	}
	cout << endl << "The upset plot generated by TP_bench for multiple user callsets is saved to: " << RscriptfilePNGnamePath << endl;
}

void findCommonFN(string& outputFile, vector<string>& inputFiles) {
	string outputFilename = "bench_low_quality_variant";
	outputFile = outputFile + "/" + outputFilename;
	unordered_map<string, int> lineCounts;
	vector<string> lineOrder;
	int fileCount = inputFiles.size();

	for (const auto& inputFile : inputFiles) {
		ifstream infile(inputFile);
		if (!infile.is_open()) {
			cerr << "Could not open the file " << inputFile << endl;
			continue;
		}
		string line;
		while (getline(infile, line)) {
			// Use the entire line as the key
			if (lineCounts[line] == 0) {
				lineOrder.push_back(line);
			}
			lineCounts[line]++;
		}
		infile.close();
	}
	ofstream outfile(outputFile);
	if (!outfile.is_open()) {
		cerr << "Could not open the output file " << outputFile << endl;
		return;
	}
	for (const auto& key : lineOrder) {
		if ((int)lineCounts[key] == fileCount) {
			outfile << key << endl;
		}
	}
	outfile.close();
}

// Helper function to extract chromosome and position
int chromToInt(const string& chrom) {
	string chromStr = chrom;
	if (chromStr.rfind("chr", 0) == 0) { // Check if "chr" prefix exists
		chromStr = chromStr.substr(3); // Remove "chr" prefix
	}

	if (chromStr == "X") return 23;
	if (chromStr == "Y") return 24;
	if (chromStr == "MT" || chromStr == "M") return 25; // Handle mitochondrial DNA
	try {
		return stoi(chromStr); // Convert numeric chromosome strings directly
	} catch (...) {
		return 26; // Unknown chromosomes will be sorted to the end
	}
}
//pair<string, int> parseChromPos(const string& variantInfo) {
//	size_t firstTab = variantInfo.find('\t');
//	size_t secondTab = variantInfo.find('\t', firstTab + 1);
//
//	// Extract chromosome and position
//	string chrom = variantInfo.substr(0, firstTab); // Chromosome is a string
//	int pos = stoi(variantInfo.substr(firstTab + 1, secondTab - firstTab - 1)); // Position is an integer
//
//	return {chrom, pos};
//}
pair<string, int> parseChromPos(const string& variantInfo){
	size_t firstTab = variantInfo.find('\t');
	size_t secondTab = variantInfo.find('\t', firstTab + 1);

	if(firstTab == string::npos || secondTab == string::npos){
		throw invalid_argument("Invalid variantInfo format");
	}

	// Extract chromosome and position
	string chrom = variantInfo.substr(0, firstTab); // Chromosome is a string
	string posStr = variantInfo.substr(firstTab + 1, secondTab - firstTab - 1); // Position is a string

	// Check if posStr is a valid integer
	for(char c : posStr){
		if(!isdigit(c)){
			throw invalid_argument("Position is not a valid integer");
		}
	}

	int pos = stoi(posStr); // Convert position string to integer

	return {chrom, pos};
}
bool numericSort(const pair<string, int>& lhs, const pair<string, int>& rhs) {
	auto lhsParsed = parseChromPos(lhs.first);
	auto rhsParsed = parseChromPos(rhs.first);

	int lhsChrom = chromToInt(lhsParsed.first);
	int rhsChrom = chromToInt(rhsParsed.first);

	if (lhsChrom == rhsChrom) {
		return lhsParsed.second < rhsParsed.second; // Sort by position if chromosome is the same
	}
	return lhsChrom < rhsChrom; // Otherwise sort by chromosome
}

//Shared FN
void GenerateSharedFNfile(string& outputFile){
	outputSharedFNsfilename = outputFile + "/" + "bench_low_quality_variant.vcf";
	outputBenchfilename = outputFile + "/" + "Refined_benchmark.vcf";

	// Assume benchmarklineMap is of type map<string, int>
	vector<pair<string, int>> sortedVec(benchmarklineMap.begin(), benchmarklineMap.end());
	sort(sortedVec.begin(), sortedVec.end(), numericSort);

	// Writing to the output file
	ofstream outFile(outputSharedFNsfilename);
	if (!outFile) {
		cerr << "Unable to open output file!" << endl;
		return;
	}
	for (const auto& headerLine : benchmarkannotationLines) {
		outFile << headerLine << endl;
	}
	for (const auto& pair : sortedVec) {
		if (pair.second == int(SVcallernames.size())) {
			outFile << pair.first << endl;
		}
	}
	outFile.close();
	cout << endl << "The result of the shared FN records has been output to: " << outputSharedFNsfilename << endl;

	// Writing to the benchmark file
	ofstream outFile1(outputBenchfilename);
	if (!outFile1) {
		cerr << "Unable to open output file!" << endl;
		return;
	}
	for (const auto& headerLine : benchmarkannotationLines) {
		outFile1 << headerLine << endl;
	}
	for (const auto& pair : sortedVec) {
		if (pair.second != int(SVcallernames.size())) {
			outFile1 << pair.first << endl;
		}
	}
	outFile1.close();
	cout << endl << "The refined benchmark set has been saved to: " << outputBenchfilename << endl;
}

//sv distribution
void SvNumberDistributionGraph(const vector<string> &fileNamesPaths, string& svdistributionDirPath) {
	string path;
	SVdistributionPath = svdistributionDirPath + "/" + "sv_distribution.png";
	// Check if each file can be opened
	for (const auto &fileNamePath : fileNamesPaths) {
		ifstream inputFile(fileNamePath.c_str(), ios::in);
		if (!inputFile) {
			cerr << "Unable to open data file: " << fileNamePath << endl;
			return;
		}
		inputFile.close();
	}

	FILE* gnuplotPipe = popen("gnuplot -persist", "w");
	if (!gnuplotPipe) {
		cerr << "Unable to start GNUplot." << endl;
		return;
	}

	// Set the GNUplot options
	fprintf(gnuplotPipe, "set title 'SV distribution'\n");
	fprintf(gnuplotPipe, "set xlabel 'SV size'\n");
	fprintf(gnuplotPipe, "set ylabel 'Count'\n");
	fprintf(gnuplotPipe, "set xrange [-100:100]\n");
	fprintf(gnuplotPipe, "set terminal png size 1000,800\n");
//fprintf(gnuplotPipe, "set terminal pngcairo size 1000,800 dpi 300\n");
	fprintf(gnuplotPipe, "set border linecolor rgb 'black' linewidth 2\n");
	fprintf(gnuplotPipe, "set output '%s'\n", SVdistributionPath.c_str());
//fprintf(gnuplotPipe, "set output 'output_combined.png'\n");

	// Plot each file with its corresponding label
	string plotCommand = "plot ";
	for (size_t i = 0; i < fileNamesPaths.size(); ++i) {
//		plotCommand += "'" + fileNamesPaths[i] + "' every ::2 using 1:2 with linespoints pt 2 linewidth 2 title '" + SVcallernames[i] + "'";
		plotCommand += "'" + fileNamesPaths[i] + "' every 5::2 using 1:2 with lines lw 2 title '" + SVcallernames[i] + "'";
		if (i != fileNamesPaths.size() - 1) {
			plotCommand += ", ";
		}
	}
	plotCommand += "\n";
	fprintf(gnuplotPipe, "%s", plotCommand.c_str());

	// Close GNUplot
	fflush(gnuplotPipe);
	pclose(gnuplotPipe);

	//output
	cout << endl << "The SV distribution statistics of multiple user callsets are saved as: " << SVdistributionPath << endl;
	//get HTML PNGPath
	path = SVdistributionPath;
	size_t pos = path.find('/');
	if (pos != std::string::npos && pos + 1 < path.length()) {
		path = path.substr(pos + 1);
	}
	SVdistributionPath = path;
}
